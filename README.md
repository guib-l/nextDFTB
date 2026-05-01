# nextDFTB

Implémentation minimale et modulaire d'un moteur **DFTB** (Density Functional
Tight-Binding) en **Fortran 2018**. L'objectif est la lisibilité, la cohérence
de l'architecture et un code de calcul fonctionnel sans sur-ingénierie.

Le code lit un fichier `input.dat` décrivant le calcul et un fichier
`geometry.dat` (format xyz étendu), puis exécute le calculateur correspondant
(DFTB ou SKF).

## Architecture

```
input.dat
geometry.dat
   │
   ▼
 main ──► cli ──► parse_input / parse_geometry
                          │
                          ▼
                       driver ──► calculateur (dftb | skf | dft)
                                          │
                                          ▼
                                    write_output
```

### Arborescence

```
src/
├── base/         kind, constants, units, atoms, molecule, structure, globals,
│                 method_basis, method_calc
├── dftb/
│   ├── core/     matel, gamma, density, skrot
│   ├── compute/  dftb_energy, dftb_grad, coulomb, charges
│   ├── solver/   scc, mixer (simple/Broyden), linalg
│   ├── repulsif/ repulsif
│   ├── disp/, pola/        (placeholders)
│   ├── dftbstate.f90       état du dernier calcul
│   ├── write_dftb.f90      sortie spécialisée DFTB
│   └── dftb.f90            façade publique unique
├── skf/          skf, readskf, interp, electronic, repulsive, slakos, write_skf
├── dft/          dft.f90, write_dft.f90 (placeholders)
├── driver/
│   ├── driver.f90, single_point.f90
│   └── opt/, md/           (placeholders)
├── parser/       parse_input, parse_geometry, keywords, property
├── utils/
│   ├── logger/   timer, logger
│   └── output/   output_base, write_output, write_matrix
├── errors/       errors (gestion centralisée)
├── cli.f90       CLI minimaliste
└── main.f90      entrée principale

test/             tests unitaires (CTest, cf. plus bas)
```

### Interface publique du module `dftb`

Le module `dftb.f90` est le **seul** point d'entrée du calculateur. Il expose
exactement :

```fortran
subroutine init(geometry, basis)
subroutine execute(dograd)
function   get_total_energy()
function   get_repulsive_energy()
function   get_coulomb_energy()
function   get_band_energy()
function   get_gradient()           ! (3, natoms)
function   get_charges()            ! (natoms) — Mulliken
```

### Interface publique du module `skf`

```fortran
subroutine init(geometry, basis)
subroutine readslako
subroutine build_repulsion
subroutine build_electronic
function   get_repulsive(atom_A, atom_B, r)
function   get_overlaps(atom_A, atom_B, r, binding)
function   get_hamiltonian(atom_A, atom_B, r, binding)
function   get_hubbard(atom_A)
function   get_eps(atom_A)
```

Aucune autre méthode publique n'est autorisée.

## Prérequis

- CMake ≥ 3.20
- `gfortran` ≥ 11 (ou un compilateur Fortran 2018 équivalent)
- LAPACK (`liblapack-dev` sur Debian/Ubuntu)
- *(optionnel)* OpenMP

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

L'exécutable est produit dans `build/src/nextdftb`.

### Options CMake

| Option                   | Défaut  | Description                            |
|--------------------------|---------|----------------------------------------|
| `CMAKE_BUILD_TYPE`       | Release | `Debug`, `Release`, `RelWithDebInfo`   |
| `NEXTDFTB_ENABLE_OPENMP` | ON      | Active OpenMP dans les noyaux          |
| `NEXTDFTB_ENABLE_TESTS`  | ON      | Compile les tests unitaires `test/`    |

## Utilisation

```bash
build/src/nextdftb input.dat
```

## Format des fichiers d'entrée

Deux fichiers sont nécessaires :

- `input.dat` — description du calcul (sections `GEOMETRY`, `BASIS`, `CALC`,
  `DRIVER`, `OUTPUT`, `SYSTEM`, `OPTION`).
- `geometry.dat` — géométrie au format xyz étendu :

  ```
  natoms
  commentaire
  symbol  x  y  z  charge
  ...
  ```

  Coordonnées en angström, converties en bohr en interne. La colonne
  `charge` est la charge atomique initiale.

### Mots-clés principaux

| Section     | Clé        | Défaut          | Description                             |
|-------------|------------|-----------------|-----------------------------------------|
| `GEOMETRY`  | `GEO`      | `geometry.dat`  | Fichier de géométrie                    |
| `BASIS`     | `SRC`      | `.`             | Dossier des fichiers SKF                |
| `BASIS`     | `EXT`      | `.skf`          | Extension des fichiers SKF              |
| `BASIS`     | `SEP`      | `-`             | Séparateur entre symboles               |
| `BASIS`     | `TYPE`     | `spd`           | `spd`, `spdf` ou `custom`               |
| `BASIS`     | `ORBITALS` | —               | Orbitales de valence par atome          |
| `CALC.DFTB` | `SCC`      | `False`         | Active SCC-DFTB                         |
| `CALC.DFTB` | `MAXSCC`   | `100`           | Itérations SCC max                      |
| `CALC.DFTB` | `TOLSCC`   | `1e-5`          | Tolérance SCC                           |
| `DRIVER`    | `TYPE`     | `SINGLE`        | Type de driver                          |
| `OUTPUT`    | `OUT`      | `results.out`   | Fichier de sortie                       |
| `OUTPUT`    | `LOG`      | (désactivé)     | Fichier de log                          |

### Exemple — `input.dat`

```
GEOMETRY = {
    GEO = geometry.dat
}

BASIS = {
    SRC = ./basis/mio-1-1/
    SEP = -
    TYPE = spd
    ORBITALS = {
        H = 1s
        O = 2s 2p
    }
}

CALC = {
    DFTB = {
        SCC = True
        MAXSCC = 100
        TOLSCC = 1e-7
    }
}

DRIVER = {
    TYPE = SINGLE
}

OUTPUT = {
    OUT = "results.out"
}
```

### Exemple — `geometry.dat`

```
3
water
 O 0.00 0.00 0.00 0.0
 H 0.00 0.96 0.00 0.0
 H 0.93 -0.24 0.00 0.0
```

## Format du fichier de sortie

Le fichier `OUT` (par défaut `results.out`) suit toujours la même structure :

1. **Entête** — bannière, version et horodatage de démarrage.
2. **Calculation properties** — rappel du driver, du calculateur, des
   options SCC, de la géométrie (xyz en Å) et de la base.
3. **Sortie du calculateur** — pour DFTB : entête de la boucle SCF
   (`scc_enabled`, `max_iterations`, `tolerance`), une ligne par itération
   `iter / max|dq|`, statut de convergence.
4. **Résultat final** — pour DFTB : énergies (`E_band`, `E_coulomb`,
   `E_repulsive`, `E_total`) et charges Mulliken par atome.
5. **Execution timings (s)** — temps CPU des phases enregistrées
   (`dftb_init`, `dftb_execute`, `total`, …) via le registre du module
   `timer`.
6. **Footer** — bannière de clôture avec horodatage de fin.

Extrait :

```
======================================================================
nextDFTB — single point
======================================================================
  version             : 0.1.0
  started             : 2026-04-30  07:47:30

----------------------------------------------------------------------
 SCF cycle
----------------------------------------------------------------------
  scc_enabled         : True
  max_iterations      : 100
  tolerance           :     1.0000000000E-07

    iter        max|dq|
       1        7.823478E-01
      ...
      55        7.733238E-08

  converged in        : 55

----------------------------------------------------------------------
 DFTB final result
----------------------------------------------------------------------

  > Energies (Hartree)
  E_total             :    -3.6303661258E+00
  ...

----------------------------------------------------------------------
 Execution timings (s)
----------------------------------------------------------------------
  dftb_init           :       0.026428
  dftb_execute        :       0.014326
  total               :       0.041802
```

Côté code :

- `utils/output/write_output.f90` — entête, sections, helpers `kv_*`,
  `write_timings`, `write_footer`.
- `dftb/write_dftb.f90` — `write_dftb_scc_header`, `write_dftb_scc_iter`,
  `write_dftb_scc_status`, `write_dftb_final`.
- `utils/logger/timer.f90` — `tic`/`toc`/`elapsed` et registre nommé
  (`timer_record`, `timer_count`, `timer_get`, `timer_reset`) consommé par
  `write_timings`.

## Tests unitaires

Les tests sont dans `test/`, sans framework externe (CTest pilote l'exécution).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
(cd build && ctest --output-on-failure)
```

Suites courantes :

| Test                     | Cible                                            |
|--------------------------|--------------------------------------------------|
| `test_units`             | Conversions Bohr↔Å, Ha↔eV (round-trip)           |
| `test_atoms`             | Objet `atoms_t` et `atoms_set`                   |
| `test_structure`         | `structure_init`, `structure_set_atom`, molécule |
| `test_parse_geometry`    | Round-trip `write_geometry` / `read_geometry`    |

## Pipeline de calcul

1. `parse_input` lit `input.dat` ; `parse_geometry` lit `geometry.dat` et
   construit les objets `Atoms` / `Molecule` / `Structure`.
2. Le driver (par défaut `single_point`) appelle l'interface publique du
   calculateur sélectionné.
3. Le calculateur DFTB charge les fichiers SKF nécessaires via le module
   `skf`, construit `H₀`, `S` et la matrice `gamma`, exécute la boucle SCC,
   calcule l'énergie répulsive et assemble l'énergie totale.
4. Le driver écrit l'entête et les propriétés du calcul, le calculateur
   pousse ses étapes (boucle SCF, résultat final) via `write_dftb`, puis
   le driver clôt le fichier avec les timings (`write_timings`) et le
   footer.

## Règles de conception

- Un fichier = un module, une responsabilité par module.
- Pas de logique métier dans `main.f90`.
- Pas de couplage fort : les données circulent explicitement.
- Aucune dépendance externe non listée.
- Erreurs centralisées dans `errors.f90` (pas de `STOP` brutal).

## Licence

[MIT](LICENSE).
