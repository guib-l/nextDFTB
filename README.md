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
├── base/         atoms, molecule, structure, neighbor
├── dftb/
│   ├── core/     matel, gamma, density, skrot
│   ├── compute/  dftb_energy, dftb_grad, coulomb, charges
│   ├── solver/
│   │   ├── scc.f90, linalg.f90
│   │   └── mixer/  mixer (interface), simple, random, broyden, factory
│   ├── repulsif/ repulsif
│   ├── disp/, pola/        (placeholders)
│   ├── orbitals.f90        construction de la base d'orbitales
│   ├── dftbstate.f90       état du dernier calcul
│   ├── write_dftb.f90      sortie spécialisée DFTB
│   └── dftb.f90            façade publique (type dftb_calc_t)
├── skf/          skf, readskf, polint, interp, electronic, repulsive,
│                 slakos, write_skf
├── dft/          dft.f90, write_dft.f90 (placeholders)
├── driver/
│   ├── abstract_calc.f90, abstract_basis.f90  interfaces abstraites
│   ├── driver.f90, single_point.f90
│   └── opt/, md/           (placeholders)
├── parser/       parse_input, parse_geometry, keywords, property
├── utils/        kind, constants, units
│   ├── logger/   timer, logger
│   └── output/   output_base, write_output, write_matrix
├── errors/       errors (gestion centralisée)
├── globals.f90   variables globales partagées
├── cli.f90       CLI minimaliste
└── main.f90      entrée principale

test/             tests unitaires (CTest, cf. plus bas)
```

### Interface publique du module `dftb`

Le module `dftb.f90` expose le type `dftb_calc_t`, qui étend l'interface
abstraite `method_calc_t`. Le calculateur est manipulé exclusivement via
ce type :

```fortran
type, extends(method_calc_t) :: dftb_calc_t
contains
    procedure :: init               ! (struct, basis, method)
    procedure :: build              ! charge SKF, alloue état
    procedure :: execute            ! SCF + énergies
    procedure :: get_total_energy
    procedure :: get_total_gradient ! (3, natoms)
    procedure :: write_output
    ! méthodes spécifiques DFTB
    procedure :: get_band_energy
    procedure :: get_coulomb_energy
    procedure :: get_repulsive_energy
    procedure :: get_charges        ! Mulliken (natoms)
end type dftb_calc_t
```

### Interface publique du module `skf`

```fortran
subroutine init(geometry, basis)
subroutine readslako
subroutine build_repulsion
subroutine build_electronic

! Évaluation Slater-Koster (H, S et leurs dérivées par rapport à r)
function get_hamiltonian (atom_A, atom_B, r, binding)
function get_overlaps    (atom_A, atom_B, r, binding)
function get_dhamiltonian(atom_A, atom_B, r, binding)
function get_doverlaps   (atom_A, atom_B, r, binding)

! Énergie répulsive
function get_repulsive       (atom_A, atom_B, r)   ! alias spline
function get_repulsion_spline(atom_A, atom_B, r)
function get_repulsion_poly  (atom_A, atom_B, r)
function get_drepulsive      (atom_A, atom_B, r)

! Paramètres atomiques
function get_hubbard    (atom_A)
function get_eps        (atom_A)
function get_occupations(atom_A)
function get_mass       (atom_A)

! Index / introspection
function nelements()
function element_symbol(i)
function element_index (sym)
subroutine write_summary(...)
```

L'argument `binding` est le libellé d'une intégrale SK ; sans cet argument,
l'API renvoie le vecteur des 10 intégrales (convention DFTB+ :
`1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ`).

L'interpolation H/S utilise `polint` (Neville, ordre fixe) sur une fenêtre
centrée autour de r. Les dérivées analytiques sont fournies pour la
répulsion ; celles de H/S sont obtenues via `dpolint` sur la même fenêtre.

## Installation

Voir [INSTALL.md](INSTALL.md) pour les prérequis, la compilation, les options
CMake et l'exécution des tests.

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

| Section            | Clé            | Défaut         | Description                                   |
|--------------------|----------------|----------------|-----------------------------------------------|
| `GEOMETRY`         | `GEO`          | `geometry.dat` | Fichier de géométrie                          |
| `BASIS`            | `SRC`          | `.`            | Dossier des fichiers SKF                      |
| `BASIS`            | `EXT`          | `.skf`         | Extension des fichiers SKF                    |
| `BASIS`            | `SEP`          | `-`            | Séparateur entre symboles                     |
| `BASIS`            | `TYPE`         | `spd`          | `spd`, `spdf` ou `custom`                     |
| `BASIS`            | `ORBITALS`     | —              | Orbitales de valence par atome                |
| `CALC.DFTB`        | `WRITE_MATRIX` | `False`        | Sauvegarde H, S, γ sur disque                 |
| `CALC.DFTB`        | `GAMMA`        | `MEAN`         | `BASE` \| `MEAN` \| `STDR` (Standard DFTB+)   |
| `CALC.DFTB.SCHEM`  | `TYPE`         | `SCC`          | `BASIC` \| `NOSCC` \| `SCC`                   |
| `CALC.DFTB.SCHEM`  | `MAXSCC`       | `100`          | Itérations SCC max                            |
| `CALC.DFTB.SCHEM`  | `TOLSCC`       | `1e-5`         | Tolérance SCC                                 |
| `CALC.DFTB.MIXING` | `TYPE`         | `SIMPLE`       | `SIMPLE` \| `RANDOM` \| `BROYDEN`             |
| `CALC.DFTB.MIXING` | `FACTOR`       | `0.1`          | Facteur de mélange                            |
| `CALC.DFTB.MIXING` | `HISTORY`      | `6`            | Profondeur d'historique (Broyden)             |
| `CALC.DFTB.MIXING` | `OMEGA0`       | `0.01`         | Régularisation Broyden                        |
| `DRIVER`           | `TYPE`         | `SINGLE`       | Type de driver                                |
| `OUTPUT`           | `OUT`          | `results.out`  | Fichier de sortie                             |
| `OUTPUT`           | `LOG`          | (désactivé)    | Fichier de log                                |

### Schéma γ (mot-clé `GAMMA`)

| Valeur | Description                                                            |
|--------|------------------------------------------------------------------------|
| `BASE` | γ_AB = 1/r (Coulomb pur, sans correction de Hubbard)                   |
| `MEAN` | γ_AB analytique avec un Hubbard moyen `(U_a + U_b)/2`                  |
| `STDR` | γ_AB DFTB+ standard (Slater, expressions distinctes pour U_a ≠ U_b)    |

### Mixers SCF

- `SIMPLE` — α p_new + (1-α) p_old (paramètre `FACTOR`)
- `RANDOM` — perturbation aléatoire (debug)
- `BROYDEN` — Broyden modifié (paramètres `HISTORY`, `OMEGA0`, `FACTOR`)

### Exemple — `input.dat`

```
GEOMETRY = {
    GEO = geometry.dat
}

BASIS = {
    SRC = ./mio-1-1/
    SEP = -
    TYPE = spd
    ORBITALS = {
        H = 1s
        O = 2s 2p
    }
}

CALC = {
    DFTB = {
        WRITE_MATRIX = False
        GAMMA = STDR
        SCHEM = {
            TYPE   = SCC
            MAXSCC = 100
            TOLSCC = 1e-7
        }
        MIXING = {
            TYPE    = BROYDEN
            FACTOR  = 0.1
            HISTORY = 6
            OMEGA0  = 0.01
        }
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
2. **Calculation properties** — rappel du driver, du calculateur, de la
   géométrie (xyz en Å) et de la base.
3. **SCF cycle** — entête (`scc_enabled`, `max_iterations`, `tolerance`),
   une ligne par itération `iter / dE_elec / max|dq|`, statut de
   convergence.
4. **Résultat final** — pour DFTB : énergies (`E_band`, `E_coulomb`,
   `E_repulsive`, `E_total`) et charges Mulliken par atome.
5. **TIMER** — tableau récapitulatif (boucle, temps, % du total) issu du
   registre du module `timer`.
6. **Footer** — bannière de clôture avec horodatage de fin.

Extrait :

```
================================================================================
nextDFTB — Calculation
================================================================================
  version             : 0.1.0
  started             : 2026-05-05  19:38:48

--------------------------------------------------------------------------------
 SCF cycle
--------------------------------------------------------------------------------
  scc_enabled         : True
  max_iterations      : 100
  tolerance           :     1.0000000000E-07

    iter         dE_elec         max|dq|
  --------  ----------------  ----------------
       1        ---                 8.350189E-01
       2        2.770352E-01        5.814383E-01
      ...
      11        1.372030E-07        1.964562E-09

  converged in        : 11

--------------------------------------------------------------------------------
 DFTB final result
--------------------------------------------------------------------------------

  > Energies (Hartree)
  E_band              :    -6.9195224461E+00
  E_coulomb           :     3.4247076625E-02
  E_repulsive         :     2.6842198734E-02
  E_total             :    -6.9269273240E+00

-------------------------------------------------------------------------------
 TIMER ::
 ┌───────────────────────────────────────────┬────────┬───────────┬───────────┐
 |           Name                            |  Loop  |  Time [s] | Relat [%] |
 ├───────────────────────────────────────────┼────────┼───────────┼───────────┤
 | CALC_INIT_BUILD                           |       1|     0.0043|     85.937|
 | CALC_EXECUTE                              |       1|     0.0005|      9.735|
 ├───────────────────────────────────────────┼────────┼───────────┼───────────┤
 | TOTAL                                     |       1|     0.0050|    100.000|
 └───────────────────────────────────────────┴────────┴───────────┴───────────┘
```

Côté code :

- `utils/output/write_output.f90` — entête, sections, helpers `kv_*`,
  `write_footer`.
- `dftb/write_dftb.f90` — `write_dftb_scc_header`, `write_dftb_scc_iter`,
  `write_dftb_scc_status`, `write_dftb_final`.
- `utils/logger/timer.f90` — `tic`/`toc`/`elapsed` et registre nommé
  (`timer_record`, `timer_count`, `timer_get`, `timer_reset`) consommé par
  l'écriture du tableau TIMER.

## Pipeline de calcul

1. `parse_input` lit `input.dat` ; `parse_geometry` lit `geometry.dat` et
   construit les objets `Atoms` / `Molecule` / `Structure`.
2. Le driver (par défaut `single_point`) instancie le calculateur sélectionné
   (ici `dftb_calc_t`) puis appelle successivement `init`, `build`,
   `execute` et `write_output`.
3. Le calculateur DFTB charge les fichiers SKF nécessaires via le module
   `skf`, construit `H₀`, `S` (jusqu'aux orbitales d) et la matrice γ
   (schéma `BASE`, `MEAN` ou `STDR`), exécute la boucle SCC avec le mixer
   choisi, calcule l'énergie répulsive et assemble l'énergie totale.
4. Le driver écrit l'entête et les propriétés du calcul, le calculateur
   pousse ses étapes (boucle SCF, résultat final) via `write_dftb`, puis
   le driver clôt le fichier avec le tableau TIMER et le footer.

## Règles de conception

- Un fichier = un module, une responsabilité par module.
- Pas de logique métier dans `main.f90`.
- Pas de couplage fort : les données circulent explicitement.
- Aucune dépendance externe non listée.
- Erreurs centralisées dans `errors.f90` (pas de `STOP` brutal).

## Licence

[MIT](LICENSE).
