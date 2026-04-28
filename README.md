# nextDFTB

Implémentation minimale et modulaire d'un moteur **DFTB** (Density Functional
Tight-Binding) en **Fortran 2018**. L'objectif est la lisibilité, la cohérence
de l'architecture et un code de calcul fonctionnel sans sur-ingénierie.

Le code lit un input texte décrivant une géométrie, une base SKF et un
calculateur, puis produit un calcul DFTB simple-point (énergie + charges
Mulliken).

## Architecture

```
input.in
   │
   ▼
 main.f90 ──► cli ──► parse_input ──► driver ──► dftb (façade)
                                                     │
                                                     ├─ skf       (lecture/interp SKF)
                                                     ├─ build     (H, S, gamma, density)
                                                     ├─ compute   (diag, charges, coulomb, énergie)
                                                     ├─ repulsif  (énergie répulsive)
                                                     └─ resolved  (SCC, mixer, linalg)
                                                     │
                                                     ▼
                                                 write_output
```

### Arborescence

```
src/
├── base/         constants, units, kind, default, keywords, globals
├── dftb/
│   ├── skf/      readskf, interp
│   ├── build/    dftb_matrix (H/S, Slater-Koster), gamma, density
│   ├── compute/  diag, charges (Mulliken), coulomb, dftb_energy, dftb_grad
│   ├── repulsif/ repulsif
│   ├── resolved/ scc, mixer (simple/Anderson), linalg
│   ├── qmmm/ disp/ pola/   (placeholders)
│   └── dftb.f90  façade publique unique du calculateur
├── skf/, dft/    (placeholders)
├── driver/
│   ├── driver.f90  driver simple-point par défaut
│   └── opt/, md/   (placeholders)
├── io/
│   ├── parser/   parse_input, parse_param
│   ├── logger/   timer, logger
│   └── output/   write_matrix, write_output, dftb_results
├── errors/       errors (gestion centralisée)
├── cli.f90       CLI minimaliste
└── main.f90      entrée principale
```

### Interface publique du module `dftb`

Le module `dftb.f90` est le **seul** point d'entrée du calculateur. Il expose
exactement :

```fortran
subroutine init(geometry, basis, calc)
subroutine execute(dograd)
function   get_total_energy()       ! Hartree
function   get_repulsive_energy()
function   get_coulomb_energy()
function   get_band_energy()
function   get_gradient()           ! (3, natoms)
function   get_charges()            ! (natoms) — Mulliken
```

Aucune autre méthode publique n'est autorisée. Le driver et tout code client
n'utilisent que ces routines — pas d'accès direct aux modules internes.

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

## Utilisation

```bash
build/src/nextdftb input.in
```

Si aucun fichier n'est passé en argument, `input.in` est utilisé par défaut.

## Format d'input

L'input est un fichier texte structuré en sections `SECTION = { ... }`.
Les balises autorisées sont : `GEOMETRY`, `BASIS`, `CALC`, `DRIVER`, `OUTPUT`.

### Mots-clés

| Section    | Clé      | Type     | Défaut         | Description                                |
|------------|----------|----------|----------------|--------------------------------------------|
| `GEOMETRY` | `NATOMS` | int      | —              | Nombre d'atomes (suivi des lignes atomes)  |
| `BASIS`    | `SRC`    | string   | —              | Dossier contenant les fichiers SKF         |
| `BASIS`    | `EXT`    | string   | `.skf`         | Extension des fichiers SKF                 |
| `BASIS`    | `SEP`    | string   | —              | Séparateur entre symboles (ex: `-`)        |
| `BASIS`    | `TYPE`   | string   | `VALENCE`      | Type des paramètres (seul `VALENCE` géré)  |
| `CALC.DFTB`| `SCC`    | bool     | `False`        | Calcul SCC-DFTB                            |
| `CALC.DFTB`| `MAXSCC` | int      | `100`          | Itérations SCC max                         |
| `CALC.DFTB`| `TOLSCC` | float    | `1e-5`         | Tolérance SCC sur max\|Δq\|                |
| `OUTPUT`   | `OUT`    | string   | `results.out`  | Fichier d'output                           |
| `OUTPUT`   | `LOG`    | string/bool | `False`     | Fichier de log (ou `False` pour désactiver) |

Lignes atomes (dans `GEOMETRY`, en angström) :

```
SYMBOLE  X  Y  Z  [GROUPE]
```

Le `GROUPE` est facultatif (entier ≥ 1).

### Exemple

```
GEOMETRY = {
    NATOMS = 3
    O 0.00 0.00 0.00
    H 0.00 0.96 0.00
    H 0.93 -0.24 0.00
}

BASIS = {
    SRC = "./basis/mio-1-1/"
    SEP = "-"
    TYPE = VALENCE
}

CALC = {
    DFTB = {
        SCC = True
        MAXSCC = 100
        TOLSCC = 1e-5
    }
}

OUTPUT = {
    OUT = "results.out"
    LOG = False
}
```

## Pipeline de calcul

1. `parse_input` lit l'input et remplit un `input_t` (géométrie en bohr, base, calc, output).
2. `driver%run_default` initialise `dftb` et appelle `execute`.
3. `dftb%init` charge tous les fichiers SKF nécessaires (toutes les paires
   `A-B.skf` pour les éléments présents) et construit le `basis_system_t`
   (orbitales s/p par atome, mapping atome → orbitales globales).
4. `dftb%execute` :
   - construit `H₀` (Slater-Koster, orbitales s/p) et `S` ;
   - construit la matrice `gamma` (forme Klopman-Ohno) ;
   - exécute la boucle SCC :
     `H = H₀ + ½ S (V_A + V_B)`, diag `H C = S C ε`, occupations Aufbau,
     densité `P`, charges Mulliken, mélange (simple ou Anderson) jusqu'à
     `max|Δq_new − Δq_old| < TOLSCC` ;
   - calcule l'énergie répulsive depuis les splines/polynômes des SKF ;
   - assemble `E_total = E_band − E_coul + E_rep`.
5. `driver` écrit le résumé via `write_output` / `dftb_results`.

## Limites actuelles

- Orbitales `d` non supportées dans la transformation Slater-Koster
  (erreur explicite si rencontrées). Seul s/p est implémenté.
- Gradient analytique non implémenté : `get_gradient()` retourne 0.
- Le mixer Broyden complet n'est pas finalisé ; le mode `broyden` se replie
  sur un mixing simple avec historique.
- `DRIVER`, `DFT`, génération de SKF, QMMM, dispersion, polarisation,
  optimisation et MD : modules placeholders uniquement.
- Format SKF étendu (préfixe `@`) non supporté.

## Règles de conception

- Un fichier = un module, une responsabilité par module.
- Pas de logique métier dans `main.f90`.
- Pas de couplage fort : les données circulent explicitement par arguments.
- Aucune dépendance externe non listée.
- Erreurs centralisées dans `errors.f90` (pas de `STOP` brutal).

## Licence

[MIT](LICENSE).
