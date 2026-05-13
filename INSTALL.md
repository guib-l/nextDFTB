# Installation

Instructions de compilation, d'exécution et de test de **nextDFTB**.

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

Voir le [README](README.md) pour le format des fichiers d'entrée et la
description des sections du fichier `input.dat`.

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
| `test_utils`             | Helpers utilitaires divers                       |
