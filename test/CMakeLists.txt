# Tests unitaires nextDFTB — minimaux, sans framework externe, pilotés par CTest.

# Bibliothèque interne réunissant les modules testables (réutilise les sources
# déjà déclarées dans src/CMakeLists.txt — chemins relatifs au répertoire
# parent).
set(NEXTDFTB_TEST_SOURCES
    ${CMAKE_SOURCE_DIR}/src/utils/kind.f90
    ${CMAKE_SOURCE_DIR}/src/utils/constants.f90
    ${CMAKE_SOURCE_DIR}/src/utils/units.f90
    ${CMAKE_SOURCE_DIR}/src/base/atoms.f90
    ${CMAKE_SOURCE_DIR}/src/base/molecule.f90
    ${CMAKE_SOURCE_DIR}/src/base/neighbor.f90
    ${CMAKE_SOURCE_DIR}/src/base/structure.f90
    ${CMAKE_SOURCE_DIR}/src/errors/errors.f90
    ${CMAKE_SOURCE_DIR}/src/parser/parse_geometry.f90
)

add_library(nextdftb_test_core STATIC ${NEXTDFTB_TEST_SOURCES})
set_target_properties(nextdftb_test_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/test_modules)
target_include_directories(nextdftb_test_core PUBLIC ${CMAKE_BINARY_DIR}/test_modules)

# Module d'assertions partagé.
add_library(nextdftb_test_utils STATIC test_utils.f90)
set_target_properties(nextdftb_test_utils PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/test_modules)
target_include_directories(nextdftb_test_utils PUBLIC ${CMAKE_BINARY_DIR}/test_modules)
target_link_libraries(nextdftb_test_utils PUBLIC nextdftb_test_core)

# Helper : déclare un test et l'enregistre dans CTest.
function(nextdftb_add_test name)
    add_executable(${name} ${name}.f90)
    set_target_properties(${name} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/test_modules)
    target_link_libraries(${name} PRIVATE nextdftb_test_core nextdftb_test_utils)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

nextdftb_add_test(test_units)
nextdftb_add_test(test_atoms)
nextdftb_add_test(test_structure)
nextdftb_add_test(test_parse_geometry)
