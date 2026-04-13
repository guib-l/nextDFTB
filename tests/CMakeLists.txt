# Tests are driven externally via pytest; CMake only registers CTest entries
# that invoke the runner scripts. No .f90 test sources live in this tree.

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(NEXTDFTB_TEST_RUNNER_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/runners)
set(NEXTDFTB_TEST_REFERENCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/reference_data)

# Example registration (enable once the binary and a reference case exist):
#
# add_test(
#     NAME h2_scc
#     COMMAND ${Python3_EXECUTABLE} -m pytest
#             ${NEXTDFTB_TEST_RUNNER_DIR}/test_h2_scc.py
#             --binary=$<TARGET_FILE:nextdftb>
#             --reference=${NEXTDFTB_TEST_REFERENCE_DIR}/h2
# )
# set_tests_properties(h2_scc PROPERTIES LABELS "integration")
