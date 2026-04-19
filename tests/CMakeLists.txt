# nextDFTB — test suite.
#
# CTest drives pytest-based runners. Each runner is a small Python script
# that exercises one slice of the stack.

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(NEXTDFTB_TEST_RUNNER_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/runners)
set(NEXTDFTB_TEST_REFERENCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/reference_data)

# ---- CLI smoke test: nextdftb runs the infra test -------------------
if(TARGET nextdftb_cli)
    add_test(
        NAME cli_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${NEXTDFTB_TEST_RUNNER_DIR}/test_cli_smoke.py
                --binary $<TARGET_FILE:nextdftb_cli>
    )
    set_tests_properties(cli_smoke PROPERTIES LABELS "smoke;cli")
endif()

# ---- Python import + infra test round-trip --------------------------
if(NEXTDFTB_ENABLE_PYTHON)
    add_test(
        NAME python_import
        COMMAND ${Python3_EXECUTABLE}
                ${NEXTDFTB_TEST_RUNNER_DIR}/test_python_import.py
    )
    set_tests_properties(python_import PROPERTIES
        LABELS      "smoke;python"
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
    )
endif()
