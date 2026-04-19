# nextDFTB — Python bindings (pybind11).
#
# Produces the compiled extension module `nextdftb._core` and stages the
# pure-Python package files next to it in ${CMAKE_BINARY_DIR}/python, so
# the in-tree layout works with `PYTHONPATH=build/python python`.

if(NOT NEXTDFTB_ENABLE_PYTHON)
    return()
endif()

find_package(pybind11 REQUIRED CONFIG)

pybind11_add_module(nextdftb_py_core MODULE python/module.cpp)
add_library(nextdftb::py_core ALIAS nextdftb_py_core)

target_link_libraries(nextdftb_py_core PRIVATE nextdftb::core_cpp)

set_target_properties(nextdftb_py_core PROPERTIES
    OUTPUT_NAME              _core
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/nextdftb"
)

# Under scikit-build-core the extension lands in <site-packages>/nextdftb/
# and libnextdftb_core_cpp.so lands in <site-packages>/nextdftb/lib/ — make
# sure the runtime loader finds it relative to the extension itself.
if(SKBUILD)
    set_target_properties(nextdftb_py_core PROPERTIES
        INSTALL_RPATH "$ORIGIN/lib"
        BUILD_WITH_INSTALL_RPATH ON
    )
endif()

# Stage Python sources.
set(_nextdftb_py_pkg_dir "${CMAKE_BINARY_DIR}/python/nextdftb")
file(MAKE_DIRECTORY "${_nextdftb_py_pkg_dir}")

configure_file(python/nextdftb/__init__.py     "${_nextdftb_py_pkg_dir}/__init__.py"     COPYONLY)
configure_file(python/nextdftb/_log_handler.py "${_nextdftb_py_pkg_dir}/_log_handler.py" COPYONLY)

# Install layout.
#   * Under scikit-build-core (`pip install .`), SKBUILD is defined and
#     CMAKE_INSTALL_PREFIX already points at the wheel's platlib + the
#     `wheel.install-dir = "nextdftb"` from pyproject.toml, so we install
#     relative to that.
#   * Under a plain `cmake --install`, we fall back to Python3_SITEARCH.
if(SKBUILD)
    install(TARGETS nextdftb_py_core LIBRARY DESTINATION ".")
    install(DIRECTORY python/nextdftb/
            DESTINATION   "."
            FILES_MATCHING PATTERN "*.py")
elseif(Python3_SITEARCH)
    install(TARGETS nextdftb_py_core
            LIBRARY DESTINATION "${Python3_SITEARCH}/nextdftb")
    install(DIRECTORY python/nextdftb/
            DESTINATION   "${Python3_SITEARCH}/nextdftb"
            FILES_MATCHING PATTERN "*.py")
endif()
