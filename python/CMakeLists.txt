# =============================================================================
# Python Bindings CMake Configuration
# =============================================================================
#
# This file configures the pybind11 Python module for CppFM.
#
# Usage:
#   From the root CppFM directory:
#     mkdir build_python && cd build_python
#     cmake -DBUILD_PYTHON_BINDINGS=ON ..
#     make -j
#
#   The module will be built as: build_python/python/cppfm.cpython-*.so
#
# To use in Python:
#     import sys
#     sys.path.insert(0, '/path/to/CppFM/build_python/python')
#     import cppfm
#
# Or install with pip (after adding setup.py):
#     pip install -e .
# =============================================================================

# Find pybind11
# Try config mode first (if installed via pip or conda)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    # Fallback: try to find via Python
    message(STATUS "pybind11 not found via CONFIG, trying to find via Python...")

    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE PYBIND11_FIND_RESULT
    )

    if(PYBIND11_FIND_RESULT EQUAL 0)
        list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        # Last resort: FetchContent
        message(STATUS "pybind11 not found, fetching from GitHub...")
        include(FetchContent)
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG v2.11.1
        )
        FetchContent_MakeAvailable(pybind11)
    endif()
endif()

message(STATUS "pybind11 found: ${pybind11_VERSION}")

# =============================================================================
# Python Module Target
# =============================================================================

pybind11_add_module(_core
    ${CMAKE_SOURCE_DIR}/src/bindings/module_main.cpp
    ${CMAKE_SOURCE_DIR}/src/bindings/bind_market.cpp
    ${CMAKE_SOURCE_DIR}/src/bindings/bind_models.cpp
    ${CMAKE_SOURCE_DIR}/src/bindings/bind_simulators.cpp
    ${CMAKE_SOURCE_DIR}/src/bindings/bind_pricers.cpp
)

# Link against the core CppFM library
target_link_libraries(_core PRIVATE CppFM_lib)

# Include directories
target_include_directories(_core PRIVATE
    ${CMAKE_SOURCE_DIR}/src/bindings
)

# C++20 standard
target_compile_features(_core PRIVATE cxx_std_20)

# Compiler-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(_core PRIVATE
        -O3
        -fvisibility=hidden
    )
endif()

# Output to python directory for easy import
set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# =============================================================================
# Installation (for scikit-build-core)
# =============================================================================

install(TARGETS _core LIBRARY DESTINATION .)

# Print helpful message
message(STATUS "")
message(STATUS "=== Python Bindings Configuration ===")
message(STATUS "Module will be built at: ${CMAKE_CURRENT_BINARY_DIR}/_core*.so")
message(STATUS "Install with: pip install .")
message(STATUS "======================================")
