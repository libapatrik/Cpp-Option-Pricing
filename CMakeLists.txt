cmake_minimum_required(VERSION 3.18)

# Set CMake policies IMMEDIATELY after cmake_minimum_required   // must be set here! Order matters
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()


project(CppFM)

set(CMAKE_CXX_STANDARD 20)

# Required for building shared libraries (e.g., Python bindings) that link static libs
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ============================================================================
# Dependencies
# ============================================================================

# Google Test via FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Fetch and build googletest
FetchContent_MakeAvailable(googletest)

# Make Google Test headers available for IDE indexing
# This helps CLion find the headers for code completion and error checking
include_directories(${googletest_SOURCE_DIR}/googletest/include)
include_directories(${googletest_SOURCE_DIR}/googlemock/include)

# Boost
set(Boost_ROOT "/opt/homebrew/opt/boost")
find_package(Boost REQUIRED)

# OpenMP - for parallel Monte Carlo simulation
# macOS Apple Clang needs explicit configuration for Homebrew's libomp
if(APPLE)
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${OpenMP_ROOT}/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OpenMP_ROOT}/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${OpenMP_ROOT}/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

# TBB - for C++17 parallel algorithms (std::execution::par)
# Install with: brew install tbb
find_package(TBB REQUIRED)

# Enable testing
enable_testing()

# ============================================================================
# Core Library
# ============================================================================
add_library(CppFM_lib STATIC
    # Core Models
    Model.cpp
    Model.h
    Utils.cpp
    Utils.h
    
    # Path Simulators
    PathSimulator.cpp
    PathSimulator.h
    PathSimulator2D.cpp
    PathSimulator2D.h
    
    # Financial Instruments
    FinancialInstrument.cpp
    FinancialInstrument.h
    Pricer.cpp
    Pricer.h
    
    # Market Data
    DiscountCurve.cpp
    DiscountCurve.h
    VolatilitySurface.cpp
    VolatilitySurface.h
    InterpolationSchemes.cpp
    InterpolationSchemes.h
    BlackScholesFormulas.cpp
    BlackScholesFormulas.h
    
    # PDE Solver
    PDEs/PDE.cpp
    PDEs/PDE.h
    PDEs/Grid.cpp
    PDEs/Grid.h
    PDEs/BoundaryConditions.cpp
    PDEs/BoundaryConditions.h
    PDEs/Solver.cpp
    PDEs/Solver.h
    PDEs/PDEGreeksCalculator.cpp
    PDEs/PDEGreeksCalculator.h
)

# Link Boost to the library
# Use PUBLIC SYSTEM to:
#   1. Propagate Boost headers to targets linking CppFM_lib
#   2. Suppress warnings from Boost headers (-isystem)
target_include_directories(CppFM_lib SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})

# Link OpenMP to the library
target_link_libraries(CppFM_lib PUBLIC OpenMP::OpenMP_CXX)

# Link TBB for parallel algorithms
target_link_libraries(CppFM_lib PUBLIC TBB::tbb)

# ============================================================================
# Main Executable
# ============================================================================
add_executable(CppFM main.cpp)
target_link_libraries(CppFM CppFM_lib)

# ============================================================================
# Test Executables
# ============================================================================

# Volatility Surface Tests
add_executable(test_volatility_surface_construction tests/test_volatility_surface_construction.cpp tests/test_utils.h)
target_link_libraries(test_volatility_surface_construction CppFM_lib gtest_main)

add_executable(test_interpolation tests/test_interpolation.cpp tests/test_utils.h)
target_link_libraries(test_interpolation CppFM_lib gtest_main)

add_executable(test_dupire tests/test_dupire.cpp tests/test_utils.h)
target_link_libraries(test_dupire CppFM_lib gtest_main)

add_executable(test_dupire_pricer tests/test_dupire_pricer.cpp tests/test_utils.h)
target_link_libraries(test_dupire_pricer CppFM_lib gtest_main)

add_executable(test_black_scholes tests/test_black_scholes.cpp tests/test_utils.h)
target_link_libraries(test_black_scholes CppFM_lib gtest_main)

add_executable(test_greeks tests/test_greeks.cpp tests/test_utils.h)
target_link_libraries(test_greeks CppFM_lib gtest_main)

add_executable(test_integration tests/test_integration.cpp tests/test_utils.h)
target_link_libraries(test_integration CppFM_lib gtest_main)

# PDE Solver Tests
add_executable(test_pde_solver tests/test_pde_solver.cpp)
target_link_libraries(test_pde_solver CppFM_lib gtest_main)

add_executable(test_pde_stability tests/test_pde_stability.cpp)
target_link_libraries(test_pde_stability CppFM_lib gtest_main)

add_executable(test_tridiagonal_solver tests/test_tridiagonal_solver.cpp)
target_link_libraries(test_tridiagonal_solver CppFM_lib gtest_main)

add_executable(test_pde_greeks tests/test_pde_greeks.cpp)
target_link_libraries(test_pde_greeks CppFM_lib gtest_main)

add_executable(test_cos_method tests/test_cos_method.cpp)
target_link_libraries(test_cos_method CppFM_lib gtest_main)

add_executable(test_newton_methods tests/test_newton_methods.cpp)
target_link_libraries(test_newton_methods CppFM_lib gtest_main)

add_executable(test_heston_model tests/test_heston_model.cpp)
target_link_libraries(test_heston_model CppFM_lib gtest_main)

add_executable(test_heston_scheme tests/test_heston_scheme.cpp)
target_link_libraries(test_heston_scheme CppFM_lib gtest_main)

add_executable(test_heston_slv tests/test_heston_slv.cpp)
target_link_libraries(test_heston_slv CppFM_lib gtest_main)

add_executable(test_pde_grid tests/test_pde_grid.cpp)
target_link_libraries(test_pde_grid CppFM_lib gtest_main)

add_executable(test_bessel_function tests/test_bessel_function.cpp)
target_link_libraries(test_bessel_function CppFM_lib gtest_main)

# ============================================================================
# CTest Configuration
# ============================================================================
add_test(NAME VolatilitySurfaceConstruction COMMAND test_volatility_surface_construction)
add_test(NAME VolatilitySurfaceInterpolation COMMAND test_interpolation)
add_test(NAME VolatilitySurfaceDupire COMMAND test_dupire)
add_test(NAME DupirePricer COMMAND test_dupire_pricer)
add_test(NAME VolatilitySurfaceBlackScholes COMMAND test_black_scholes)
add_test(NAME VolatilitySurfaceGreeks COMMAND test_greeks)
add_test(NAME VolatilitySurfaceIntegration COMMAND test_integration)
add_test(NAME PDESolver COMMAND test_pde_solver)
add_test(NAME PDEStability COMMAND test_pde_stability)
add_test(NAME TridiagonalSolver COMMAND test_tridiagonal_solver)
add_test(NAME PDEGreeks COMMAND test_pde_greeks)
add_test(NAME COSMethod COMMAND test_cos_method)
add_test(NAME NewtonMethods COMMAND test_newton_methods)
add_test(NAME HestonModel COMMAND test_heston_model)
add_test(NAME HestonScheme COMMAND test_heston_scheme)
add_test(NAME HestonSLV COMMAND test_heston_slv)
add_test(NAME PDEGrid COMMAND test_pde_grid)
add_test(NAME BesselFunction COMMAND test_bessel_function)

# ============================================================================
# Python Bindings - to use in QFB
# ============================================================================
# Enable with: cmake -DBUILD_PYTHON_BINDINGS=ON ..

option(BUILD_PYTHON_BINDINGS "Build Python bindings using pybind11" OFF)

if(BUILD_PYTHON_BINDINGS)
    message(STATUS "Python bindings enabled")

    # Find Python interpreter and development libraries
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    message(STATUS "Python found: ${Python_EXECUTABLE} (${Python_VERSION})")

    # Add the python subdirectory
    add_subdirectory(python)
endif()