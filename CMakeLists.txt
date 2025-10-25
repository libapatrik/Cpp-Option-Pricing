cmake_minimum_required(VERSION 3.10)
project(CppFM)

set(CMAKE_CXX_STANDARD 20)

# Enable testing
enable_testing()

# Add Google Test using FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Main executable
add_executable(CppFM main.cpp
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
)

# Test executables
add_executable(test_volatility_surface_construction tests/test_volatility_surface_construction.cpp
        tests/test_utils.h
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
        tests/test_volatility_surface_construction.cpp
)

add_executable(test_interpolation tests/test_interpolation.cpp
        tests/test_utils.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
)

add_executable(test_dupire tests/test_dupire.cpp
        tests/test_utils.h
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
        tests/test_volatility_surface_construction.cpp
)

add_executable(test_black_scholes tests/test_black_scholes.cpp
        tests/test_utils.h
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
        tests/test_volatility_surface_construction.cpp
)

add_executable(test_greeks tests/test_greeks.cpp
        tests/test_utils.h
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
        tests/test_volatility_surface_construction.cpp
)

add_executable(test_integration tests/test_integration.cpp
        tests/test_utils.h
        Model.cpp
        Model.h
        Utils.h
        Utils.cpp
        PathSimulator.h
        PathSimulator.cpp
        FinancialInstrument.h
        FinancialInstrument.cpp
        Pricer.cpp
        Pricer.h
        DiscountCurve.h
        DiscountCurve.cpp
        VolatilitySurface.cpp
        VolatilitySurface.h
        InterpolationSchemes.cpp
        InterpolationSchemes.h
        BlackScholesFormulas.cpp
        BlackScholesFormulas.h
        tests/test_volatility_surface_construction.cpp
)

# Link test executables with Google Test
target_link_libraries(test_volatility_surface_construction gtest_main)
target_link_libraries(test_interpolation gtest_main)
target_link_libraries(test_dupire gtest_main)
target_link_libraries(test_black_scholes gtest_main)
target_link_libraries(test_greeks gtest_main)
target_link_libraries(test_integration gtest_main)

# Add tests to CTest
add_test(NAME VolatilitySurfaceConstruction COMMAND test_volatility_surface_construction)
add_test(NAME VolatilitySurfaceInterpolation COMMAND test_interpolation)
add_test(NAME VolatilitySurfaceDupire COMMAND test_dupire)
add_test(NAME VolatilitySurfaceBlackScholes COMMAND test_black_scholes)
add_test(NAME VolatilitySurfaceGreeks COMMAND test_greeks)
add_test(NAME VolatilitySurfaceIntegration COMMAND test_integration)
