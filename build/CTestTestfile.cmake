# CMake generated Testfile for 
# Source directory: /Users/patrikliba/CLionProjects/CppFM
# Build directory: /Users/patrikliba/CLionProjects/CppFM/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(VolatilitySurfaceConstruction "/Users/patrikliba/CLionProjects/CppFM/build/test_volatility_surface_construction")
set_tests_properties(VolatilitySurfaceConstruction PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;238;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(VolatilitySurfaceInterpolation "/Users/patrikliba/CLionProjects/CppFM/build/test_interpolation")
set_tests_properties(VolatilitySurfaceInterpolation PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;239;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(VolatilitySurfaceDupire "/Users/patrikliba/CLionProjects/CppFM/build/test_dupire")
set_tests_properties(VolatilitySurfaceDupire PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;240;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(VolatilitySurfaceBlackScholes "/Users/patrikliba/CLionProjects/CppFM/build/test_black_scholes")
set_tests_properties(VolatilitySurfaceBlackScholes PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;241;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(VolatilitySurfaceGreeks "/Users/patrikliba/CLionProjects/CppFM/build/test_greeks")
set_tests_properties(VolatilitySurfaceGreeks PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;242;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(VolatilitySurfaceIntegration "/Users/patrikliba/CLionProjects/CppFM/build/test_integration")
set_tests_properties(VolatilitySurfaceIntegration PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;243;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(PDESolver "/Users/patrikliba/CLionProjects/CppFM/build/test_pde_solver")
set_tests_properties(PDESolver PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;244;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(PDEStability "/Users/patrikliba/CLionProjects/CppFM/build/test_pde_stability")
set_tests_properties(PDEStability PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;245;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(TridiagonalSolver "/Users/patrikliba/CLionProjects/CppFM/build/test_tridiagonal_solver")
set_tests_properties(TridiagonalSolver PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;246;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
add_test(PDEGreeks "/Users/patrikliba/CLionProjects/CppFM/build/test_pde_greeks")
set_tests_properties(PDEGreeks PROPERTIES  _BACKTRACE_TRIPLES "/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;247;add_test;/Users/patrikliba/CLionProjects/CppFM/CMakeLists.txt;0;")
subdirs("_deps/googletest-build")
