#include <gtest/gtest.h>
#include "../Utils.h"
#include <cmath>
#include <iostream>

/**
 * Test suite for ThomasAlgorithm
 */

// ============================================================================
// Test 1: Simple 3x3 system
// ============================================================================
TEST(TridiagonalSolverTest, Simple3x3System)
{
    // System:
    //   2x_0 + 1x_1         = 1
    //   1x_0 + 2x_1 + 1x_2  = 2
    //         1x_1  + 2x_2  = 1
    //
    // Exact solution: x = [0, 1, 0]
    
    std::vector<double> lower = {0.0, 1.0, 1.0};  // a[0] unused
    std::vector<double> diag  = {2.0, 2.0, 2.0};
    std::vector<double> upper = {1.0, 1.0, 0.0};  // c[2] unused
    std::vector<double> rhs   = {1.0, 2.0, 1.0};
    
    auto x = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], 0.0, 1e-10);
    EXPECT_NEAR(x[1], 1.0, 1e-10);
    EXPECT_NEAR(x[2], 0.0, 1e-10);
    
    std::cout << "Simple 3x3 test passed: x = [" 
              << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
}

// ============================================================================
// Test 2: Identity system
// ============================================================================
TEST(TridiagonalSolverTest, IdentitySystem)
{
    // System: I·x = d, so x = d
    
    std::vector<double> lower = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> diag  = {1.0, 1.0, 1.0, 1.0}; // Diagonal elements
    std::vector<double> upper = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> rhs   = {1.0, 2.0, 3.0, 4.0};
    
    auto x = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    ASSERT_EQ(x.size(), 4);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(x[i], rhs[i], 1e-10);
    }
}

// ============================================================================
// Test 3: Tridiagonal Laplacian (from finite differences)
// ============================================================================
TEST(TridiagonalSolverTest, LaplacianMatrix)
{
    // Discrete Laplacian: -u_{i-1} + 2u_i - u_{i+1} = f_i
    // With boundary conditions u_0 = 0, u_5 = 0
    
    size_t n = 4;  // Interior points
    std::vector<double> lower(n, -1.0);
    std::vector<double> diag(n, 2.0);
    std::vector<double> upper(n, -1.0);
    std::vector<double> rhs(n, 1.0);  // Constant forcing
    
    lower[0] = 0.0;  // Boundary
    upper[n-1] = 0.0;  // Boundary
    
    auto x = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    ASSERT_EQ(x.size(), n);
    
    // Solution should be symmetric and peaked in the middle
    EXPECT_NEAR(x[0], x[3], 1e-10) << "Solution should be symmetric";
    EXPECT_NEAR(x[1], x[2], 1e-10) << "Solution should be symmetric";
    EXPECT_GT(x[1], x[0]) << "Solution should peak in middle";
    
    std::cout << "Laplacian test solution: [";
    for (size_t i = 0; i < n; ++i) {
        std::cout << x[i] << (i < n-1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}

// ============================================================================
// Test 4: Large system
// ============================================================================
TEST(TridiagonalSolverTest, LargeSystem)
{
    // Test with N=100 to ensure algorithm scales well
    size_t N = 100;
    
    // Create a tridiagonal system with known solution
    // A·x = b where x = [1, 2, 3, ..., N]
    std::vector<double> x_exact(N);
    for (size_t i = 0; i < N; ++i) {
        x_exact[i] = static_cast<double>(i + 1);
    }
    
    // Construct A as a simple tridiagonal matrix
    std::vector<double> lower(N, -1.0);
    std::vector<double> diag(N, 3.0);
    std::vector<double> upper(N, -1.0);
    
    lower[0] = 0.0;
    upper[N-1] = 0.0;
    
    // Compute b = A·x_exact
    std::vector<double> rhs(N, 0.0);
    for (size_t i = 0; i < N; ++i) {
        rhs[i] = diag[i] * x_exact[i];
        if (i > 0) {
            rhs[i] += lower[i] * x_exact[i-1];
        }
        if (i < N-1) {
            rhs[i] += upper[i] * x_exact[i+1];
        }
    }
    
    // Solve A·x = b
    auto x = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    ASSERT_EQ(x.size(), N);
    
    // Verify solution matches x_exact
    double max_error = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double error = std::abs(x[i] - x_exact[i]);
        max_error = std::max(max_error, error);
    }
    
    EXPECT_LT(max_error, 1e-10) << "Max error in large system: " << max_error;
    std::cout << "Large system (N=" << N << ") max error: " << max_error << std::endl;
}

// ============================================================================
// Test 5: Error handling - singular matrix
// ============================================================================
TEST(TridiagonalSolverTest, SingularMatrix)
{
    // Matrix with zero diagonal element
    std::vector<double> lower = {0.0, 1.0, 1.0};
    std::vector<double> diag  = {0.0, 2.0, 2.0};  // First element is zero 
    std::vector<double> upper = {1.0, 1.0, 0.0};
    std::vector<double> rhs   = {1.0, 2.0, 1.0};
    
    EXPECT_THROW(
        ThomasAlgorithm::solve(lower, diag, upper, rhs),
        std::runtime_error
    ) << "Should throw on singular matrix";
}

// ============================================================================
// Test 6: Error handling - size mismatch
// ============================================================================
TEST(TridiagonalSolverTest, SizeMismatch)
{
    std::vector<double> lower = {0.0, 1.0};      // Wrong size
    std::vector<double> diag  = {2.0, 2.0, 2.0};
    std::vector<double> upper = {1.0, 1.0, 0.0};
    std::vector<double> rhs   = {1.0, 2.0, 1.0};
    
    EXPECT_THROW(
        ThomasAlgorithm::solve(lower, diag, upper, rhs),
        std::invalid_argument) << "Should throw on size mismatch";
}

// ============================================================================
// Test 7: Single equation
// ============================================================================
TEST(TridiagonalSolverTest, SingleEquation)
{
    // Trivial case: 1x1 system
    // 2*x_0 = 4  ->  x_0 = 2
    
    std::vector<double> lower = {0.0};
    std::vector<double> diag  = {2.0};
    std::vector<double> upper = {0.0};
    std::vector<double> rhs   = {4.0};
    
    auto x = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 2.0, 1e-10);
}

// int main(int argc, char** argv)
// {
//     testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }

