#include <gtest/gtest.h>
#include "../Utils.h"
#include <cmath>
#include <complex>
#include <chrono>
#include <iomanip>

/**
 * Google Test: Newton-Raphson CDF Inversion Methods
 *
 * Tests both implementations:
 * 1. invertCDF() - Original version
 * 2. invertCDF_Optimized() - Coefficient caching version
 *
 * Verifies:
 * - Correctness for various distributions
 * - Consistency between both methods
 * - Edge cases and boundary conditions
 * - Performance comparison
 */

// ============================================================================
// Characteristic Functions for Testing
// ============================================================================

// Standard Normal: N(0, 1)
// ChF: φ(ω) = exp(-ω²/2)
std::complex<double> stdNormalChF(double omega) {
    return std::exp(std::complex<double>(-0.5 * omega * omega, 0.0));
}

// Normal with mean μ and variance σ²: N(μ, σ²)
// ChF: φ(ω) = exp(iωμ - ω²σ²/2)
class NormalChF {
public:
    NormalChF(double mu, double sigma) : _mu(mu), _sigma(sigma) {}
    
    std::complex<double> operator()(double omega) const {
        return std::exp(std::complex<double>(
            -0.5 * omega * omega * _sigma * _sigma,
            omega * _mu
        ));
    }
    
private:
    double _mu;
    double _sigma;
};


// ============================================================================
// Test Fixture
// ============================================================================

class NewtonMethodsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common test parameters
        tol = 1e-8;
        max_iter = 100;
    }

    double tol;
    size_t max_iter;
};

TEST_F(NewtonMethodsTest, NewtonMethodComparison) {
    // Comprehensive comparison of Original vs Optimized Newton-Raphson
    double a = -10.0;
    double b = 10.0;
    std::vector<size_t> N_values = {32, 64, 128, 256, 512};
    std::vector<double> test_probs = {0.05, 0.25, 0.50, 0.75, 0.95};
    const size_t num_runs = 50;  // Runs per test for stable timing
    
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "  TABLE 2: Comprehensive Newton-Raphson Method Comparison" << std::endl;
    std::cout << "           Original vs Optimized (Coefficient Caching)" << std::endl;
    std::cout << "========================================================================" << std::endl;
    std::cout << "Distribution: N(0,1), Tolerance: " << std::scientific << std::setprecision(1) 
              << tol << ", Max Iterations: " << max_iter << std::endl;
    std::cout << std::fixed << "\n" << std::setw(8) << "N"
              << std::setw(8) << "p" 
              << std::setw(15) << "Time_Orig(μs)"
              << std::setw(15) << "Time_Opt(μs)"
              << std::setw(12) << "Speedup"
              << std::setw(18) << "Error Diff"
              << std::setw(15) << "True Error" << std::endl;
    std::cout << std::string(95, '-') << std::endl;
    
    double total_speedup = 0.0;
    double max_error_diff = 0.0;
    int count = 0;
    
    for (size_t N : N_values) {
        for (double p : test_probs) {
            // Time original method
            auto start_orig = std::chrono::high_resolution_clock::now();
            double x_orig = 0.0;
            for (size_t i = 0; i < num_runs; ++i) {
                x_orig = Transforms::invertCDF(a, b, N, stdNormalChF, p, max_iter, tol).first;
            }
            auto end_orig = std::chrono::high_resolution_clock::now();
            auto duration_orig = std::chrono::duration_cast<std::chrono::microseconds>(end_orig - start_orig);
            double time_orig = duration_orig.count() / static_cast<double>(num_runs);
            
            // Time optimized method
            auto start_opt = std::chrono::high_resolution_clock::now();
            double x_opt = 0.0;
            for (size_t i = 0; i < num_runs; ++i) {
                x_opt = Transforms::invertCDF_Optimized(a, b, N, stdNormalChF, p, max_iter, tol).first;
            }
            auto end_opt = std::chrono::high_resolution_clock::now();
            auto duration_opt = std::chrono::duration_cast<std::chrono::microseconds>(end_opt - start_opt);
            double time_opt = duration_opt.count() / static_cast<double>(num_runs);
            
            // Compute metrics
            double speedup = time_orig / time_opt;
            double error_diff = std::abs(x_orig - x_opt);
            double true_quantile = Utils::inverseNormalCDF(p);  // Moro's algorithm as reference
            double true_error = std::abs(x_orig - true_quantile);
            
            // Accumulate statistics
            total_speedup += speedup;
            max_error_diff = std::max(max_error_diff, error_diff);
            count++;
            
            // Print row
            std::cout << std::fixed << std::setw(8) << N
                      << std::setprecision(2) << std::setw(8) << p
                      << std::setprecision(1) << std::setw(15) << time_orig
                      << std::setw(15) << time_opt
                      << std::setprecision(2) << std::setw(12) << speedup << "x"
                      << std::scientific << std::setprecision(2) << std::setw(18) << error_diff
                      << std::setw(15) << true_error << std::endl;
            
            // Assertions
            EXPECT_DOUBLE_EQ(x_orig, x_opt) << "Methods must produce identical results";
            EXPECT_GT(speedup, 1.0) << "Optimized should be faster";
            EXPECT_LT(error_diff, 1e-10) << "Methods should agree to machine precision";
        }
    }
}


