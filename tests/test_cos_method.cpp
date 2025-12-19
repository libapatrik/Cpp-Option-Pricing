/**
 * Google Test Suite for COS Method
 * 
 * Tests:
 * 1. CDF Recovery Accuracy
 * 2. PDF Recovery Accuracy
 * 3. CDF Inversion Accuracy
 * 4. Performance Benchmarks for different N
 */

#include <gtest/gtest.h>
#include "../Utils.h"
#include <cmath>
#include <chrono>
#include <iomanip>


/**
 *  NOTES:
 *  (1) Newton's Method Comparison: Original vs. Optimized
 *  (2) Maximum error when recovering f(x) from φ(ω) by Fourier-cosine expansion (PDF Recovery)
 *  (3)
 *  
           
 */
// ============================================================================
// Test Fixture
// ============================================================================

class COSMethodTest : public ::testing::Test {
protected:
    // COS parameters as suggested in the paper by Fang & Oosterlee (2008)
    // For N(0,1): c1 = mean, c2 = variance, c4 = fourth CUMULANT (not moment!)
    // Note: Fourth cumulant κ₄ = 0 for normal distribution (excess kurtosis)
    //       Fourth moment E[X⁴] = 3, but we need the CUMULANT here!
    double c1 = 0.0;      // mean
    double c2 = 1.0;      // variance  
    double c4 = 0.0;      // fourth cumulant (NOT fourth moment!)
    double L = 10.0;      // truncation parameter
    double a = c1 - L * std::sqrt(c2 + std::sqrt(std::max(c4, 1e-10)));  // avoid sqrt of negative
    double b = c1 + L * std::sqrt(c2 + std::sqrt(std::max(c4, 1e-10)));
    
    // Standard Normal N(0,1) characteristic function - use Utils implementation
    // std::function<std::complex<double>(double)>  Utils::stdNormChF = Utils::stdNormChF;
};


/// TODO: Count the number of iterations taken until convergence
/// TODO: Compare the invertCDF vs. invertCDF_Optimized methods!
TEST_F(COSMethodTest, Comparenewtons) {
    std::vector<size_t> N_values = {64, 128, 256, 512};
    // Test at various quantiles including median and tails
    std::vector<double> test_probs = {0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};
    
    std::cout << "\n════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  Newton's Method Comparison: Original vs. Optimized\n";
    std::cout << "  Testing CDF Inversion at " << test_probs.size() << " probability points\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    std::cout << "ORIGINAL METHOD (recomputes coefficients each iteration):\n";
    std::cout << std::setw(8) << "N" 
              << std::setw(18) << "Avg Time (ms)" 
              << std::setw(18) << "Avg Iterations" << std::endl;
    std::cout << std::string(44, '-') << std::endl;
    
    // Store results for comparison table
    std::vector<double> original_times;
    std::vector<double> original_iters;
    
    for (size_t N : N_values) {
        size_t total_numIters = 0;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (double p : test_probs) {
            auto [quantile, numIters] = Transforms::invertCDF(a, b, N,  Utils::stdNormChF, p);
            (void)quantile;  // Suppress unused warning
            total_numIters += numIters;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        
        double ms_avg = duration.count() / test_probs.size();
        double numIters_avg = static_cast<double>(total_numIters) / test_probs.size();
        
        original_times.push_back(ms_avg);
        original_iters.push_back(numIters_avg);
        
        std::cout << std::setw(7) << N << " "
                  << std::setw(18) << std::fixed << std::setprecision(5) << ms_avg
                  << std::setw(18) << std::setprecision(2) << numIters_avg << std::endl;
    }
    
    std::cout << "\n";
    std::cout << "OPTIMIZED METHOD (caches coefficients):\n";
    std::cout << std::setw(8) << "N" 
              << std::setw(18) << "Avg Time (ms)" 
              << std::setw(18) << "Avg Iterations" << std::endl;
    std::cout << std::string(44, '-') << std::endl;
    
    // Store results for comparison
    std::vector<double> optimized_times;
    std::vector<double> optimized_iters;
    
    for (size_t N : N_values) {
        size_t total_numIters = 0;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (double p : test_probs) {
            auto [quantile, numIters] = Transforms::invertCDF_Optimized(a, b, N,  Utils::stdNormChF, p);
            (void)quantile;  // Suppress unused warning
            total_numIters += numIters;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        
        double ms_avg = duration.count() / test_probs.size();
        double numIters_avg = static_cast<double>(total_numIters) / test_probs.size();
        
        optimized_times.push_back(ms_avg);
        optimized_iters.push_back(numIters_avg);
        
        std::cout << std::setw(7) << N << " "
                  << std::setw(18) << std::fixed << std::setprecision(5) << ms_avg
                  << std::setw(18) << std::setprecision(2) << numIters_avg << std::endl;
    }
    
    std::cout << "\n";
    std::cout << "COMPARISON SUMMARY:\n";
    std::cout << std::setw(8) << "N"
              << std::setw(15) << "Speedup"
              << std::setw(20) << "Iter. Difference" << std::endl;
    std::cout << std::string(43, '-') << std::endl;
    
    for (size_t i = 0; i < N_values.size(); ++i) {
        double speedup = original_times[i] / optimized_times[i];
        double iter_diff = original_iters[i] - optimized_iters[i];
        
        std::cout << std::setw(7) << N_values[i] << " "
                  << std::setw(14) << std::fixed << std::setprecision(2) << speedup << "x"
                  << std::setw(20) << std::showpos << std::setprecision(2) << iter_diff 
                  << std::noshowpos << std::endl;
    }
}

/**
 *  INTERPRETATION:
 *   • Speedup:         How many times faster Optimized vs. Original Newton method
 *   • Iter. Difference: Change in iteration count (should be ~0 if same algorithm)
 *   • Observed:        Speedup of 2.8-3.5x due to coefficient caching
 *   • Note:            Utils.h claims 5-10x, but actual speedup depends on:
 *                      - Number of Newton iterations (fewer = less caching benefit)
 *                      - ChF complexity (simple N(0,1) is fast to evaluate)
 *   • Both methods:    Converge in identical number of iterations
 * ════════════════════════════════════════════════════════════════════════
 */
// =============================================================================================
// Test: Accuracy vs N
// Replicating Table 1 in Fang & Oosterlee (2008), "A Novel Pricing Method for European Options"
// =============================================================================================
/// TEST: Compare the COS PDF recovery at each point

TEST_F(COSMethodTest, AccuracyVsN) {
    std::vector<size_t> N_values = {4, 8, 16, 32, 64, 128, 256};
    
    // Test points should span the truncation range [a, b]
    const size_t num_test_points = 1000;
    std::vector<double> test_points(num_test_points);
    double step = (b - a) / (num_test_points - 1);
    for (size_t i = 0; i < num_test_points; ++i) {
        test_points[i] = a + i * step;
    }
    
    std::cout << "\n════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  TABLE 1: Maximum error when recovering f(x) from φ(ω)\n";
    std::cout << "           by Fourier-cosine expansion (PDF Recovery)\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "Distribution: N(0,1) - Standard Normal\n";
    std::cout << "Truncation:   [" << a << ", " << b << "]\n";
    std::cout << "Test points:  " << num_test_points << " evenly spaced in [" << a << ", " << b << "]\n";
    std::cout << "\n";
    std::cout << std::setw(6) << "N" 
              << std::setw(20) << "Error" 
              << std::setw(22) << "CPU time (ms)"
              << std::setw(25) << "Diff. in CPU (ms)" << std::endl;
    std::cout << std::string(73, '-') << std::endl;
    
    double prev_time = 0.0;
    
    for (size_t N : N_values) {
        // Run multiple times to get stable timing
        const int num_runs = 100;
        auto start = std::chrono::high_resolution_clock::now();
        
        std::vector<double> pdf_cos_values;
        for (int run = 0; run < num_runs; ++run) {
            pdf_cos_values = Transforms::recoverPDF(a, b, N,  Utils::stdNormChF, test_points);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        // Use duration<double, milli> for fractional milliseconds directly
        std::chrono::duration<double, std::milli> duration = end - start;
        double time_ms = duration.count() / num_runs;  // Average time per run in ms
        
        // Compute maximum error over all test points
        double max_error = 0.0;
        
        for (size_t i = 0; i < test_points.size(); ++i) {
            double x = test_points[i];
            double true_pdf_val = Utils::stdNormPdf(x);
            double error = std::abs(pdf_cos_values[i] - true_pdf_val);
            
            if (error > max_error) {
                max_error = error;
            }
        }
        
        // Calculate time difference from previous N
        double time_diff = (prev_time > 0.0) ? (time_ms - prev_time) : 0.0;
        
        std::cout << std::fixed << std::setw(5) << N 
                  << std::scientific << std::setprecision(6) << std::setw(20) << max_error
                  << std::fixed << std::setprecision(6) << std::setw(22) << time_ms;
        
        if (prev_time > 0.0) {
            std::cout << std::setw(25) << std::showpos << time_diff << std::noshowpos;
        } else {
            std::cout << std::setw(25) << "NaN";
        }
        
        std::cout << std::endl;
        prev_time = time_ms;
    }
}

/**
*  Notes:
*   • Error:          Maximum absolute error over all 1000 test points
*   • CPU time:       Average time to recover PDF at 1000 points (milliseconds)
*   • Diff. in CPU:   Change in CPU time compared to previous N value
*   • Expected:       Exponential convergence in error
*                     Linear increase in CPU time with N
* ════════════════════════════════════════════════════════════════════════
* 
*/


// We can make it even faster by precomputingCoefficients
TEST_F(COSMethodTest, AccuracyVsN_faster) {
    std::vector<size_t> N_values = {4, 8, 16, 32, 64, 128, 256};

    const size_t num_test_points = 200;
    std::vector<double> test_points(num_test_points);
    double step = (b - a) / (num_test_points - 1);
    for (size_t i = 0; i < num_test_points; ++i) {
        test_points[i] = a + i * step;
    }

    std::cout << "\n════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  TABLE 1: Maximum error when recovering f(x) from φ(ω)\n";
    std::cout << "           by Fourier-cosine expansion (PDF Recovery)\n";
    std::cout << "           with precomputed coefficients\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "Distribution: N(0,1) - Standard Normal\n";
    std::cout << "Truncation:   [" << a << ", " << b << "]\n";
    std::cout << "Test points:  " << num_test_points << " evenly spaced in [" << a << ", " << b << "]\n";
    std::cout << "\n";
    std::cout << std::setw(6) << "N"
              << std::setw(20) << "Error"
              << std::setw(22) << "CPU time (ms)"
              << std::setw(25) << "Diff. in CPU (ms)" << std::endl;
    std::cout << std::string(73, '-') << std::endl;

    double prev_time = 0.0;

    for (size_t N : N_values) {
        // PRECOMPUTE COEFFICIENTS ONCE (outside the timing loop)
        auto coeffs = Transforms::precomputeCoefficients(a, b, N, Utils::stdNormChF);

        // Run multiple times to get stable timing
        const int num_runs = 10;
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<double> pdf_cos_values(test_points.size());
        for (int run = 0; run < num_runs; ++run) {
            // Fast evaluation using precomputed coefficients
            for (size_t i = 0; i < test_points.size(); ++i) {
                pdf_cos_values[i] = Transforms::evaluatePDF(coeffs, test_points[i]);
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        double time_ms = duration.count() / num_runs;

        // Compute maximum error
        double max_error = 0.0;
        for (size_t i = 0; i < test_points.size(); ++i) {
            double true_pdf_val = Utils::stdNormPdf(test_points[i]);
            double error = std::abs(pdf_cos_values[i] - true_pdf_val);
            max_error = std::max(max_error, error);
        }

        double time_diff = (prev_time > 0.0) ? (time_ms - prev_time) : 0.0;

        std::cout << std::fixed << std::setw(5) << N
                  << std::scientific << std::setprecision(6) << std::setw(20) << max_error
                  << std::fixed << std::setprecision(6) << std::setw(22) << time_ms;

        if (prev_time > 0.0) {
            std::cout << std::setw(25) << std::showpos << time_diff << std::noshowpos;
        } else {
            std::cout << std::setw(25) << "NaN";
        }

        std::cout << std::endl;
        prev_time = time_ms;
    }
}

/// TODO: Add test for the bessel_function.py - just compare the modifiedBessel in Utils, against the data generated in py-file