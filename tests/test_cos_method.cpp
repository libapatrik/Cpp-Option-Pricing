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
#include <cppfm/cos/COS.h>
#include <cppfm/pricers/COSPricer.h>
#include <cppfm/utils/Utils.h>  // Utils::stdNormChF, stdNormPdf
#include <cmath>
#include <chrono>
#include <iomanip>
#include <tuple>


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
    // COS parameters - Fang & Oosterlee (2008)
    // N(0,1): c1 = mean, c2 = variance, c4 = fourth cumulant
    double L = 10.0;
    double a, b;
    void SetUp() override {
        std::tie(a, b) = computeTruncationBounds(0.0, 1.0, 0.0, L);
    }
};


// benchmark invertCDF convergence vs N
TEST_F(COSMethodTest, InvertCDFBenchmark) {
    std::vector<size_t> N_values = {64, 128, 256, 512};
    std::vector<double> test_probs = {0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};

    std::cout << "\n════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  CDF Inversion Benchmark\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";

    std::cout << std::setw(8) << "N"
              << std::setw(18) << "Avg Time (ms)"
              << std::setw(18) << "Avg Iterations" << std::endl;
    std::cout << std::string(44, '-') << std::endl;

    for (size_t N : N_values) {
        size_t total_numIters = 0;

        auto start = std::chrono::high_resolution_clock::now();

        for (double p : test_probs) {
            auto [quantile, numIters] = Transforms::invertCDF(a, b, N, Utils::stdNormChF, p);
            (void)quantile;
            total_numIters += numIters;
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;

        double ms_avg = duration.count() / test_probs.size();
        double numIters_avg = static_cast<double>(total_numIters) / test_probs.size();

        std::cout << std::setw(7) << N << " "
                  << std::setw(18) << std::fixed << std::setprecision(5) << ms_avg
                  << std::setw(18) << std::setprecision(2) << numIters_avg << std::endl;
    }
}
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

// ============================================================================
// Truncation Bounds Tests
// ============================================================================

TEST(TruncationBounds, NormalDistribution) {
    // N(0,1): c1=0, c2=1, c4=0
    double L = 10.0;
    auto [a, b] = computeTruncationBounds(0.0, 1.0, 0.0, L);

    EXPECT_NEAR(a, -L, 0.01);
    EXPECT_NEAR(b,  L, 0.01);
    EXPECT_LT(a, 0.0);
    EXPECT_GT(b, 0.0);
    EXPECT_NEAR(a + b, 0.0, 1e-10);
}

TEST(TruncationBounds, WithKurtosis) {
    double L = 10.0;
    auto [a1, b1] = computeTruncationBounds(0.0, 1.0, 0.0, L);
    auto [a2, b2] = computeTruncationBounds(0.0, 1.0, 9.0, L);

    // c4=9: sqrt(c2 + sqrt(c4)) = sqrt(1 + 3) = 2
    EXPECT_NEAR(a2, -20.0, 0.01);
    EXPECT_NEAR(b2,  20.0, 0.01);
    EXPECT_LT(a2, a1);
    EXPECT_GT(b2, b1);
}

TEST(TruncationBounds, NonzeroMean) {
    double L = 10.0;
    auto [a, b] = computeTruncationBounds(2.0, 1.0, 0.0, L);
    EXPECT_NEAR((a + b) / 2.0, 2.0, 0.01);
}

// ============================================================================
// Heston Analytical Cumulant Tests
// ============================================================================

TEST(HestonCumulants, AnalyticalVsNumerical) {
    // F&O Table 11 params
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.025, T = 1.0;

    HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
    auto cum = cf.cumulants();

    // validate against numerical cumulants at h=7e-3
    double h = 7e-3;
    auto lnphi = [&cf](double u) { return std::log(cf(u)); };

    std::complex<double> lp1 = lnphi(h);
    std::complex<double> lm1 = lnphi(-h);
    std::complex<double> lp2 = lnphi(2.0 * h);
    std::complex<double> lm2 = lnphi(-2.0 * h);

    double num_c1 = std::imag(lp1 - lm1) / (2.0 * h);
    double num_c2 = -std::real(lp1 + lm1) / (h * h);
    double num_c4 = std::real(lp2 - 4.0 * lp1 - 4.0 * lm1 + lm2) / (h * h * h * h);

    // c1 and c2 should match to high precision
    EXPECT_NEAR(cum.c1, num_c1, 1e-6) << "c1 mismatch";
    EXPECT_NEAR(cum.c2, num_c2, 1e-6) << "c2 mismatch";
    // c4 numerical at h=7e-3 has ~1% relative error
    EXPECT_NEAR(cum.c4, num_c4, std::abs(cum.c4) * 0.05) << "c4 mismatch";
}

TEST(HestonCumulants, SanityCheck) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.0, T = 1.0;

    HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
    auto cum = cf.cumulants();

    // c1 ~ (r - vbar/2)*T + correction, should be small negative
    EXPECT_LT(cum.c1, 0.0);
    EXPECT_GT(cum.c1, -0.1);

    // c2 = variance of log-returns, positive
    EXPECT_GT(cum.c2, 0.0);
    EXPECT_LT(cum.c2, 1.0);

    EXPECT_TRUE(std::isfinite(cum.c4));
}

TEST(HestonCumulants, BoundsAreReasonable) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.0, T = 10.0;

    HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
    auto cum = cf.cumulants();

    double x0 = 0.0; // ATM
    auto [a, b] = computeTruncationBounds(x0 + cum.c1, cum.c2, cum.c4, 10.0);

    EXPECT_TRUE(std::isfinite(a));
    EXPECT_TRUE(std::isfinite(b));
    EXPECT_LT(a, 0.0);
    EXPECT_GT(b, 0.0);
    EXPECT_GT(b - a, 0.1);
}

TEST(HestonCumulants, LongMaturityCumulantsGrowCorrectly) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.0;

    HestonCF cf1(kappa, vbar, sigma, rho, v0, r, 1.0);
    HestonCF cf10(kappa, vbar, sigma, rho, v0, r, 10.0);
    auto cum1 = cf1.cumulants();
    auto cum10 = cf10.cumulants();

    // variance grows with T
    EXPECT_GT(cum10.c2, cum1.c2);
}

// ============================================================================
// Integrated Variance Cumulant Tests
// ============================================================================

TEST(IntVarCumulants, PositiveMeanAndVariance) {
    double kappa = 1.5, vbar = 0.04, sigma = 0.3;
    double v_s = 0.04, v_t = 0.05, tau = 0.1;

    auto cum = ChFIntegratedVariance::cumulants(kappa, vbar, sigma, v_s, v_t, tau);

    // integrated variance has positive mean ~ V_avg * tau
    double V_avg = 0.5 * (v_s + v_t);
    EXPECT_GT(cum.c1, 0.0);
    EXPECT_NEAR(cum.c1, V_avg * tau, V_avg * tau * 0.5);

    EXPECT_GT(cum.c2, 0.0);
    EXPECT_EQ(cum.c4, 0.0);  // explicitly set to 0
}

// ============================================================================
// COSPricer Cumulant Bounds Tests
// ============================================================================

TEST(COSPricerCumulantBounds, HestonATMCall) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.0, T = 1.0;
    double S0 = 100.0, K = 100.0;

    HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
    auto chfFunc = [&cf](double u) { return cf(u); };
    auto cum = cf.cumulants();

    double price_cum = COSPricer::callPrice(S0, K, r, T, chfFunc, 128, 10.0,
                                             cum.c1, cum.c2, cum.c4);

    // reference: sigma-hint with many terms
    double sigmaHint = std::sqrt(std::max(v0, vbar));
    double price_ref = COSPricer::callPrice(S0, K, r, T, chfFunc, 512, 10.0, sigmaHint);

    EXPECT_NEAR(price_cum, price_ref, 0.01);
    EXPECT_GT(price_cum, 0.0);
}

TEST(COSPricerCumulantBounds, HestonMultiStrike) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double rho = -0.5711, v0 = 0.0175, r = 0.0, T = 1.0;
    double S0 = 100.0;
    std::vector<double> strikes = {80, 90, 100, 110, 120};

    HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
    auto chfFunc = [&cf](double u) { return cf(u); };
    auto cum = cf.cumulants();

    auto prices_cum = COSPricer::callPrices(S0, strikes, r, T, chfFunc, 128, 10.0,
                                             cum.c1, cum.c2, cum.c4);

    double sigmaHint = std::sqrt(std::max(v0, vbar));
    auto prices_ref = COSPricer::callPrices(S0, strikes, r, T, chfFunc, 512, 10.0, sigmaHint);

    for (size_t i = 0; i < strikes.size(); ++i) {
        EXPECT_NEAR(prices_cum[i], prices_ref[i], 0.02) << "Strike " << strikes[i];
        EXPECT_GT(prices_cum[i], 0.0);
    }
}