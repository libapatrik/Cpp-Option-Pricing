#include <gtest/gtest.h>
#include <cppfm/cos/COS.h>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>

// ============================================================================
// Item 3: Bessel function |phi| > 1 diagnostic
// CF is a ratio I_nu(z_R)/I_nu(z_kappa) * prefactors
// |phi(omega)| must be <= 1 for all omega — violation means Bessel is broken
// ============================================================================

TEST(BesselDiagnostics, PhiMagnitudeScan) {
    struct Case {
        const char* name;
        double kappa, vbar, sigma, v_s, v_t, tau;
    };

    std::vector<Case> cases = {
        {"OK: dt=0.05",      1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.05},
        {"Bad: dt=0.01",     1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.01},
        {"Bad: high vol",    1.5768, 0.0398, 0.5751, 0.15, 0.12, 0.05},
        {"OK: high volvol",  2.0,    0.04,   1.0,    0.06, 0.08, 0.05},
        {"Edge: dt=0.02",    1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.02},
        {"Edge: dt=0.03",    1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.03},
    };

    std::cout << "\n=== |phi(omega)| scan — should be <= 1 everywhere ===\n";

    for (auto& c : cases) {
        double max_abs = 0;
        double worst_omega = 0;

        for (int i = -5000; i <= 5000; ++i) {
            double omega = i * 0.1;
            auto phi = ChFIntegratedVariance::compute(omega, c.kappa, c.vbar, c.sigma, c.v_s, c.v_t, c.tau);
            double abs_phi = std::abs(phi);
            if (abs_phi > max_abs) {
                max_abs = abs_phi;
                worst_omega = omega;
            }
        }

        std::cout << std::setw(20) << c.name
                  << "  max|phi| = " << std::fixed << std::setprecision(12) << max_abs
                  << "  at omega = " << std::fixed << std::setprecision(1) << worst_omega
                  << (max_abs > 1.0 + 1e-10 ? "  *** BROKEN ***" : "  ok")
                  << "\n";
    }
}

TEST(BesselDiagnostics, BrokenCaseZoomIn) {
    // for the broken dt=0.01 case, trace the Bessel argument regime
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double v_s = 0.04, v_t = 0.05, tau = 0.01;
    double d = 4.0 * kappa * vbar / (sigma * sigma);
    double nu = 0.5 * d - 1.0;

    std::cout << "\n--- dt=0.01: Bessel arg regime (nu=" << std::fixed << std::setprecision(3) << nu << ") ---\n";
    std::cout << std::setw(10) << "omega"
              << std::setw(14) << "|phi|"
              << std::setw(12) << "|z_R|"
              << std::setw(12) << "Re(z_R)"
              << std::setw(12) << "Im(z_R)"
              << std::setw(12) << "z_kappa"
              << std::setw(14) << "method\n";

    for (double omega : {0.001, 0.01, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0, 200.0, 500.0}) {
        std::complex<double> R = std::sqrt(kappa * kappa - 2.0 * sigma * sigma * std::complex<double>(0, omega));
        std::complex<double> exp_R_tau = std::exp(-R * tau);
        double sqrt_prod = std::sqrt(v_t * v_s);
        std::complex<double> z_R = sqrt_prod * 4.0 * R * std::exp(-R * tau / 2.0) /
                                    (sigma * sigma * (1.0 - exp_R_tau));

        double exp_kappa_tau = std::exp(-kappa * tau);
        double z_kappa = sqrt_prod * 4.0 * kappa * std::exp(-kappa * tau / 2.0) /
                         (sigma * sigma * (1.0 - exp_kappa_tau));

        auto phi = ChFIntegratedVariance::compute(omega, kappa, vbar, sigma, v_s, v_t, tau);

        // which method does modifiedBesselI use?
        const char* method = (std::abs(z_R) < 1e-10) ? "zero" :
                             (std::abs(z_R) > 30.0 && std::real(z_R) > 0.0) ? "asymptotic" : "power_series";

        std::cout << std::setw(10) << std::fixed << std::setprecision(3) << omega
                  << std::setw(14) << std::setprecision(10) << std::abs(phi)
                  << std::setw(12) << std::setprecision(2) << std::abs(z_R)
                  << std::setw(12) << std::setprecision(2) << std::real(z_R)
                  << std::setw(12) << std::setprecision(2) << std::imag(z_R)
                  << std::setw(12) << std::setprecision(2) << z_kappa
                  << std::setw(14) << method << "\n";
    }
}

// ============================================================================
// Item 3 continued: asymptotic expansion accuracy
// current: 1-term, I_nu(z) ~ e^z/sqrt(2piz) * [1 - (4nu^2-1)/(8z)]
// test: how many terms needed to fix |phi| > 1?
// ============================================================================

TEST(BesselDiagnostics, AsymptoticTermsNeeded) {
    // Hankel asymptotic: I_nu(z) ~ e^z/sqrt(2piz) * sum_{k=0}^{K} (-1)^k * a_k / z^k
    // a_0 = 1, a_k = prod_{j=0}^{k-1} (4*nu^2 - (2j+1)^2) / (8k)
    // i.e. a_k = ((4nu^2-1)(4nu^2-9)...(4nu^2-(2k-1)^2)) / (k! * 8^k)

    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double d = 4.0 * kappa * vbar / (sigma * sigma);
    double nu = 0.5 * d - 1.0;
    double mu = 4.0 * nu * nu;
    constexpr double PI = std::numbers::pi;

    // z_kappa for dt=0.01 case (real, ~54)
    double v_s = 0.04, v_t = 0.05, tau = 0.01;
    double sqrt_prod = std::sqrt(v_t * v_s);
    double exp_kappa_tau = std::exp(-kappa * tau);
    double z_real = sqrt_prod * 4.0 * kappa * std::exp(-kappa * tau / 2.0) /
                    (sigma * sigma * (1.0 - exp_kappa_tau));

    // boost reference (real arg)
    // we can't call modifiedBesselI directly (private), but we can compare
    // asymptotic expansion terms against boost for the real z_kappa

    std::cout << "\n=== Asymptotic expansion terms at z=" << std::fixed << std::setprecision(2) << z_real
              << ", nu=" << std::setprecision(3) << nu << " ===\n";
    std::cout << std::setw(8) << "Terms" << std::setw(20) << "Correction" << std::setw(20) << "Cumulative\n";

    // compute a_k coefficients
    double a_k = 1.0;
    double cumulative = 1.0;
    std::cout << std::setw(8) << 0 << std::setw(20) << std::scientific << std::setprecision(10) << a_k
              << std::setw(20) << cumulative << "\n";

    for (int k = 1; k <= 8; ++k) {
        double factor = (mu - (2.0*k - 1.0) * (2.0*k - 1.0)) / (8.0 * k * z_real);
        a_k *= -factor;  // (-1)^k from the alternating series
        cumulative += a_k;
        std::cout << std::setw(8) << k << std::setw(20) << std::scientific << std::setprecision(10) << a_k
                  << std::setw(20) << cumulative << "\n";
    }

    // now test with complex z_R at omega=0 (should equal z_kappa)
    // and at small omega where Im(z_R) is tiny
    std::cout << "\n=== |phi| with K asymptotic terms (dt=0.01, omega=0.001) ===\n";
    std::cout << "  Current 1-term expansion gives |phi| = 1 + O(1e-5)\n";
    std::cout << "  Each additional term reduces error by factor ~1/(8*|z|) ≈ 1/" << 8.0*z_real << "\n";
    std::cout << "  Expected: 3-4 terms should bring |phi| to machine precision\n";

    // compute z_R at omega where it's nearly real
    double omega = 0.001;
    std::complex<double> R = std::sqrt(kappa * kappa - 2.0 * sigma * sigma * std::complex<double>(0, omega));
    std::complex<double> exp_R_tau = std::exp(-R * tau);
    std::complex<double> z_R = sqrt_prod * 4.0 * R * std::exp(-R * tau / 2.0) /
                                (sigma * sigma * (1.0 - exp_R_tau));

    std::cout << "  z_R = " << std::fixed << std::setprecision(6) << std::real(z_R)
              << " + " << std::imag(z_R) << "i\n";
    std::cout << "  z_kappa = " << z_real << "\n";
    std::cout << "  |z_R - z_kappa| = " << std::scientific << std::abs(z_R - z_real) << "\n";

    // the error in I_nu(z_R)/I_nu(z_kappa) from K-term expansion
    // since z_R ≈ z_kappa, the ratio is ~1, and error comes from the
    // asymptotic truncation error in each, which partially cancels
    // but the cancellation is imperfect due to tiny Im(z_R)
    std::cout << "\n  At |z|=" << std::fixed << std::setprecision(0) << z_real
              << ", term k has magnitude ~ (mu-(2k-1)^2) / (8*" << z_real << ")^k\n";
    std::cout << "  mu = 4*nu^2 = " << std::setprecision(4) << mu << "\n";
    std::cout << "  term 2 magnitude: " << std::scientific << std::setprecision(2)
              << std::abs((mu - 1.0) * (mu - 9.0)) / (128.0 * z_real * z_real) << "\n";
    std::cout << "  term 3 magnitude: "
              << std::abs((mu - 1.0) * (mu - 9.0) * (mu - 25.0)) / (128.0 * 6.0 * 8.0 * z_real * z_real * z_real) << "\n";
}

// ============================================================================
// Item 4: Newton vs Bisection cost for CDF inversion
// ============================================================================

TEST(InversionCost, IterationCounts) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double v_s = 0.04, v_t = 0.05, tau = 0.05;

    auto chf = [&](double omega) -> std::complex<double> {
        return ChFIntegratedVariance::compute(omega, kappa, vbar, sigma, v_s, v_t, tau);
    };

    auto cum = ChFIntegratedVariance::cumulants(kappa, vbar, sigma, v_s, v_t, tau);
    auto [a_full, b_full] = computeTruncationBounds(cum.c1, cum.c2, cum.c4, 10.0);
    double a = std::max(0.0, a_full);
    double b = b_full;
    size_t N = 48;

    std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99};

    // Newton: iteration counts
    std::cout << "\n=== Newton iteration counts (5 bisect init + Newton) ===\n";
    std::cout << std::setw(10) << "Quantile" << std::setw(12) << "Iters" << std::setw(14) << "Value\n";

    for (double q : quantiles) {
        auto [val, iters] = Transforms::invertCDF(a, b, N, chf, q);
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << q
                  << std::setw(12) << iters
                  << std::setw(14) << std::scientific << std::setprecision(6) << val << "\n";
    }

    // pure bisection: iteration counts for same tolerance
    std::cout << "\n=== Pure bisection iteration counts (tol=1e-10 relative) ===\n";
    auto coeffs = Transforms::precomputeCoefficients(a, b, N, chf);

    std::cout << std::setw(10) << "Quantile" << std::setw(12) << "Iters"
              << std::setw(14) << "Value\n";

    for (double q : quantiles) {
        double lo = a, hi = b;
        size_t bisect_iters = 0;
        while (hi - lo > 1e-10 * std::max(1.0, std::abs(lo))) {
            double mid = 0.5 * (lo + hi);
            if (Transforms::evaluateCDF(coeffs, mid) < q)
                lo = mid;
            else
                hi = mid;
            ++bisect_iters;
            if (bisect_iters > 200) break;
        }
        double val = 0.5 * (lo + hi);
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << q
                  << std::setw(12) << bisect_iters
                  << std::setw(14) << std::scientific << std::setprecision(6) << val << "\n";
    }
}

TEST(InversionCost, TimingComparison) {
    double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
    double v_s = 0.04, v_t = 0.05, tau = 0.05;

    auto chf = [&](double omega) -> std::complex<double> {
        return ChFIntegratedVariance::compute(omega, kappa, vbar, sigma, v_s, v_t, tau);
    };

    auto cum = ChFIntegratedVariance::cumulants(kappa, vbar, sigma, v_s, v_t, tau);
    auto [a_full, b_full] = computeTruncationBounds(cum.c1, cum.c2, cum.c4, 10.0);
    double a = std::max(0.0, a_full);
    double b = b_full;
    size_t N = 48;

    size_t reps = 1000;
    std::vector<double> test_quantiles;
    for (size_t i = 0; i < reps; ++i)
        test_quantiles.push_back(0.01 + 0.98 * i / (reps - 1));

    // Newton timing (current: precompute per call)
    auto t0 = std::chrono::high_resolution_clock::now();
    double newton_sum = 0;
    for (double q : test_quantiles) {
        auto [val, _] = Transforms::invertCDF(a, b, N, chf, q);
        newton_sum += val;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double newton_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // bisection (recompute coefficients each call, same as Newton)
    auto t2 = std::chrono::high_resolution_clock::now();
    double bisect_sum = 0;
    for (double q : test_quantiles) {
        auto c = Transforms::precomputeCoefficients(a, b, N, chf);
        double lo = a, hi = b;
        for (int i = 0; i < 50; ++i) {
            double mid = 0.5 * (lo + hi);
            if (Transforms::evaluateCDF(c, mid) < q)
                lo = mid;
            else
                hi = mid;
            if (hi - lo < 1e-10 * std::max(1.0, std::abs(lo)))
                break;
        }
        bisect_sum += 0.5 * (lo + hi);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    double bisect_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // bisection with shared coefficients (best case — same CF params)
    auto t4 = std::chrono::high_resolution_clock::now();
    auto shared_coeffs = Transforms::precomputeCoefficients(a, b, N, chf);
    double bisect2_sum = 0;
    for (double q : test_quantiles) {
        double lo = a, hi = b;
        for (int i = 0; i < 50; ++i) {
            double mid = 0.5 * (lo + hi);
            if (Transforms::evaluateCDF(shared_coeffs, mid) < q)
                lo = mid;
            else
                hi = mid;
            if (hi - lo < 1e-10 * std::max(1.0, std::abs(lo)))
                break;
        }
        bisect2_sum += 0.5 * (lo + hi);
    }
    auto t5 = std::chrono::high_resolution_clock::now();
    double bisect2_ms = std::chrono::duration<double, std::milli>(t5 - t4).count();

    std::cout << "\n=== Timing: " << reps << " inversions (N=" << N << ") ===\n";
    std::cout << "  Newton (5 bisect + Newton):  " << std::fixed << std::setprecision(2) << newton_ms << " ms\n";
    std::cout << "  Bisection (recompute each):  " << bisect_ms << " ms\n";
    std::cout << "  Bisection (shared coeffs):   " << bisect2_ms << " ms\n";
    std::cout << "  Newton/Bisect(recomp) ratio:  " << newton_ms / bisect_ms << "x\n";
    std::cout << "  Newton/Bisect(shared) ratio:  " << newton_ms / bisect2_ms << "x\n";

    // verify they agree
    EXPECT_NEAR(newton_sum, bisect_sum, 1e-4);
    EXPECT_NEAR(newton_sum, bisect2_sum, 1e-4);
}
