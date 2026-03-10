/**
 * Empirical N sweep: accuracy vs N for cumulant-based COS bounds
 *
 * Two questions:
 * 1 COSPricer (Heston calibration) — how low can N go?
 * 2 CDF inversion (BK scheme integrated variance) — how low can N go?


 * Takeaways:
	1 Reduce N from 128 to 48
	2 In LM keep N=256 for precision low kappa regimes require it!
	3 Grid search: from N=128 to 96

 */

#include <chrono>
#include <cmath>
#include <cppfm/cos/COS.h>
#include <cppfm/pricers/COSPricer.h>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>

// ============================================================================
// Q1: COSPricer — Heston option pricing
// ============================================================================

TEST(NSweep, COSPricerSingleStrike)
{
	// test multiple parameter regimes
	struct Regime
	{
		const char *name;
		double kappa, vbar, sigma, rho, v0, r, T;
		double S0, K;
	};

	std::vector<Regime> regimes = {
		// standard Heston (Fang & Oosterlee benchmark)
		{"Standard ATM", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.0, 1.0, 100.0, 100.0},
		{"Standard OTM", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.0, 1.0, 100.0, 120.0},
		{"Standard ITM", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.0, 1.0, 100.0, 80.0},
		// high vol-of-vol
		{"High volvol ATM", 2.0, 0.04, 1.0, -0.7, 0.04, 0.02, 1.0, 100.0, 100.0},
		{"High volvol OTM", 2.0, 0.04, 1.0, -0.7, 0.04, 0.02, 1.0, 100.0, 130.0},
		// long maturity
		{"Long T=5 ATM", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.02, 5.0, 100.0, 100.0},
		{"Long T=5 OTM", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.02, 5.0, 100.0, 150.0},
		// low mean reversion (hard case)
		{"Low kappa ATM", 0.3, 0.06, 0.6, -0.6, 0.04, 0.0, 2.0, 100.0, 100.0},
		// short maturity
		{"Short T=0.1", 1.5768, 0.0398, 0.5751, -0.5711, 0.0175, 0.0, 0.1, 100.0, 100.0},
	};

	std::vector<size_t> Ns = {32, 48, 64, 96, 128, 256};

	std::cout << "\n=== COSPricer: Call Price Error vs N (cumulant bounds, L=10) ===\n";
	std::cout << std::setw(20) << "Regime";
	for (auto N : Ns)
		std::cout << std::setw(12) << ("N=" + std::to_string(N));
	std::cout << std::setw(12) << "Ref(512)" << "\n";
	std::cout << std::string(20 + 12 * (Ns.size() + 1), '-') << "\n";

	for (auto &reg : regimes)
	{
		HestonCF cf(reg.kappa, reg.vbar, reg.sigma, reg.rho, reg.v0, reg.r, reg.T);
		auto chfFunc = [&cf](double u)
		{ return cf(u); };
		auto cum = cf.cumulants();

		// reference: N=512 with cumulant bounds
		double ref = COSPricer::callPrice(reg.S0, reg.K, reg.r, reg.T, chfFunc, 512, 10.0,
										  cum.c1, cum.c2, cum.c4);

		std::cout << std::setw(20) << reg.name;
		for (auto N : Ns)
		{
			double p = COSPricer::callPrice(reg.S0, reg.K, reg.r, reg.T, chfFunc, N, 10.0,
											cum.c1, cum.c2, cum.c4);
			double err = std::abs(p - ref);
			std::cout << std::setw(12) << std::scientific << std::setprecision(2) << err;
		}
		std::cout << std::setw(12) << std::fixed << std::setprecision(6) << ref << "\n";
	}
}

TEST(NSweep, COSPricerMultiStrike)
{
	// calibration context: multiple strikes per slice
	double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
	double rho = -0.5711, v0 = 0.0175, r = 0.0, S0 = 100.0;
	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> strikes = {80, 85, 90, 95, 100, 105, 110, 115, 120};

	std::vector<size_t> Ns = {32, 48, 64, 96, 128, 256};

	std::cout << "\n=== COSPricer Multi-Strike: Max Abs Error vs N ===\n";
	std::cout << std::setw(10) << "T";
	for (auto N : Ns)
		std::cout << std::setw(12) << ("N=" + std::to_string(N));
	std::cout << "\n"
			  << std::string(10 + 12 * Ns.size(), '-') << "\n";

	for (double T : maturities)
	{
		HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
		auto chfFunc = [&cf](double u)
		{ return cf(u); };
		auto cum = cf.cumulants();

		auto ref = COSPricer::callPrices(S0, strikes, r, T, chfFunc, 512, 10.0,
										 cum.c1, cum.c2, cum.c4);

		std::cout << std::setw(10) << std::fixed << std::setprecision(2) << T;
		for (auto N : Ns)
		{
			auto prices = COSPricer::callPrices(S0, strikes, r, T, chfFunc, N, 10.0,
												cum.c1, cum.c2, cum.c4);
			double maxErr = 0;
			for (size_t i = 0; i < strikes.size(); ++i)
				maxErr = std::max(maxErr, std::abs(prices[i] - ref[i]));

			std::cout << std::setw(12) << std::scientific << std::setprecision(2) << maxErr;
		}
		std::cout << "\n";
	}
}

// ============================================================================
// Q2: CDF Inversion — BK scheme integrated variance
// ============================================================================

TEST(NSweep, CDFInversionIntegratedVariance)
{
	// BK scheme parameters
	struct BKCase
	{
		const char *name;
		double kappa, vbar, sigma, v_s, v_t, tau;
	};

	std::vector<BKCase> cases = {
		{"Typical dt=0.01", 1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.01},
		{"Typical dt=0.05", 1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.05},
		{"Typical dt=0.25", 1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.25},
		{"High volvol", 2.0, 0.04, 1.0, 0.06, 0.08, 0.05},
		{"Low vol state", 1.5768, 0.0398, 0.5751, 0.005, 0.008, 0.05},
		{"High vol state", 1.5768, 0.0398, 0.5751, 0.15, 0.12, 0.05},
	};

	// test at multiple quantiles
	std::vector<double> quantiles = {0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99};
	std::vector<size_t> Ns = {32, 48, 64, 96, 128};

	std::cout << "\n=== BK CDF Inversion: Max Abs Error vs N (across quantiles) ===\n";
	std::cout << std::setw(20) << "Case";
	for (auto N : Ns)
		std::cout << std::setw(12) << ("N=" + std::to_string(N));
	std::cout << std::setw(14) << "Ref range" << "\n";
	std::cout << std::string(20 + 12 * Ns.size() + 14, '-') << "\n";

	for (auto &c : cases)
	{
		auto chf = [&](double omega) -> std::complex<double>
		{
			return ChFIntegratedVariance::compute(omega, c.kappa, c.vbar, c.sigma, c.v_s, c.v_t, c.tau);
		};

		auto cum = ChFIntegratedVariance::cumulants(c.kappa, c.vbar, c.sigma, c.v_s, c.v_t, c.tau);
		auto [a_full, b_full] = computeTruncationBounds(cum.c1, cum.c2, cum.c4, 10.0);
		double a = std::max(0.0, a_full);
		double b = b_full;

		// reference: N=256
		std::vector<double> ref_vals;
		for (double q : quantiles)
		{
			auto [val, _] = Transforms::invertCDF(a, b, 256, chf, q);
			ref_vals.push_back(std::max(0.0, val));
		}

		double ref_min = *std::min_element(ref_vals.begin(), ref_vals.end());
		double ref_max = *std::max_element(ref_vals.begin(), ref_vals.end());

		std::cout << std::setw(20) << c.name;
		for (auto N : Ns)
		{
			double maxErr = 0;
			for (size_t i = 0; i < quantiles.size(); ++i)
			{
				auto [val, iters] = Transforms::invertCDF(a, b, N, chf, quantiles[i]);
				val = std::max(0.0, val);
				maxErr = std::max(maxErr, std::abs(val - ref_vals[i]));
			}
			std::cout << std::setw(12) << std::scientific << std::setprecision(2) << maxErr;
		}
		std::cout << std::setw(14) << std::scientific << std::setprecision(2)
				  << "[" << ref_min << ", " << ref_max << "]"
				  << "  bounds=[" << a << ", " << b << "]"
				  << "  c1=" << cum.c1 << " c2=" << cum.c2 << "\n";
	}
}

TEST(NSweep, CFDiagnostics)
{
	struct Case
	{
		const char *name;
		double kappa, vbar, sigma, v_s, v_t, tau;
	};

	std::vector<Case> cases = {
		{"Typical dt=0.05", 1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.05},
		{"Typical dt=0.01", 1.5768, 0.0398, 0.5751, 0.04, 0.05, 0.01},
		{"High vol state", 1.5768, 0.0398, 0.5751, 0.15, 0.12, 0.05},
	};

	std::cout << "\n=== CF diagnostics: ln(phi(h)) at various h ===\n";
	for (auto &c : cases)
	{
		std::cout << "\n--- " << c.name << " ---\n";
		std::cout << std::setw(10) << "h"
				  << std::setw(22) << "Re(ln phi(h))"
				  << std::setw(22) << "Im(ln phi(h))"
				  << std::setw(22) << "|phi(h)|"
				  << std::setw(16) << "c2 estimate\n";

		for (double h : {1e-1, 5e-2, 1e-2, 7e-3, 5e-3, 1e-3, 1e-4})
		{
			auto phi_p = ChFIntegratedVariance::compute(h, c.kappa, c.vbar, c.sigma, c.v_s, c.v_t, c.tau);
			auto phi_m = ChFIntegratedVariance::compute(-h, c.kappa, c.vbar, c.sigma, c.v_s, c.v_t, c.tau);
			auto lnp = std::log(phi_p);
			auto lnm = std::log(phi_m);
			double c2_est = -std::real(lnp + lnm) / (h * h);

			std::cout << std::setw(10) << std::scientific << std::setprecision(1) << h
					  << std::setw(22) << std::scientific << std::setprecision(8) << std::real(lnp)
					  << std::setw(22) << std::scientific << std::setprecision(8) << std::imag(lnp)
					  << std::setw(22) << std::scientific << std::setprecision(8) << std::abs(phi_p)
					  << std::setw(22) << std::scientific << std::setprecision(4) << c2_est << "\n";
		}
	}
}

// also measure relative errors
TEST(NSweep, CDFInversionRelativeError)
{
	double kappa = 1.5768, vbar = 0.0398, sigma = 0.5751;
	double v_s = 0.04, v_t = 0.05, tau = 0.05;
	std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99};
	std::vector<size_t> Ns = {32, 48, 64, 96, 128};

	auto chf = [&](double omega) -> std::complex<double>
	{
		return ChFIntegratedVariance::compute(omega, kappa, vbar, sigma, v_s, v_t, tau);
	};

	auto cum = ChFIntegratedVariance::cumulants(kappa, vbar, sigma, v_s, v_t, tau);
	auto [a_full, b_full] = computeTruncationBounds(cum.c1, cum.c2, cum.c4, 10.0);
	double a = std::max(0.0, a_full);
	double b = b_full;

	std::cout << "\n=== BK CDF Inversion: Relative Error by Quantile (typical case) ===\n";
	std::cout << std::setw(10) << "Quantile";
	for (auto N : Ns)
		std::cout << std::setw(12) << ("N=" + std::to_string(N));
	std::cout << std::setw(14) << "Ref value" << "\n";
	std::cout << std::string(10 + 12 * Ns.size() + 14, '-') << "\n";

	for (double q : quantiles)
	{
		auto [ref, _] = Transforms::invertCDF(a, b, 256, chf, q);
		ref = std::max(1e-15, ref);

		std::cout << std::setw(10) << std::fixed << std::setprecision(2) << q;
		for (auto N : Ns)
		{
			auto [val, iters] = Transforms::invertCDF(a, b, N, chf, q);
			val = std::max(0.0, val);
			double relErr = std::abs(val - ref) / ref;
			std::cout << std::setw(12) << std::scientific << std::setprecision(2) << relErr;
		}
		std::cout << std::setw(14) << std::scientific << std::setprecision(4) << ref << "\n";
	}
}
