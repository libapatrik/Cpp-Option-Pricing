#include <chrono>
#include <cmath>
#include <cppfm/calibration/HestonCalibration.h>
#include <cppfm/calibration/ImpliedVolSolver.h>
#include <cppfm/cos/COS.h>
#include <cppfm/pricers/COSPricer.h>
#include <gtest/gtest.h>
#include <iostream>

// generate synthetic market IVs from known Heston params
static std::vector<HestonSliceData> generateSyntheticData(
	const HestonParams &trueParams, double S0, double r,
	const std::vector<double> &maturities,
	const std::vector<double> &strikes)
{
	std::vector<HestonSliceData> slices;
	double sigmaHint = std::sqrt(trueParams.vbar);

	for (double T : maturities)
	{
		HestonSliceData slice;
		slice.T = T;
		slice.strikes = strikes;

		HestonCF cf(trueParams.kappa, trueParams.vbar, trueParams.sigma_v,
					trueParams.rho, trueParams.v0, r, T);
		auto chfFunc = [&cf](double u)
		{ return cf(u); };

		auto prices = COSPricer::callPrices(S0, strikes, r, T,
											chfFunc, 256, 10.0, sigmaHint);

		for (size_t i = 0; i < strikes.size(); ++i)
		{
			double iv = ImpliedVolSolver::solve(prices[i], S0, strikes[i], r, T, true);
			slice.marketIVs.push_back(iv);
		}
		slices.push_back(std::move(slice));
	}
	return slices;
}

TEST(HestonCalibration, ParamsVectorRoundTrip)
{
	HestonParams p(0.05, 2.0, 0.06, 0.4, -0.65);
	auto v = p.toVector();

	ASSERT_EQ(v.size(), 5);
	EXPECT_DOUBLE_EQ(v[0], 0.05);
	EXPECT_DOUBLE_EQ(v[1], 2.0);
	EXPECT_DOUBLE_EQ(v[2], 0.06);
	EXPECT_DOUBLE_EQ(v[3], 0.4);
	EXPECT_DOUBLE_EQ(v[4], -0.65);

	HestonParams q(v);
	EXPECT_DOUBLE_EQ(q.v0, p.v0);
	EXPECT_DOUBLE_EQ(q.kappa, p.kappa);
	EXPECT_DOUBLE_EQ(q.vbar, p.vbar);
	EXPECT_DOUBLE_EQ(q.sigma_v, p.sigma_v);
	EXPECT_DOUBLE_EQ(q.rho, p.rho);
}

TEST(HestonCalibration, FellerCondition)
{
	// 2*1.5*0.04 = 0.12 > 0.09 = 0.3^2  -> Feller satisfied
	HestonParams p1(0.04, 1.5, 0.04, 0.3, -0.7);
	EXPECT_TRUE(p1.satisfiesFellerCondition());

	// 2*0.5*0.01 = 0.01 < 1.0 = 1.0^2  -> Feller violated
	HestonParams p2(0.01, 0.5, 0.01, 1.0, -0.5);
	EXPECT_FALSE(p2.satisfiesFellerCondition());
}

TEST(HestonCalibration, SyntheticRoundTrip)
{
	// true params
	HestonParams trueParams(0.04, 1.5, 0.04, 0.3, -0.7);
	double S0 = 100.0, r = 0.02;

	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> strikes = {80, 90, 95, 100, 105, 110, 120};

	auto slices = generateSyntheticData(trueParams, S0, r, maturities, strikes);

	// perturbed guess
	HestonParams guess(0.06, 2.5, 0.06, 0.5, -0.5);

	LMOptions opts;
	opts.maxIter = 200;
	opts.tol = 1e-10;

	auto result = calibrateHeston(slices, S0, r, guess, opts);

	EXPECT_TRUE(result.converged) << result.message;
	EXPECT_LT(result.rmse, 1e-4); // sub-bp

	// param recovery (loose on kappa — flat direction)
	EXPECT_NEAR(result.params.v0, trueParams.v0, 0.01);
	EXPECT_NEAR(result.params.kappa, trueParams.kappa, 0.5);
	EXPECT_NEAR(result.params.vbar, trueParams.vbar, 0.01);
	EXPECT_NEAR(result.params.sigma_v, trueParams.sigma_v, 0.05);
	EXPECT_NEAR(result.params.rho, trueParams.rho, 0.05);
}

TEST(HestonCalibration, PerSliceRmse)
{
	HestonParams trueParams(0.04, 1.5, 0.04, 0.3, -0.7);
	double S0 = 100.0, r = 0.02;

	std::vector<double> maturities = {0.5, 1.0};
	std::vector<double> strikes = {90, 100, 110};

	auto slices = generateSyntheticData(trueParams, S0, r, maturities, strikes);

	auto result = calibrateHeston(slices, S0, r);

	// per-slice diagnostics populated
	ASSERT_EQ(result.sliceRmse.size(), 2);
	for (auto rmse : result.sliceRmse)
		EXPECT_LT(rmse, 1e-3);
}

TEST(HestonCalibration, GridSearchInit)
{
	HestonParams trueParams(0.04, 1.5, 0.04, 0.3, -0.7);
	double S0 = 100.0, r = 0.02;

	std::vector<double> maturities = {0.5, 1.0};
	std::vector<double> strikes = {90, 100, 110};

	auto slices = generateSyntheticData(trueParams, S0, r, maturities, strikes);

	// grid search for starting point
	auto gridGuess = hestonGridSearch(slices, S0, r, 3);

	// refine with LM
	LMOptions opts;
	opts.maxIter = 200;
	auto result = calibrateHeston(slices, S0, r, gridGuess, opts);

	EXPECT_TRUE(result.converged) << result.message;
	EXPECT_LT(result.rmse, 1e-4);
}

TEST(HestonCalibration, GridSearchParallelSpeedup)
{
	HestonParams trueParams(0.04, 1.5, 0.04, 0.3, -0.7);
	double S0 = 100.0, r = 0.02;

	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> strikes = {80, 90, 95, 100, 105, 110, 120};

	auto slices = generateSyntheticData(trueParams, S0, r, maturities, strikes);

	int gp = 5; // 5^5 = 3125 combos

	// sequential
	auto t0 = std::chrono::high_resolution_clock::now();
	auto seqResult = hestonGridSearch(slices, S0, r, gp);
	auto t1 = std::chrono::high_resolution_clock::now();

	// parallel
	auto t2 = std::chrono::high_resolution_clock::now();
	auto parResult = hestonGridSearchParallel(slices, S0, r, gp);
	auto t3 = std::chrono::high_resolution_clock::now();

	double seqMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
	double parMs = std::chrono::duration<double, std::milli>(t3 - t2).count();
	double speedup = seqMs / parMs;

	std::cout << "Grid search (" << gp << "^5 = " << (int)std::pow(gp, 5) << " combos):\n"
			  << "  Sequential: " << seqMs << " ms\n"
			  << "  Parallel:   " << parMs << " ms\n"
			  << "  Speedup:    " << speedup << "x\n";

	// both should find the same optimum
	EXPECT_DOUBLE_EQ(seqResult.v0, parResult.v0);
	EXPECT_DOUBLE_EQ(seqResult.kappa, parResult.kappa);
	EXPECT_DOUBLE_EQ(seqResult.vbar, parResult.vbar);
	EXPECT_DOUBLE_EQ(seqResult.sigma_v, parResult.sigma_v);
	EXPECT_DOUBLE_EQ(seqResult.rho, parResult.rho);

	// parallel should be faster (conservative: at least 1.5x on any multi-core machine)
	EXPECT_GT(speedup, 1.5) << "Parallel version not significantly faster";
}
