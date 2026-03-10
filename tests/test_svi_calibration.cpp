#include <cmath>
#include <cppfm/market/DiscountCurve.h>
#include <cppfm/market/SviCalibration.h>
#include <gtest/gtest.h>

TEST(SviCalibrationTest, TotalVarianceTest)
{
	// Test w(k) = a + b * (rho*(k-m) + sqrt((k-m)^2 + sigma^2))
	// a=0.04, b=0.1, rho=-0.5, m=0, sigma=0.2, k=0.1
	// dkm = 0.1
	// rho*dkm = -0.05
	// sqrt(0.01 + 0.04) = sqrt(0.05) = 0.223607
	// w = 0.04 + 0.1*(-0.05 + 0.223607) = 0.0573607
	SviParams params(0.04, 0.01, -0.5, 0.0, 0.2);
	double dkm = 0.1; // k-m
	double expected = 0.04 + 0.01 * (-0.5 * dkm + std::sqrt(dkm * dkm + 0.2 * 0.2));

	EXPECT_NEAR(params.totalVariance(dkm), expected, 1e-12);
}

TEST(SviCalibrationTest, toFromVectorTest)
{
	// {a, b, rho, m, sigma}
	SviParams original(0.04, 0.01, -0.5, 0.0, 0.2);
	std::vector<double> vector = original.toVector();

	ASSERT_EQ(vector.size(), 5);
	EXPECT_DOUBLE_EQ(vector[0], 0.04);
	EXPECT_DOUBLE_EQ(vector[1], 0.01);
	EXPECT_DOUBLE_EQ(vector[2], -0.5);
	EXPECT_DOUBLE_EQ(vector[3], 0.0);
	EXPECT_DOUBLE_EQ(vector[4], 0.2);

	SviParams reconstructed;
	reconstructed.fromVector(vector);
	EXPECT_DOUBLE_EQ(reconstructed.a, original.a);
	EXPECT_DOUBLE_EQ(reconstructed.b, original.b);
	EXPECT_DOUBLE_EQ(reconstructed.rho, original.rho);
	EXPECT_DOUBLE_EQ(reconstructed.m, original.m);
	EXPECT_DOUBLE_EQ(reconstructed.sigma, original.sigma);

	// vector constructor
	SviParams fromVec(vector);
	EXPECT_DOUBLE_EQ(fromVec.a, original.a);
	EXPECT_DOUBLE_EQ(fromVec.b, original.b);
	EXPECT_DOUBLE_EQ(fromVec.rho, original.rho);
	EXPECT_DOUBLE_EQ(fromVec.m, original.m);
	EXPECT_DOUBLE_EQ(fromVec.sigma, original.sigma);
}

TEST(SviCalibrationTest, SyntheticCalibrationTest)
{
	// test: known SVI -> market vols -> calibrateSmile -> recover params match
	// if this works, then the rest of the pipeline works
	SviParams trueParams(0.04, 0.08, -0.4, 0.02, 0.15);
	double forward = 100.0;
	double T = 1.0;
	// generate strikes
	std::vector<double> strikes = {80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0};
	// generate market vols
	std::vector<double> marketVols(strikes.size());

	for (size_t i = 0; i < strikes.size(); ++i)
	{
		double k = std::log(strikes[i] / forward);
		marketVols[i] = trueParams.impliedVol(k, T);
	}
	// calibrate SVI
	SviParams guess;
	auto result = calibrateSmile(guess, strikes, marketVols, forward, T);

	EXPECT_TRUE(result.converged);
	EXPECT_LT(result.rmse, 1e-6);

	auto *fitted = dynamic_cast<SviParams *>(result.params.get());
	ASSERT_NE(fitted, nullptr);

	EXPECT_NEAR(fitted->a, trueParams.a, 1e-6);
	EXPECT_NEAR(fitted->b, trueParams.b, 1e-6);
	EXPECT_NEAR(fitted->rho, trueParams.rho, 1e-6);
	EXPECT_NEAR(fitted->m, trueParams.m, 1e-6);
	EXPECT_NEAR(fitted->sigma, trueParams.sigma, 1e-6);

	// fitted smile should reproduce the market vols exactly
	for (size_t i = 0; i < strikes.size(); ++i)
	{
		double k = std::log(strikes[i] / forward);
		EXPECT_NEAR(fitted->impliedVol(k, T), marketVols[i], 1e-6)
			<< "vol mismatch at strike K = " << strikes[i];
	}
}

// ============================================================================
// TEST SSVI
// ============================================================================

TEST(SsviTest, PhiPowerLaw)
{
	// phi(theta) = eta / (theta^gamma * (1+theta)^(1-gamma))
	SsviParams p(-0.5, 0.8, 0.5);
	p.theta = 0.04;
	double expected = 0.8 / (std::pow(0.04, 0.5) * std::pow(1.04, 0.5));
	EXPECT_NEAR(p.phi(), expected, 1e-12);
}

TEST(SsviTest, TotalVarianceFormula)
{
	// w(k) = (theta/2) * (1 + rho*phi*k + sqrt((phi*k + rho)^2 + 1 - rho^2))
	SsviParams p(-0.4, 0.6, 0.5);
	p.theta = 0.04;
	double k = 0.1;
	double ph = p.phi();
	double pk = ph * k;
	double expected = (0.04 / 2.0) * (1.0 + (-0.4) * pk + std::sqrt((pk - 0.4) * (pk - 0.4) + 1.0 - 0.16));
	EXPECT_NEAR(p.totalVariance(k), expected, 1e-12);
}

TEST(SsviTest, ToFromVector)
{
	SsviParams original(-0.4, 0.6, 0.5);
	auto v = original.toVector();
	ASSERT_EQ(v.size(), 3);

	SsviParams copy;
	copy.fromVector(v);
	EXPECT_DOUBLE_EQ(copy.rho, original.rho);
	EXPECT_DOUBLE_EQ(copy.eta, original.eta);
	EXPECT_DOUBLE_EQ(copy.gamma, original.gamma);
}

TEST(SsviTest, Constraints)
{
	// valid params: eta*(1+|rho|) <= 2
	SsviParams good(-0.5, 0.8, 0.5);
	EXPECT_TRUE(good.satisfiesConstraints());

	// violate wing bound: eta*(1+|rho|) = 1.5*(1+0.5) = 2.25 > 2
	SsviParams bad(-0.5, 1.5, 0.5);
	EXPECT_FALSE(bad.satisfiesConstraints());
}

TEST(SsviTest, ButterflyArbitrageFree)
{
	// g(k) >= 0 across a range of strikes => no butterfly arb
	SsviParams p(-0.4, 0.6, 0.5);
	p.theta = 0.04;

	for (double k = -1.0; k <= 1.0; k += 0.05)
	{
		EXPECT_GE(p.gFunction(k), 0.0)
			<< "butterfly arb at k=" << k;
	}
}

TEST(SsviTest, CalendarSpreadArbitrageFree)
{
	// total variance must be non-decreasing in theta for each k --- TODO still needs fix in src
	// if theta1 < theta2, then w(k; theta1) <= w(k; theta2)
	SsviParams p(-0.4, 0.6, 0.5);

	std::vector<double> thetas = {0.01, 0.04, 0.09, 0.16, 0.25};

	for (double k = -0.5; k <= 0.5; k += 0.1)
	{
		double prevW = 0.0;
		for (double th : thetas)
		{
			p.theta = th;
			double w = p.totalVariance(k);
			EXPECT_GE(w, prevW)
				<< "calendar arb at k=" << k << ", theta=" << th;
			prevW = w;
		}
	}
}

TEST(SsviTest, DerivativesVsNumerical)
{
	// finite-difference check on dw, d2w
	SsviParams p(-0.4, 0.6, 0.5);
	p.theta = 0.04;

	double k = 0.1;
	double h = 1e-6;
	double numDw = (p.totalVariance(k + h) - p.totalVariance(k - h)) / (2.0 * h);
	double numD2w = (p.totalVariance(k + h) - 2.0 * p.totalVariance(k) + p.totalVariance(k - h)) / (h * h);

	EXPECT_NEAR(p.dw(k), numDw, 1e-5);
	EXPECT_NEAR(p.d2w(k), numD2w, 1e-4);
}

TEST(SsviTest, SyntheticCalibration)
{
	// generate market data from known SSVI, calibrate back, check recovery
	SsviParams trueParams(-0.4, 0.6, 0.5);
	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	std::vector<double> strikes = {80, 85, 90, 95, 100, 105, 110, 115, 120};

	std::vector<std::vector<double>> strikesPerMat(maturities.size(), strikes);
	std::vector<std::vector<double>> volsPerMat(maturities.size());

	// generate synthetic vols from the true SSVI
	for (size_t i = 0; i < maturities.size(); ++i)
	{
		double T = maturities[i];
		// ATM vol ~ 20% => theta ~ 0.04*T
		// use consistent theta: ATM total var at each maturity
		double atmVol = 0.20;
		trueParams.theta = atmVol * atmVol * T;

		volsPerMat[i].resize(strikes.size());
		for (size_t j = 0; j < strikes.size(); ++j)
		{
			double k = std::log(strikes[j] / forwards[i]);
			double w = trueParams.totalVariance(k);
			volsPerMat[i][j] = std::sqrt(w / T);
		}
	}

	auto result = calibrateSsvi(strikesPerMat, volsPerMat, forwards, maturities);

	EXPECT_TRUE(result.converged);
	EXPECT_LT(result.rmse, 1e-4);

	// fitted params should be close
	EXPECT_NEAR(result.params.rho, trueParams.rho, 0.05);
	EXPECT_NEAR(result.params.eta, trueParams.eta, 0.1);
	EXPECT_NEAR(result.params.gamma, trueParams.gamma, 0.1);

	// reproduced vols should match
	SsviParams fitted = result.params;
	for (size_t i = 0; i < maturities.size(); ++i)
	{
		fitted.theta = result.thetas[i];
		double T = maturities[i];
		for (size_t j = 0; j < strikes.size(); ++j)
		{
			double k = std::log(strikes[j] / forwards[i]);
			double fittedVol = std::sqrt(fitted.totalVariance(k) / T);
			EXPECT_NEAR(fittedVol, volsPerMat[i][j], 1e-3)
				<< "T=" << T << ", K=" << strikes[j];
		}
	}
}

TEST(SsviTest, ThetaMonotonicity)
{
	// non-monotone thetas -- PAVA should fix them
	SsviParams trueParams(-0.3, 0.5, 0.5);
	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	// theta[1] < theta[0] violates monotonicity
	std::vector<double> trueThetas = {0.03, 0.02, 0.04, 0.08};

	std::vector<double> strikeGrid = {80, 85, 90, 95, 100, 105, 110, 115, 120};
	std::vector<std::vector<double>> strikes(4, strikeGrid);
	std::vector<std::vector<double>> vols(4);

	for (int i = 0; i < 4; ++i)
	{
		trueParams.theta = trueThetas[i];
		vols[i].resize(strikeGrid.size());
		for (size_t j = 0; j < strikeGrid.size(); ++j)
		{
			double k = std::log(strikeGrid[j] / forwards[i]);
			vols[i][j] = std::sqrt(trueParams.totalVariance(k) / maturities[i]);
		}
	}

	SsviCalibrationOptions opts;
	opts.enforceMonotonicity = true;
	opts.lmOpts.gradTol = 1e-6;
	auto result = calibrateSsvi(strikes, vols, forwards, maturities, opts);

	EXPECT_TRUE(result.converged);
	for (size_t i = 1; i < result.thetas.size(); ++i)
		EXPECT_GE(result.thetas[i], result.thetas[i - 1]);
}

TEST(SsviTest, AnalyticalJacobianConsistency)
{
	// verify analytical partials match finite differences
	SsviParams p(-0.3, 0.5, 0.5);
	p.theta = 0.04;
	double k = 0.15;
	double h = 1e-7;

	double w0 = p.totalVariance(k);

	SsviParams pr = p;
	pr.rho += h;
	double dw_drho_fd = (pr.totalVariance(k) - w0) / h;

	SsviParams pe = p;
	pe.eta += h;
	double dw_deta_fd = (pe.totalVariance(k) - w0) / h;

	SsviParams pg = p;
	pg.gamma += h;
	double dw_dgamma_fd = (pg.totalVariance(k) - w0) / h;

	// analytical
	double ph = p.phi();
	double pk = ph * k;
	double D = std::sqrt((pk + p.rho) * (pk + p.rho) + (1.0 - p.rho * p.rho));
	double dw_dk = (p.theta / 2.0) * ph * (p.rho + (pk + p.rho) / D);

	double dw_drho_an = (p.theta / 2.0) * pk * (1.0 + 1.0 / D);
	double dw_deta_an = (k / p.eta) * dw_dk;
	double dw_dgamma_an = k * std::log((1.0 + p.theta) / p.theta) * dw_dk;

	EXPECT_NEAR(dw_drho_an, dw_drho_fd, 1e-4);
	EXPECT_NEAR(dw_deta_an, dw_deta_fd, 1e-4);
	EXPECT_NEAR(dw_dgamma_an, dw_dgamma_fd, 1e-4);
}

TEST(SsviTest, GridSearchImprovesFit)
{
	// strong skew far from default guess -- grid search should help
	SsviParams trueParams(-0.7, 0.8, 0.7);
	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	std::vector<double> trueThetas = {0.01, 0.02, 0.04, 0.08};

	std::vector<double> strikeGrid = {80, 85, 90, 95, 100, 105, 110, 115, 120};
	std::vector<std::vector<double>> strikes(4, strikeGrid);
	std::vector<std::vector<double>> vols(4);

	for (int i = 0; i < 4; ++i)
	{
		trueParams.theta = trueThetas[i];
		vols[i].resize(strikeGrid.size());
		for (size_t j = 0; j < strikeGrid.size(); ++j)
		{
			double k = std::log(strikeGrid[j] / forwards[i]);
			vols[i][j] = std::sqrt(trueParams.totalVariance(k) / maturities[i]);
		}
	}

	SsviCalibrationOptions optsGrid;
	optsGrid.useGridSearch = true;
	auto resultGrid = calibrateSsvi(strikes, vols, forwards, maturities, optsGrid);

	SsviCalibrationOptions optsNoGrid;
	optsNoGrid.useGridSearch = false;
	auto resultNoGrid = calibrateSsvi(strikes, vols, forwards, maturities, optsNoGrid);

	EXPECT_LE(resultGrid.rmse, resultNoGrid.rmse + 1e-10);
	EXPECT_TRUE(resultGrid.converged);
}

// ============================================================================
// SsviSurface tests
// ============================================================================

// helper: calibrate synthetic SSVI + build surface
static SsviCalibrationResult makeSyntheticResult()
{
	SsviParams trueParams(-0.4, 0.6, 0.5);
	std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	std::vector<double> strikes = {80, 85, 90, 95, 100, 105, 110, 115, 120};

	std::vector<std::vector<double>> strikesPerMat(maturities.size(), strikes);
	std::vector<std::vector<double>> volsPerMat(maturities.size());

	for (size_t i = 0; i < maturities.size(); ++i)
	{
		double T = maturities[i];
		double atmVol = 0.20;
		trueParams.theta = atmVol * atmVol * T;

		volsPerMat[i].resize(strikes.size());
		for (size_t j = 0; j < strikes.size(); ++j)
		{
			double k = std::log(strikes[j] / forwards[i]);
			volsPerMat[i][j] = std::sqrt(trueParams.totalVariance(k) / T);
		}
	}

	return calibrateSsvi(strikesPerMat, volsPerMat, forwards, maturities);
}

TEST(SsviSurfaceTest, DwdthetaVsFiniteDifference)
{
	SsviParams p(-0.4, 0.6, 0.5);
	p.theta = 0.04;
	double k = 0.15;
	double h = 1e-7;

	double w0 = p.totalVariance(k);
	p.theta = 0.04 + h;
	double w1 = p.totalVariance(k);
	double fd = (w1 - w0) / h;

	p.theta = 0.04;
	EXPECT_NEAR(p.dwdtheta(k), fd, 1e-4);
}

TEST(SsviSurfaceTest, DwdthetaAtKZero)
{
	// at k=0, w(0) = theta, so dw/dtheta = 1
	SsviParams p(-0.4, 0.6, 0.5);
	p.theta = 0.04;
	EXPECT_NEAR(p.dwdtheta(0.0), 1.0, 1e-12);
}

TEST(SsviSurfaceTest, LocalVolPositive)
{
	auto result = makeSyntheticResult();
	ASSERT_TRUE(result.converged);

	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	auto surface = result.buildAnalyticalSurface(forwards);

	for (double spot = 80.0; spot <= 120.0; spot += 5.0)
	{
		for (double T = 0.3; T <= 1.5; T += 0.1)
		{
			double lv = surface.localVolatility(spot, T);
			EXPECT_GT(lv, 0.0)
				<< "non-positive local vol at spot=" << spot << ", T=" << T;
			EXPECT_LT(lv, 2.0)
				<< "unreasonably large local vol at spot=" << spot << ", T=" << T;
		}
	}
}

TEST(SsviSurfaceTest, ImpliedVolMatchesSsviParams)
{
	auto result = makeSyntheticResult();
	ASSERT_TRUE(result.converged);

	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	auto surface = result.buildAnalyticalSurface(forwards);
	SsviParams fitted = result.params;

	// at calibrated maturities, IV should match direct SsviParams computation
	for (size_t i = 0; i < result.maturities.size(); ++i)
	{
		double T = result.maturities[i];
		fitted.theta = result.thetas[i];

		for (double K = 85.0; K <= 115.0; K += 5.0)
		{
			double k = std::log(K / forwards[i]);
			double ivDirect = std::sqrt(fitted.totalVariance(k) / T);
			double ivSurface = surface.impliedVolatility(K, T);
			EXPECT_NEAR(ivSurface, ivDirect, 1e-10)
				<< "IV mismatch at K=" << K << ", T=" << T;
		}
	}
}

TEST(SsviSurfaceTest, CompareToGriddedDupire)
{
	auto result = makeSyntheticResult();
	ASSERT_TRUE(result.converged);

	std::vector<double> forwards = {100.0, 100.0, 100.0, 100.0};
	auto surface = result.buildAnalyticalSurface(forwards);

	// build gridded surface + numerical Dupire for comparison
	std::vector<double> strikeGrid;
	for (double K = 70.0; K <= 130.0; K += 1.0)
		strikeGrid.push_back(K);

	FlatDiscountCurve dc(0.0); // zero rates for simplicity
	auto griddedSurface = result.buildSurface(strikeGrid, forwards, dc);
	ASSERT_NE(griddedSurface, nullptr);

	// compare at a few points -- analytical vs numerical shouldn't diverge too much
	// tolerance is loose since numerical Dupire has FD noise
	for (double spot : {90.0, 95.0, 100.0, 105.0, 110.0})
	{
		for (double T : {0.5, 1.0})
		{
			double lvAnalytical = surface.localVolatility(spot, T);
			double lvNumerical = griddedSurface->localVolatility(spot, T);

			EXPECT_NEAR(lvAnalytical, lvNumerical, 0.05)
				<< "spot=" << spot << ", T=" << T
				<< " analytical=" << lvAnalytical << " numerical=" << lvNumerical;
		}
	}
}
