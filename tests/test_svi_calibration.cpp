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
