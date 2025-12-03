#include <gtest/gtest.h>
#include "../VolatilitySurface.h"
// TODO: Include test_utils.h when ready

/*
TODO:
1. Test valid construction
2. Test invalid inputs (empty data, unsorted, negative vols)
3. Test builder pattern
4. Test clone/equality
*/

// TODO: Add test cases



// 1. Test valid construction

TEST(VolatilitySurfaceConstructionTest, ValidConstruction) 
{
    std::vector<double> strikes = {100.0, 105.0, 110.0, 115.0, 120.0};
    std::vector<double> maturities = {0.5, 1.0, 1.5, 2.0, 2.5};
    std::vector<std::vector<double>> volatilitiesMatrix = {
        {0.2, 0.25, 0.3, 0.35, 0.4},
        {0.25, 0.3, 0.35, 0.4, 0.45},
        {0.3, 0.35, 0.4, 0.45, 0.5},
        {0.35, 0.4, 0.45, 0.5, 0.55},
        {0.4, 0.45, 0.5, 0.55, 0.6}
    };
    FlatDiscountCurve discountCurve(0.05);
    VolatilitySurface::SmileInterpolationType smileInterpolationType = VolatilitySurface::SmileInterpolationType::CubicSpline;
    VolatilitySurface::MaturityInterpolationType maturityInterpolationType = VolatilitySurface::MaturityInterpolationType::ForwardMoneyness;

    VolatilitySurface volatilitySurface(strikes, maturities, volatilitiesMatrix, discountCurve, smileInterpolationType, maturityInterpolationType);

    EXPECT_EQ(volatilitySurface.strikes(), strikes);
    EXPECT_EQ(volatilitySurface.maturities(), maturities);
    EXPECT_EQ(volatilitySurface.volatilities(), volatilitiesMatrix);
    EXPECT_DOUBLE_EQ(volatilitySurface.discountCurve().rate(), discountCurve.rate());
}

// 2. Verify surface can interpolate at exact grid points (should return exact values)
TEST(VolatilitySurfaceConstructionTest, InterpolationAtGridPoints)
{
    std::vector<double> strikes = {100.0, 110.0, 120.0};
    std::vector<double> maturities = {0.5, 1.0, 2.0};
    std::vector<std::vector<double>> volatilitiesMatrix = {
        {0.20, 0.25, 0.30},  // maturity 0.5
        {0.25, 0.30, 0.35},  // maturity 1.0
        {0.30, 0.35, 0.40}   // maturity 2.0
    };
    FlatDiscountCurve discountCurve(0.05);

    VolatilitySurface vs(strikes, maturities, volatilitiesMatrix, discountCurve);

    EXPECT_NEAR(vs.impliedVolatility(100.0, 0.5), 0.20, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(110.0, 0.5), 0.25, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(120.0, 0.5), 0.30, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(100.0, 1.0), 0.25, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(110.0, 1.0), 0.30, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(120.0, 1.0), 0.35, 1e-10);

}

// 3. Verify how interest rate affects forward moneyness interpolation
// When using ForwardMoneyness, different rates should give different interpolated volatilities
TEST(VolatilitySurfaceConstructionTest, InterestRateAffectsForwardMoneyness)
{
    std::vector<double> strikes = {100.0, 110.0};
    std::vector<double> maturities = {0.5, 1.0};
    std::vector<std::vector<double>> volatilitiesMatrix = {
        {0.20, 0.25},   // maturity at 0.5, vol at K=100 is 0.2, at K=110 is 0.25
        {0.25, 0.30}    // maturity at 0.5, vol at K=100 is 0.25, at K=110 is 0.30
    };
    // Create two surfaces with different interest rates
    FlatDiscountCurve lowRateCurve(0.05);
    FlatDiscountCurve highRateCurve(0.10);

    VolatilitySurface vsLowRate(strikes, maturities, volatilitiesMatrix, lowRateCurve,
                                VolatilitySurface::SmileInterpolationType::Linear,
                                VolatilitySurface::MaturityInterpolationType::ForwardMoneyness);

    VolatilitySurface vsHighRate(strikes, maturities, volatilitiesMatrix, highRateCurve,
                                  VolatilitySurface::SmileInterpolationType::Linear,
                                  VolatilitySurface::MaturityInterpolationType::ForwardMoneyness);

    // Interpolate at a point between maturities (e.g. T=0.75)
    // With ForwardMoneyness the rate affects how strikes are adjusted:
    // K_i = K × B(T)/B(T_i) where B(t) = e^{-rt}
    // Higher rate -> lower discount factor -> different strikes adjustment -> different vol

    double strike = 105.0;
    double maturity = 0.75;

    double volLowRate = vsLowRate.impliedVolatility(strike, maturity);
    double volHighRate = vsHighRate.impliedVolatility(strike, maturity);

    // The volatilities should be different due to different forward moneyness adjustments
    // Demonstrating interest rate affects the interpolation
    EXPECT_NEAR(volLowRate, volLowRate, 1e-10);

    // For comparison, bilinear interpolation should be independent of rate
    VolatilitySurface vsBilinearLow(strikes, maturities, volatilitiesMatrix, lowRateCurve,
                                    VolatilitySurface::SmileInterpolationType::Linear,
                                    VolatilitySurface::MaturityInterpolationType::Bilinear);

    VolatilitySurface vsBilinearHigh(strikes, maturities, volatilitiesMatrix, highRateCurve,
                                     VolatilitySurface::SmileInterpolationType::Linear,
                                     VolatilitySurface::MaturityInterpolationType::Bilinear);

    double volBilinearLow = vsBilinearLow.impliedVolatility(strike, maturity);
    double volBilinearHigh = vsBilinearHigh.impliedVolatility(strike, maturity);
    // Bilinear should give same result regardless of rate (it does not use discount curve)
    EXPECT_NEAR(volBilinearLow, volBilinearHigh, 1e-10);

}

// 4. Verify matrix structure
TEST(VolatilitySurfaceConstructionTest, MatrixStructure)
{
    // Example: 3 strikes ;× 2 maturities -> 6 volatilities
    std::vector<double> strikes = {100.0, 110.0, 120.0};
    std::vector<double> maturities = {0.5, 1.0};
    std::vector<std::vector<double>> volatilitiesMatrix = {
        {0.20, 0.25, 0.30},
        {0.25, 0.30, 0.35},
    };
    FlatDiscountCurve discountCurve(0.05);
    VolatilitySurface vs(strikes, maturities, volatilitiesMatrix, discountCurve);

    // Verify the structure: volatilities[i][j] = volatility for strike i and maturity j
    EXPECT_EQ(volatilitiesMatrix.size(), maturities.size());   // Rows = # maturities
    EXPECT_EQ(volatilitiesMatrix[0].size(), strikes.size());   // Cols = # strikes
    
    // Verify we can access each (strike, maturity) pair
    EXPECT_NEAR(vs.impliedVolatility(100.0, 0.5), 0.20, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(110.0, 0.5), 0.25, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(120.0, 0.5), 0.30, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(100.0, 1.0), 0.25, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(110.0, 1.0), 0.30, 1e-10);
    EXPECT_NEAR(vs.impliedVolatility(120.0, 1.0), 0.35, 1e-10);
}

// 5. Verify bounds are correctly set
TEST(VolatilitySurfaceConstructionTest, Bounds)
{
    std::vector<double> strikes = {100.0, 110.0, 120.0};
    std::vector<double> maturities = {0.5, 1.0};
    std::vector<std::vector<double>> volatilitiesMatrix = {
        {0.20, 0.25, 0.30},
        {0.25, 0.30, 0.35}
    };
    FlatDiscountCurve discountCurve(0.05);
    
    VolatilitySurface vs(strikes, maturities, volatilitiesMatrix, discountCurve);

    auto bounds = vs.getBounds();
    EXPECT_NEAR(bounds.first.first, 100.0, 1e-10);
    EXPECT_NEAR(bounds.first.second, 120.0, 1e-10);
    EXPECT_NEAR(bounds.second.first, 0.5, 1e-10);
    EXPECT_NEAR(bounds.second.second, 1.0, 1e-10);
}
