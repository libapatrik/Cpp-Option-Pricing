//
// Test file for Heston Model parameters and initialization
//

#include <gtest/gtest.h>
#include "../Model.h"
#include "../DiscountCurve.h"

/**
 * Test fixture for Heston Model tests
 */
class HestonModelTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Standard Heston model parameters
        S0 = 100.0;
        v0 = 0.04;      // Initial variance (σ₀² = 20%²)
        kappa = 2.0;    // Mean reversion speed
        theta = 0.04;   // Long-term variance mean
        sigma_v = 0.3;  // Volatility of volatility
        rho = -0.7;     // Correlation (typically negative for equity)
        r = 0.05;       // Risk-free rate

        // Create discount curve
        discountCurve = new FlatDiscountCurve(r);

        // Create Heston model
        hestonModel = new HestonModel(S0, *discountCurve, v0, kappa, theta, sigma_v, rho);
    }

    void TearDown() override
    {
        delete hestonModel;
        delete discountCurve;
    }

    // Model parameters
    double S0, v0, kappa, theta, sigma_v, rho, r;
    FlatDiscountCurve* discountCurve;
    HestonModel* hestonModel;
};

/**
 * Test Feller condition
 * Feller condition: 2κθ => σ²
 * If satisfied, variance process stays positive (in continuous limit)
 * TODO: EXERCISE: Investigate how Feller condition affects PV of options under Heston model
 *  - Model1: Feller violated vs. Model2: Feller satisfied
 *  - Model3: Feller ϵ close vs. Model4: Feller equality
 *  - Smaller κ higher vol-of-vol increase chance of hitting zero variance
 *  - TODO: Do cases for different schemes used for Heston path simulation
 */
TEST_F(HestonModelTest, FellerConditionTest)
{
    // Access parameters using individual accessors
    double kappa_val = hestonModel->kappa();
    double theta_val = hestonModel->vbar();  // vbar is θ (long-term variance mean)
    double sigma_v_val = hestonModel->sigma_v();
    
    double feller_value = 2.0 * kappa_val * theta_val;
    double sigma_squared = sigma_v_val * sigma_v_val;

    // Check if Feller condition is satisfied
    bool feller_satisfied = (feller_value >= sigma_squared);

    // For our default parameters: 2*2*0.04 = 0.16, σ_v² = 0.09
    // So Feller condition IS satisfied
    EXPECT_TRUE(feller_satisfied) << "Feller condition: 2κθ = " << feller_value
                                   << " should be ≥ σ_v² = " << sigma_squared;

    // Also test the built-in Feller condition checker
    EXPECT_TRUE(hestonModel->satisfiesFellerCondition());
}

/**
 * Test parameter validation
 * Validate the parameters
 * Impact of parameters on PV
 *  - Are some Heston schemes more robust? Take BK as reference.
 */
TEST_F(HestonModelTest, ParameterTest)
{
    // All parameters should be positive
    EXPECT_GT(S0, 0.0) << "Spot price must be positive";
    EXPECT_GT(v0, 0.0) << "Initial variance must be positive";
    EXPECT_GT(hestonModel->v0(), 0.0) << "Initial variance in model must be positive";
    EXPECT_GT(hestonModel->kappa(), 0.0) << "Mean reversion speed must be positive";
    EXPECT_GT(hestonModel->vbar(), 0.0) << "Long-term variance must be positive";
    EXPECT_GT(hestonModel->sigma_v(), 0.0) << "Vol of vol must be positive";

    // Correlation must be in [-1, 1]
    EXPECT_GE(hestonModel->correlation(), -1.0) << "Correlation must be >= -1";
    EXPECT_LE(hestonModel->correlation(), 1.0) << "Correlation must be <= 1";
    
    // Verify correlation is set correctly
    EXPECT_DOUBLE_EQ(hestonModel->rho(), rho);
    EXPECT_DOUBLE_EQ(hestonModel->correlation(), rho);
}

