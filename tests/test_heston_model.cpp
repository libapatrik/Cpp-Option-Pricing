#include <gtest/gtest.h>
#include <cppfm/models/Model.h>
#include <cppfm/market/DiscountCurve.h>
#include <cppfm/simulators/PathSimulator2D.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>



/** ? In practice the Feller condition is almost never satisfied, so how can we price this?
 * * Recall Grzelak: in lecture metnioned: "even if the Feller condition is not satisfied we use these methods" - but Euler breaks, so how exactly?
 *
 *
 *  OBSERVATIONS:
 *   1. Feller violated -> More zero-variance occurrences, especially with Euler
 *   2. BK Scheme handles all cases
 * ! 3. TG and QE perform well even when Feller is violated - interesting
 *   4. Euler scheme shows large biased when Feller is violated
 *   5. Higher  κ (speed of mean reversion) -> better behaved variance process.
 *
 *
 *   From the tests: QE seems robust wrt to Feller condition
 *   So, even if the Feller condition is not satisfied, we are able to provide good-enough approximate prices (relative to Broadie-Kaya)
 */






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
        // discountCurve = new FlatDiscountCurve(r);
        discountCurve = std::make_unique<FlatDiscountCurve>(r);
        // Create Heston model
        // hestonModel = new HestonModel(S0, *discountCurve, v0, kappa, theta, sigma_v, rho);
        hestonModel = std::make_unique<HestonModel>(S0, *discountCurve, v0, kappa, theta, sigma_v, rho);
    }

    // void TearDown() override
    // {
    //     delete hestonModel;
    //     delete discountCurve;
    // }

    // Model parameters
    double S0, v0, kappa, theta, sigma_v, rho, r;
    // FlatDiscountCurve* discountCurve;
    std::unique_ptr<HestonModel> hestonModel;
    std::unique_ptr<FlatDiscountCurve> discountCurve;
    // HestonModel* hestonModel;
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
 *  - NOTE: Grzelak mentions that  κ = 1/2 is standard (source videos)
 *
 *
 * 
 * Test Setup:    Fixed θ, κ tweaked
 * - Model1 (Violated):  2κθ = 0.08 < σ²_v = 0.36  (κ=1.0, θ=0.04, σ_v=0.6)
 * - Model2 (Satisfied): 2κθ = 0.32 > σ²_v = 0.09  (κ=4.0, θ=0.04, σ_v=0.3)
 * - Model3 (Near):      2κθ = 0.168 ≈ σ²_v = 0.16  (κ=2.1, θ=0.04, σ_v=0.4)
 * - Model4 (Equality):  2κθ = 0.16 = σ²_v = 0.16  (κ=2.0, θ=0.04, σ_v=0.4)
 * 
 * We compare prices across different discretization schemes:
 * - Euler: Known to have issues with negative variance
 * - BK: Exact scheme (reference)
 * - TG: Correct variance, no truncation
 * - QE: Quadratic exponential scheme
 */
 TEST_F(HestonModelTest, FellerConditionImpactOnOptionPricing)
 {
     const size_t numPaths = 10000;
     const size_t seed = 42;
     const double T = 1.0;  // 1 year maturity
     const double K = 100.0;  // ATM strike
     std::vector<double> timeSteps = {0.0, 0.25, 0.5, 0.75, 1.0};  // Euler discr error = Δt=0.25, errors expected to be large
     // When Feller is violated the Euler will need finer steps to reduce the negative variance occurences
     // std::vector<double> timeSteps = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; // Check for finer mesh
     // Noting, BK is exact.

     // Lambda to price a call option via Monte Carlo
     auto priceCall = [&](PathSimulator2D& simulator) -> std::tuple<double, int, int> {
         double payoffSum = 0.0;
         int zeroVarianceCount = 0;
         int pathsWithZero = 0; // unique 0 hits
         
         for (size_t i = 0; i < numPaths; ++i) {
             auto [assetPath, variancePath] = simulator.paths();
             bool pathHitZero = false;

             // Count how many times variance hits zero (or near-zero)
             for (double v : variancePath) {
                 if (v < 1e-10) {
                     zeroVarianceCount++;
                     pathHitZero = true;
                 }
             }
             if (pathHitZero) pathsWithZero++;
             
             // Call payoff at maturity
             double ST = assetPath.back();
             payoffSum += std::max(ST - K, 0.0);
         }
         
         double avgPayoff = payoffSum / numPaths;
         double callPrice = avgPayoff * std::exp(-r * T);  // Discount to present
         
         return {callPrice, zeroVarianceCount, pathsWithZero};
     };
     
     // Lambda to compute Feller value
     auto fellerInfo = [](double kappa, double theta, double sigma_v) {
         double lhs = 2.0 * kappa * theta;
         double rhs = sigma_v * sigma_v;
         double diff = lhs - rhs;
         return std::make_tuple(lhs, rhs, diff);
     };
     
     std::cout << "\n";
     std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
     std::cout << "FELLER CONDITION IMPACT ON OPTION PRICING\n";
     std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
     std::cout << "Setup: ATM Call, K=" << K << ", T=" << T << " year, " << numPaths << " paths\n";
     std::cout << "       S₀=" << S0 << ", v₀=0.04, r=" << r << ", ρ=" << rho << "\n\n";
     
     // =========================================================================
     // Model 1: Feller VIOLATED (high vol-of-vol, low mean reversion)
     // =========================================================================
     {
         double kappa1 = 1.0, theta1 = 0.04, sigma_v1 = 0.6;
         auto [lhs, rhs, diff] = fellerInfo(kappa1, theta1, sigma_v1);
         
         FlatDiscountCurve dc1(r);
         HestonModel model1(S0, dc1, v0, kappa1, theta1, sigma_v1, rho);
         
         std::cout << "MODEL 1: FELLER VIOLATED\n";
         std::cout << "  Parameters: κ=" << kappa1 << ", θ=" << theta1 << ", σᵥ=" << sigma_v1 << "\n";
         std::cout << "  Feller: 2κθ = " << std::setprecision(4) << lhs 
                   << " < σ²ᵥ = " << rhs << " (diff = " << std::showpos << diff << ")\n" << std::noshowpos;
         std::cout << "  Status: VIOLATED (variance may hit zero)\n\n";
         
         EulerPathSimulator2D simEuler1(timeSteps, model1, seed);
         BKExactPathSimulator2D simBK1(timeSteps, model1, seed, NewtonMethod::Optimized);
         TGPathSimulator2D simTG1(timeSteps, model1, seed);
         QEPathSimulator2D simQE1(timeSteps, model1, seed);
         
         auto [priceEuler1, zerosEuler1, pathsZeroEuler1] = priceCall(simEuler1);
         auto [priceBK1, zerosBK1, pathsZeroBK1] = priceCall(simBK1);
         auto [priceTG1, zerosTG1, pathsZeroTG1] = priceCall(simTG1);
         auto [priceQE1, zerosQE1, pathsZeroQE1] = priceCall(simQE1);
         
         std::cout << "  Scheme    │ Call Price │ Zero-Var Hits │ Paths w/Zero │ Bias vs BK\n";
         std::cout << "  ──────────┼────────────┼───────────────┼──────────────┼────────────\n";
         std::cout << "  Euler     │ " << std::setw(10) << std::setprecision(4) << priceEuler1 
                   << " │ " << std::setw(13) << zerosEuler1
                   << " │ " << std::setw(12) << pathsZeroEuler1
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceEuler1 - priceBK1) / priceBK1 * 100 << "%\n" << std::noshowpos;
         std::cout << "  BK *ref*  │ " << std::setw(10) << std::setprecision(4) << priceBK1 
                   << " │ " << std::setw(13) << zerosBK1
                   << " │ " << std::setw(12) << pathsZeroBK1 << " │     -\n";
         std::cout << "  TG        │ " << std::setw(10) << std::setprecision(4) << priceTG1
                   << " │ " << std::setw(13) << zerosTG1
                   << " │ " << std::setw(12) << pathsZeroTG1
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceTG1 - priceBK1) / priceBK1 * 100 << "%\n" << std::noshowpos;
         std::cout << "  QE        │ " << std::setw(10) << std::setprecision(4) << priceQE1 
                   << " │ " << std::setw(13) << zerosQE1
                   << " │ " << std::setw(12) << pathsZeroQE1
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceQE1 - priceBK1) / priceBK1 * 100 << "%\n\n" << std::noshowpos;

         
         EXPECT_FALSE(model1.satisfiesFellerCondition());
     }
     
     // =========================================================================
     // Model 2: Feller SATISFIED (strong mean reversion, moderate vol-of-vol)
     // =========================================================================
     {
         double kappa2 = 4.0, theta2 = 0.04, sigma_v2 = 0.3;
         auto [lhs, rhs, diff] = fellerInfo(kappa2, theta2, sigma_v2);
         
         FlatDiscountCurve dc2(r);
         HestonModel model2(S0, dc2, v0, kappa2, theta2, sigma_v2, rho);
         
         std::cout << "MODEL 2: FELLER SATISFIED\n";
         std::cout << "  Parameters: κ=" << kappa2 << ", θ=" << theta2 << ", σᵥ=" << sigma_v2 << "\n";
         std::cout << "  Feller: 2κθ = " << std::setprecision(4) << lhs 
                   << " > σ²ᵥ = " << rhs << " (diff = " << std::showpos << diff << ")\n" << std::noshowpos;
         std::cout << "  Status: SATISFIED (variance stays positive)\n\n";
         
         EulerPathSimulator2D simEuler2(timeSteps, model2, seed);
         BKExactPathSimulator2D simBK2(timeSteps, model2, seed, NewtonMethod::Optimized);
         TGPathSimulator2D simTG2(timeSteps, model2, seed);
         QEPathSimulator2D simQE2(timeSteps, model2, seed);
         
         auto [priceEuler2, zerosEuler2, pathsZeroEuler2] = priceCall(simEuler2);
         auto [priceBK2, zerosBK2, pathsZeroBK2] = priceCall(simBK2);
         auto [priceTG2, zerosTG2, pathsZeroTG2] = priceCall(simTG2);
         auto [priceQE2, zerosQE2, pathsZeroQE2] = priceCall(simQE2);
         
         std::cout << "  Scheme    │ Call Price │ Zero-Var Hits │ Paths w/Zero │ Bias vs BK\n";
         std::cout << "  ──────────┼────────────┼───────────────┼──────────────┼────────────\n";
         std::cout << "  Euler     │ " << std::setw(10) << std::setprecision(4) << priceEuler2 
                   << " │ " << std::setw(13) << zerosEuler2
                   << " │ " << std::setw(9) << pathsZeroEuler2
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceEuler2 - priceBK2) / priceBK2 * 100 << "%\n" << std::noshowpos;
         std::cout << "  BK *ref*  │ " << std::setw(10) << std::setprecision(4) << priceBK2 
                   << " │ " << std::setw(13) << zerosBK2
                   << " │ " << std::setw(9) << pathsZeroBK2 << " │     -\n";
         std::cout << "  TG        │ " << std::setw(10) << std::setprecision(4) << priceTG2 
                   << " │ " << std::setw(13) << zerosTG2
                   << " │ " << std::setw(9) << pathsZeroTG2
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceTG2 - priceBK2) / priceBK2 * 100 << "%\n" << std::noshowpos;
         std::cout << "  QE        │ " << std::setw(10) << std::setprecision(4) << priceQE2 
                   << " │ " << std::setw(13) << zerosQE2
                   << " │ " << std::setw(9) << pathsZeroQE2
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceQE2 - priceBK2) / priceBK2 * 100 << "%\n\n" << std::noshowpos;
         
         EXPECT_TRUE(model2.satisfiesFellerCondition());
     }
     
     // =========================================================================
     // Model 3: Feller NEAR BOUNDARY (slightly above)
     // =========================================================================
     {
         double kappa3 = 2.1, theta3 = 0.04, sigma_v3 = 0.4;
         auto [lhs, rhs, diff] = fellerInfo(kappa3, theta3, sigma_v3);
         
         FlatDiscountCurve dc3(r);
         HestonModel model3(S0, dc3, v0, kappa3, theta3, sigma_v3, rho);
         
         std::cout << "MODEL 3: FELLER NEAR BOUNDARY (ε > 0)\n";
         std::cout << "  Parameters: κ=" << kappa3 << ", θ=" << theta3 << ", σᵥ=" << sigma_v3 << "\n";
         std::cout << "  Feller: 2κθ = " << std::setprecision(4) << lhs 
                   << " ≈ σ²ᵥ = " << rhs << " (diff = " << std::showpos << diff << ")\n" << std::noshowpos;
         std::cout << "  Status: barely SATISFIED (may have numerical issues)\n\n";
         
         EulerPathSimulator2D simEuler3(timeSteps, model3, seed);
         BKExactPathSimulator2D simBK3(timeSteps, model3, seed, NewtonMethod::Optimized);
         TGPathSimulator2D simTG3(timeSteps, model3, seed);
         QEPathSimulator2D simQE3(timeSteps, model3, seed);
         
         auto [priceEuler3, zerosEuler3, pathsZeroEuler3] = priceCall(simEuler3);
         auto [priceBK3, zerosBK3, pathsZeroBK3] = priceCall(simBK3);
         auto [priceTG3, zerosTG3, pathsZeroTG3] = priceCall(simTG3);
         auto [priceQE3, zerosQE3, pathsZeroQE3] = priceCall(simQE3);
         
         std::cout << "  Scheme    │ Call Price │ Zero-Var Hits │ Paths w/Zero │ Bias vs BK\n";
         std::cout << "  ──────────┼────────────┼───────────────┼──────────────┼────────────\n";
         std::cout << "  Euler     │ " << std::setw(10) << std::setprecision(4) << priceEuler3 
                   << " │ " << std::setw(13) << zerosEuler3
                   << " │ " << std::setw(12) << pathsZeroEuler3
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceEuler3 - priceBK3) / priceBK3 * 100 << "%\n" << std::noshowpos;
         std::cout << "  BK *ref*  │ " << std::setw(10) << std::setprecision(4) << priceBK3 
                   << " │ " << std::setw(13) << zerosBK3
                   << " │ " << std::setw(12) << pathsZeroBK3  << " │     -\n";
         std::cout << "  TG        │ " << std::setw(10) << std::setprecision(4) << priceTG3 
                   << " │ " << std::setw(13) << zerosTG3
                   << " │ " << std::setw(12) << pathsZeroTG3
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceTG3 - priceBK3) / priceBK3 * 100 << "%\n" << std::noshowpos;
         std::cout << "  QE        │ " << std::setw(10) << std::setprecision(4) << priceQE3 
                   << " │ " << std::setw(13) << zerosQE3
                   << " │ " << std::setw(12) << pathsZeroQE3
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceQE3 - priceBK3) / priceBK3 * 100 << "%\n\n" << std::noshowpos;
         
         EXPECT_TRUE(model3.satisfiesFellerCondition());
     }
     
     // =========================================================================
     // Model 4: Feller EQUALITY (exactly at boundary)
     // =========================================================================
     {
         double kappa4 = 2.0, theta4 = 0.04, sigma_v4 = 0.4;
         auto [lhs, rhs, diff] = fellerInfo(kappa4, theta4, sigma_v4);
         
         FlatDiscountCurve dc4(r);
         HestonModel model4(S0, dc4, v0, kappa4, theta4, sigma_v4, rho);
         
         std::cout << "MODEL 4: FELLER EQUALITY (2κθ = σ²ᵥ)\n";
         std::cout << "  Parameters: κ=" << kappa4 << ", θ=" << theta4 << ", σᵥ=" << sigma_v4 << "\n";
         std::cout << "  Feller: 2κθ = " << std::setprecision(4) << lhs 
                   << " = σ²ᵥ = " << rhs << " (diff = " << std::showpos << diff << ")\n" << std::noshowpos;
         std::cout << "  Status:  BOUNDARY (variance can touch zero but won't cross)\n\n";
         
         EulerPathSimulator2D simEuler4(timeSteps, model4, seed);
         BKExactPathSimulator2D simBK4(timeSteps, model4, seed, NewtonMethod::Optimized);
         TGPathSimulator2D simTG4(timeSteps, model4, seed);
         QEPathSimulator2D simQE4(timeSteps, model4, seed);
         
         auto [priceEuler4, zerosEuler4, pathsZeroEuler4] = priceCall(simEuler4);
         auto [priceBK4, zerosBK4, pathsZeroBK4] = priceCall(simBK4);
         auto [priceTG4, zerosTG4, pathsZeroTG4] = priceCall(simTG4);
         auto [priceQE4, zerosQE4, pathsZeroQE4] = priceCall(simQE4);
         
         std::cout << "  Scheme    │ Call Price │ Zero-Var Hits │ Paths w/Zero │ Bias vs BK\n";
         std::cout << "  ──────────┼────────────┼───────────────┼──────────────┼────────────\n";
         std::cout << "  Euler     │ " << std::setw(10) << std::setprecision(4) << priceEuler4 
                   << " │ " << std::setw(13) << zerosEuler4
                   << " │ " << std::setw(9) << pathsZeroEuler4
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceEuler4 - priceBK4) / priceBK4 * 100 << "%\n" << std::noshowpos;
         std::cout << "  BK *ref*  │ " << std::setw(10) << std::setprecision(4) << priceBK4 
                   << " │ " << std::setw(13) << zerosBK4
                   << " │ " << std::setw(9) << pathsZeroBK4 << " │     -\n";
         std::cout << "  TG        │ " << std::setw(10) << std::setprecision(4) << priceTG4 
                   << " │ " << std::setw(13) << zerosTG4
                   << " │ " << std::setw(9) << pathsZeroTG4
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceTG4 - priceBK4) / priceBK4 * 100 << "%\n" << std::noshowpos;
         std::cout << "  QE        │ " << std::setw(10) << std::setprecision(4) << priceQE4 
                   << " │ " << std::setw(13) << zerosQE4
                   << " │ " << std::setw(9) << pathsZeroQE4
                   << " │ " << std::setprecision(2) << std::showpos 
                   << (priceQE4 - priceBK4) / priceBK4 * 100 << "%\n\n" << std::noshowpos;
         
         EXPECT_TRUE(model4.satisfiesFellerCondition());  // Boundary case: >= satisfied
     }
 }







/**
 * Test parameter validation
 * Validate the parameters
 * Impact of parameters on PV
 *  - Are some Heston schemes more robust?
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

/// TODO: Impact of params on IV