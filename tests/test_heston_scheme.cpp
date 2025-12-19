//
// Test file for PathSimulator2D discretization schemes
//

#include <gtest/gtest.h>
#include "../PathSimulator2D.h"
#include "../Model.h"
#include "../DiscountCurve.h"
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>


/**
 * TEST: PathSimulator2D Discretization Schemes
 * ADD: Compare the implemented methods
 * ADD: How these methods work for long timeSteps, non-equidistant timeSteps?
 *           - BK should still be the most accurate - take as reference
 *
 * For pairs (X, V) of discretization schemes used for X, V processes:
 * - (Euler, Euler) - full truncation σ²_vV(t)Δt version
 * - (BK, BK) - Reference method: AES Heston simulation scheme
 * - (Eq33, TG) - QE scheme with correct variance, no truncation required
 * - (Eq33, QE) - QE scheme with correct variance, no truncation required


 * TOFIX: Add interest rate to the methods
 * 
 * NOTES: 
 * -
 * - Why is TG somewhat similar to Euler?
 *   - Euler has full truncation σ²_vV(t)Δt
 *   - TG has correct conditional variance , no truncation required
 *   - QE correct variance, no truncation required
 * TODO: Add description of Eq11, Eq33, BKTG, BKQE - what do these mean?
 */
class PathSimulator2DTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Standard Heston model parameters
        S0 = 100.0;     // Spot price
        v0 = 0.04;      // Initial variance -> std_dev = sqrt(v0) = 0.2 = "usual sigma" = volatility
        kappa = 2.0;    // Mean reversion speed
        vbar = 0.04;    // Long-term variance mean
        sigma_v = 0.3;  // Volatility of volatility
        rho = -0.7;     // Correlation (typically negative for equity)
        r = 0.05;       // Risk-free rate

        // Create discount curve
        FlatDiscountCurve discountCurve(r);

        // Create Heston model
        hestonModel = new HestonModel(S0, discountCurve, v0, kappa, vbar, sigma_v, rho);

        // Time steps: [0, 0.25, 0.5, 0.75, 1.0] (quarterly steps over 1 year)
        timeSteps = {0.0, 0.25, 0.5, 0.75, 1.0};  // equidistant time steps
        // timeSteps = {0.0, 0.1, 0.33, 0.5, 0.7,0.88, 0.95, 1.0}; // non-equidistant time steps
    }

    void TearDown() override
    {
        delete hestonModel;
    }

    // Model parameters
    double S0, v0, kappa, vbar, sigma_v, rho, r;
    HestonModel* hestonModel;
    std::vector<double> timeSteps;
};

TEST_F(PathSimulator2DTest, WhatRiskFreeRate)
{
    std::cout << "Risk-free rate: " << hestonModel->riskFreeRate() << "\n";
}

/**
 * TEST: Simple verification that risk-free rate affects prices
 */
 TEST_F(PathSimulator2DTest, RiskFreeRateTest)
 {   // Andersen uses no drift, we want to use drift - because we want to price the options.
     const size_t numPaths = 1000;
     auto seed = 1;

     // Create two models: r=0 and r=0.05
     FlatDiscountCurve zeroCurve(0.0);
     FlatDiscountCurve normalCurve(0.05);

     HestonModel modelZero(S0, zeroCurve, v0, kappa, vbar, sigma_v, rho);
     HestonModel modelNormal(S0, normalCurve, v0, kappa, vbar, sigma_v, rho);

    // Simulate with BK scheme
    BKExactPathSimulator2D simZero(timeSteps, modelZero, seed);
    BKExactPathSimulator2D simNormal(timeSteps, modelNormal, seed);

     std::vector<double> pricesZero, pricesNormal;

     for (size_t i = 0; i < numPaths; ++i) {
         auto [pZ, vZ] = simZero.paths();
         pricesZero.push_back(pZ.back());

         auto [pN, vN] = simNormal.paths();
         pricesNormal.push_back(pN.back());
     }

     // Compute means
     auto mean = [](const std::vector<double>& v) {
         double sum = 0.0;
         for (double x : v) sum += x;
         return sum / v.size();
     };

     double meanZero = mean(pricesZero);
     double meanNormal = mean(pricesNormal);
     double actualRatio = meanNormal / meanZero;
     double expectedRatio = std::exp(0.05 * 1.0);  // e^(r*T)

     std::cout << "\n";
     std::cout << "Simple Risk-Free Rate Test (T=1 year, " << numPaths << " paths)\n";
     std::cout << "─────────────────────────────────────────────────────\n";
     std::cout << "Mean S(T) with r=0.00:  " << std::setprecision(5) << std::fixed << meanZero << "\n";
     std::cout << "Mean S(T) with r=0.05:  " << meanNormal << "\n";
     std::cout << "Actual ratio:           " << std::setprecision(5) << actualRatio << "\n";
     std::cout << "Expected ratio e^(rT):  " << expectedRatio << "\n";
     std::cout << "Difference:             " << std::setprecision(5)
               << std::abs(actualRatio - expectedRatio) / expectedRatio * 100 << "%\n";
     std::cout << "─────────────────────────────────────────────────────\n";
     std::cout << "Risk-free rate is used correctly!\n\n";

     // Verify: Mean with r=0.05 should be higher than with r=0
     EXPECT_GT(meanNormal, meanZero);

     // Verify: Ratio should be close to e^(rT) = 1.0513
     EXPECT_NEAR(actualRatio, expectedRatio, 0.02);
 }



/**
 * Test Broadie-Kaya (BK) - Original Newton vs. Optimized Newton method
 * Compares performance and correctness of two Newton implementations for CDF inversion:
 * - Original: Recomputes Fourier coefficients in each iteration
 * - Optimized: Caches coefficients for faster evaluation
 */
TEST_F(PathSimulator2DTest, BKNewtonTest)
{
    size_t seed = 1;
    const size_t numPaths = 100;
    
    // ========================================================================
    // Test Original Newton Method
    // ========================================================================
    BKExactPathSimulator2D simulatorOriginal(timeSteps, *hestonModel, seed, NewtonMethod::Original);
    
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> finalPricesOriginal;
    std::vector<double> finalVariancesOriginal;
    
    for (size_t i = 0; i < numPaths; ++i) {
        auto [assetPath, variancePath] = simulatorOriginal.paths();
        finalPricesOriginal.push_back(assetPath.back());
        finalVariancesOriginal.push_back(variancePath.back());
        
        // Sanity checks
        EXPECT_EQ(assetPath.size(), timeSteps.size());
        EXPECT_EQ(variancePath.size(), timeSteps.size());
        EXPECT_DOUBLE_EQ(assetPath[0], S0);
        EXPECT_DOUBLE_EQ(variancePath[0], v0);
        
        for (const auto& price : assetPath) {
            EXPECT_TRUE(std::isfinite(price));
            EXPECT_GT(price, 0.0);
        }
        for (const auto& var : variancePath) {
            EXPECT_TRUE(std::isfinite(var));
            EXPECT_GE(var, 0.0);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto durationOriginal = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // ========================================================================
    // Test Optimized Newton Method
    // ========================================================================
    BKExactPathSimulator2D simulatorOptimized(timeSteps, *hestonModel, seed, NewtonMethod::Optimized);
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> finalPricesOptimized;
    std::vector<double> finalVariancesOptimized;
    
    for (size_t i = 0; i < numPaths; ++i) {
        auto [assetPath, variancePath] = simulatorOptimized.paths();
        finalPricesOptimized.push_back(assetPath.back());
        finalVariancesOptimized.push_back(variancePath.back());
        
        // Sanity checks
        EXPECT_EQ(assetPath.size(), timeSteps.size());
        EXPECT_EQ(variancePath.size(), timeSteps.size());
        EXPECT_DOUBLE_EQ(assetPath[0], S0);
        EXPECT_DOUBLE_EQ(variancePath[0], v0);
        
        for (const auto& price : assetPath) {
            EXPECT_TRUE(std::isfinite(price));
            EXPECT_GT(price, 0.0);
        }
        for (const auto& var : variancePath) {
            EXPECT_TRUE(std::isfinite(var));
            EXPECT_GE(var, 0.0);
        }
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto durationOptimized = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    auto computeStats = [](const std::vector<double>& data) {
        double sum = 0.0;
        for (double x : data) sum += x;
        double mean = sum / data.size();
        
        double variance = 0.0;
        for (double x : data) variance += (x - mean) * (x - mean);
        double stddev = std::sqrt(variance / data.size());
        
        return std::make_pair(mean, stddev);
    };
    
    auto [meanPriceOriginal, stdPriceOriginal] = computeStats(finalPricesOriginal);
    auto [meanPriceOptimized, stdPriceOptimized] = computeStats(finalPricesOptimized);
    auto [meanVarOriginal, stdVarOriginal] = computeStats(finalVariancesOriginal);
    auto [meanVarOptimized, stdVarOptimized] = computeStats(finalVariancesOptimized);

    std::cout << "\n";
    std::cout << "Broadie-Kaya Scheme: 2 Newton Methods Comparison\n";

    std::cout << "    Number of paths:  " << std::setw(4) << numPaths << "\n";
    std::cout << "    Time steps:       " << std::setw(4) << timeSteps.size() << "\n";
    std::cout << "    Maturity:         " << std::setw(4) << std::setprecision(2) << std::fixed << timeSteps.back() << " years\n";
    std::cout << "\n";
    std::cout << "    Method           │  Time (ms) │   Final S(T) ± SE    │   Final V(T) ± SE      \n";
    std::cout << "    Original Newton  │ " << std::setw(9) << durationOriginal.count()
              << "   │ " << std::setw(6) << std::setprecision(4) << meanPriceOriginal
              << " ± " << std::setw(4) << std::setprecision(4) << stdPriceOriginal/sqrt(numPaths)
              << "   │ " << std::setw(6) << std::setprecision(4) << meanVarOriginal
              << " ± " << std::setw(5) << std::setprecision(4) << stdVarOriginal/sqrt(numPaths) << "\n";
    std::cout << "    Optimized Newton │ " << std::setw(9) << durationOptimized.count()
              << "   │ " << std::setw(6) << std::setprecision(4) << meanPriceOptimized
              << " ± " << std::setw(4) << std::setprecision(4) << stdPriceOptimized/sqrt(numPaths)
              << "   │ " << std::setw(6) << std::setprecision(4) << meanVarOptimized
              << " ± " << std::setw(5) << std::setprecision(4) << stdVarOptimized/sqrt(numPaths) << "\n";

    std::cout << std::endl;
    
    // ========================================================================
    // Verify that results are statistically similar (same seed -> same results)
    // ========================================================================
    EXPECT_NEAR(meanPriceOriginal, meanPriceOptimized, 0.01) 
        << "Mean prices should be nearly identical for same seed";
    EXPECT_NEAR(meanVarOriginal, meanVarOptimized, 0.001) 
        << "Mean variances should be nearly identical for same seed";
}

/**
 *  Table comparing the methods simulating the pairs (X, V)
 *  (Euler, Euler), (BK, BK), (Eq33, TG), (Eq33, QE)
 *
 */
 TEST_F(PathSimulator2DTest, CompareHestonScheme)
 {
     auto seed = 1;
     const size_t numPaths = 50000;
 
    EulerPathSimulator2D simulatorEuler(timeSteps, *hestonModel, seed);
    BKExactPathSimulator2D simulatorBK(timeSteps, *hestonModel, seed, NewtonMethod::Optimized);
    TGPathSimulator2D simulatorTG(timeSteps, *hestonModel, seed);
    QEPathSimulator2D simulatorQE(timeSteps, *hestonModel, seed);
    BKTGPathSimulator2D simulatorBKTG(timeSteps, *hestonModel, seed, 0.5, 0.5);
    BKQEPathSimulator2D simulatorBKQE(timeSteps, *hestonModel, seed, 1.5, 0.5, 0.5);

    std::vector<double> finalPricesEuler, finalVariancesEuler;
    std::vector<double> finalPricesBK, finalVariancesBK;
    std::vector<double> finalPricesTG, finalVariancesTG;
    std::vector<double> finalPricesQE, finalVariancesQE;
    std::vector<double> finalPricesBKTG, finalVariancesBKTG;
    std::vector<double> finalPricesBKQE, finalVariancesBKQE;
 
     // Euler
     auto startEuler = std::chrono::high_resolution_clock::now();
     for (size_t i = 0; i < numPaths; ++i) {
         auto [assetPath, variancePath] = simulatorEuler.paths();
         finalPricesEuler.push_back(assetPath.back());
         finalVariancesEuler.push_back(variancePath.back());
 
         // Sanity checks
         EXPECT_EQ(assetPath.size(), timeSteps.size());
         EXPECT_EQ(variancePath.size(), timeSteps.size());
         EXPECT_DOUBLE_EQ(assetPath[0], S0);
         EXPECT_DOUBLE_EQ(variancePath[0], v0);
         
         for (const auto& price : assetPath) {
             EXPECT_TRUE(std::isfinite(price));
             EXPECT_GT(price, 0.0);
         }
         for (const auto& var : variancePath) {
             EXPECT_TRUE(std::isfinite(var));
             EXPECT_GE(var, 0.0);
         }
     }
     auto endEuler = std::chrono::high_resolution_clock::now();
     auto durationEuler = std::chrono::duration_cast<std::chrono::milliseconds>(endEuler - startEuler);

     // BK
     auto startBK = std::chrono::high_resolution_clock::now();
     for (size_t i = 0; i < numPaths; ++i) {
         auto [assetPath, variancePath] = simulatorBK.paths();
         finalPricesBK.push_back(assetPath.back());
         finalVariancesBK.push_back(variancePath.back());
     }
     auto endBK = std::chrono::high_resolution_clock::now();
     auto durationBK = std::chrono::duration_cast<std::chrono::milliseconds>(endBK - startBK);

     // TG
     auto startTG = std::chrono::high_resolution_clock::now();
     for (size_t i = 0; i < numPaths; ++i) {
         auto [assetPath, variancePath] = simulatorTG.paths();
         finalPricesTG.push_back(assetPath.back());
         finalVariancesTG.push_back(variancePath.back());
     }
     auto endTG = std::chrono::high_resolution_clock::now();
     auto durationTG = std::chrono::duration_cast<std::chrono::milliseconds>(endTG - startTG);

    // QE
    auto startQE = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numPaths; ++i) {
        auto [assetPath, variancePath] = simulatorQE.paths();
        finalPricesQE.push_back(assetPath.back());
        finalVariancesQE.push_back(variancePath.back());
    }
    auto endQE = std::chrono::high_resolution_clock::now();
    auto durationQE = std::chrono::duration_cast<std::chrono::milliseconds>(endQE - startQE);

    // BKTG
    auto startBKTG = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numPaths; ++i) {
        auto [assetPath, variancePath] = simulatorBKTG.paths();
        finalPricesBKTG.push_back(assetPath.back());
        finalVariancesBKTG.push_back(variancePath.back());
    }
    auto endBKTG = std::chrono::high_resolution_clock::now();
    auto durationBKTG = std::chrono::duration_cast<std::chrono::milliseconds>(endBKTG - startBKTG);

    // BKQE
    auto startBKQE = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < numPaths; ++i) {
        auto [assetPath, variancePath] = simulatorBKQE.paths();
        finalPricesBKQE.push_back(assetPath.back());
        finalVariancesBKQE.push_back(variancePath.back());
    }
    auto endBKQE = std::chrono::high_resolution_clock::now();
    auto durationBKQE = std::chrono::duration_cast<std::chrono::milliseconds>(endBKQE - startBKQE);

     // Compute statistics: mean, std dev, and standard error
     auto computeStats = [](const std::vector<double>& data) {
         double sum = 0.0;
         for (double x : data) sum += x;
         double mean = sum / data.size();
         
         double variance = 0.0;
         for (double x : data) variance += (x - mean) * (x - mean);
         double stddev = std::sqrt(variance / data.size());
         double stderr = stddev / std::sqrt(data.size());  // Standard Error of Mean MC
         
         return std::make_tuple(mean, stddev, stderr);
     };
 
     auto [meanPriceEuler, stdPriceEuler, sePriceEuler] = computeStats(finalPricesEuler);
     auto [meanVarEuler, stdVarEuler, seVarEuler] = computeStats(finalVariancesEuler);
 
     auto [meanPriceBK, stdPriceBK, sePriceBK] = computeStats(finalPricesBK);
     auto [meanVarBK, stdVarBK, seVarBK] = computeStats(finalVariancesBK);
 
     auto [meanPriceTG, stdPriceTG, sePriceTG] = computeStats(finalPricesTG);
     auto [meanVarTG, stdVarTG, seVarTG] = computeStats(finalVariancesTG);
 
    auto [meanPriceQE, stdPriceQE, sePriceQE] = computeStats(finalPricesQE);
    auto [meanVarQE, stdVarQE, seVarQE] = computeStats(finalVariancesQE);

    auto [meanPriceBKTG, stdPriceBKTG, sePriceBKTG] = computeStats(finalPricesBKTG);
    auto [meanVarBKTG, stdVarBKTG, seVarBKTG] = computeStats(finalVariancesBKTG);

    auto [meanPriceBKQE, stdPriceBKQE, sePriceBKQE] = computeStats(finalPricesBKQE);
    auto [meanVarBKQE, stdVarBKQE, seVarBKQE] = computeStats(finalVariancesBKQE);

    // Calculate bias vs. BK (reference method)
    double biasEulerPrice = (meanPriceEuler - meanPriceBK) / meanPriceBK * 100.0;
    double biasEulerVar = (meanVarEuler - meanVarBK) / meanVarBK * 100.0;
    double biasTGPrice = (meanPriceTG - meanPriceBK) / meanPriceBK * 100.0;
    double biasTGVar = (meanVarTG - meanVarBK) / meanVarBK * 100.0;
    double biasQEPrice = (meanPriceQE - meanPriceBK) / meanPriceBK * 100.0;
    double biasQEVar = (meanVarQE - meanVarBK) / meanVarBK * 100.0;
    double biasBKTGPrice = (meanPriceBKTG - meanPriceBK) / meanPriceBK * 100.0;
    double biasBKTGVar = (meanVarBKTG - meanVarBK) / meanVarBK * 100.0;
    double biasBKQEPrice = (meanPriceBKQE - meanPriceBK) / meanPriceBK * 100.0;
    double biasBKQEVar = (meanVarBKQE - meanVarBK) / meanVarBK * 100.0;

    // Store times as doubles
    double timeEulerMs = durationEuler.count();
    double timeBKMs = durationBK.count();
    double timeTGMs = durationTG.count();
    double timeQEMs = durationQE.count();
    double timeBKTGMs = durationBKTG.count();
    double timeBKQEMs = durationBKQE.count();
     

     std::cout << "\n";
     std::cout << "═══════════════════════════════════════════════════════════════════════════════════════\n";
     std::cout << "Heston Discretization Schemes Comparison\n";
     std::cout << "═══════════════════════════════════════════════════════════════════════════════════════\n";
     std::cout << "Simulation Parameters:\n";
     std::cout << "  • Paths:    " << numPaths << "\n";
     std::cout << "  • Steps:    " << timeSteps.size() << "\n";
     std::cout << "  • Maturity: " << std::setprecision(2) << std::fixed << timeSteps.back() << " year\n";
     std::cout << "  • Model:    S₀=" << S0 << ", v₀=" << v0 << ", κ=" << kappa 
               << ", v̄=" << vbar << ", σᵥ=" << sigma_v << ", ρ=" << rho << ", r=" << r << "\n";
     std::cout << "───────────────────────────────────────────────────────────────────────────────────────\n\n";
 
     std::cout << "ASSET PRICE S(T) AT MATURITY:\n";
     std::cout << "─────────────────────────────────────────────────────────────────────────────────────────────────\n";
     std::cout << "Method         │ Time (ms) │  Mean ± SE      │ Std Dev │ Bias vs BK (%) │ Speedup vs BK\n";
     std::cout << "───────────────┼───────────┼─────────────────┼─────────┼────────────┼──────────────\n";
 
     std::cout << std::fixed;
     std::cout << "(Euler, Euler) │ " << std::setw(8) << std::setprecision(1) << timeEulerMs
               << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceEuler
               << " ± " << std::setw(4) << std::setprecision(4) << sePriceEuler
               << " │ " << std::setw(6) << std::setprecision(4) << stdPriceEuler
               << "  │ " << std::setw(7) << std::setprecision(4) << std::showpos << biasEulerPrice << "%"
               << std::noshowpos << " │ " << std::setw(10) << std::setprecision(0) << (timeBKMs / timeEulerMs) << "x\n";
 
     std::cout << "(BKExact) *ref* │ " << std::setw(8) << std::setprecision(1) << timeBKMs
               << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceBK
               << " ± " << std::setw(4) << std::setprecision(4) << sePriceBK
               << " │ " << std::setw(6) << std::setprecision(4) << stdPriceBK
               << "  │     -      │       1x\n";
 
     std::cout << "(Eq33, TG)     │ " << std::setw(8) << std::setprecision(1) << timeTGMs
               << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceTG
               << " ± " << std::setw(4) << std::setprecision(4) << sePriceTG
               << " │ " << std::setw(6) << std::setprecision(4) << stdPriceTG
               << "  │ " << std::setw(7) << std::setprecision(4) << std::showpos << biasTGPrice << "%"
               << std::noshowpos << " │ " << std::setw(10) << std::setprecision(0) << (timeBKMs / timeTGMs) << "x\n";
 
    std::cout << "(Eq33, QE)     │ " << std::setw(8) << std::setprecision(1) << timeQEMs
              << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceQE
              << " ± " << std::setw(4) << std::setprecision(4) << sePriceQE
              << " │ " << std::setw(6) << std::setprecision(4) << stdPriceQE
              << "  │ " << std::setw(7) << std::setprecision(4) << std::showpos << biasQEPrice << "%"
              << std::noshowpos << " │ " << std::setw(10) << std::setprecision(0) << (timeBKMs / timeQEMs) << "x\n";

    std::cout << "(Eq11, BKTG)   │ " << std::setw(8) << std::setprecision(1) << timeBKTGMs
              << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceBKTG
              << " ± " << std::setw(4) << std::setprecision(4) << sePriceBKTG
              << " │ " << std::setw(6) << std::setprecision(4) << stdPriceBKTG
              << "  │ " << std::setw(7) << std::setprecision(4) << std::showpos << biasBKTGPrice << "%"
              << std::noshowpos << " │ " << std::setw(10) << std::setprecision(0) << (timeBKMs / timeBKTGMs) << "x\n";

    std::cout << "(Eq11, BKQE)   │ " << std::setw(8) << std::setprecision(1) << timeBKQEMs
              << "  │ " << std::setw(7) << std::setprecision(4) << meanPriceBKQE
              << " ± " << std::setw(4) << std::setprecision(4) << sePriceBKQE
              << " │ " << std::setw(6) << std::setprecision(4) << stdPriceBKQE
              << "  │ " << std::setw(7) << std::setprecision(4) << std::showpos << biasBKQEPrice << "%"
              << std::noshowpos << " │ " << std::setw(10) << std::setprecision(0) << (timeBKMs / timeBKQEMs) << "x\n";

    std::cout << "\n";
     std::cout << "VARIANCE V(T) AT MATURITY:\n";
     std::cout << "─────────────────────────────────────────────────────────────────────────────────────\n";
     std::cout << "Method         │  Mean ± SE        │ Std Dev │ Bias vs BK\n";
     std::cout << "───────────────┼───────────────────┼─────────┼────────────\n";

     std::cout << "(Euler, Euler) │ " << std::setw(7) << std::setprecision(4) << meanVarEuler
               << " ± " << std::setw(6) << std::setprecision(4) << seVarEuler
               << " │ " << std::setw(6) << std::setprecision(4) << stdVarEuler
               << "  │ " << std::setw(7) << std::setprecision(3) << std::showpos << biasEulerVar << "%\n" << std::noshowpos;

     std::cout << "(BKExact) *ref* │ " << std::setw(7) << std::setprecision(4) << meanVarBK
               << " ± " << std::setw(6) << std::setprecision(4) << seVarBK
               << " │ " << std::setw(6) << std::setprecision(4) << stdVarBK
               << "  │     -\n";

     std::cout << "(Eq33, TG)     │ " << std::setw(7) << std::setprecision(4) << meanVarTG
               << " ± " << std::setw(6) << std::setprecision(4) << seVarTG
               << " │ " << std::setw(6) << std::setprecision(4) << stdVarTG
               << "  │ " << std::setw(7) << std::setprecision(3) << std::showpos << biasTGVar << "%\n" << std::noshowpos;

    std::cout << "(Eq33, QE)     │ " << std::setw(7) << std::setprecision(4) << meanVarQE
              << " ± " << std::setw(6) << std::setprecision(4) << seVarQE
              << " │ " << std::setw(6) << std::setprecision(4) << stdVarQE
              << "  │ " << std::setw(7) << std::setprecision(3) << std::showpos << biasQEVar << "%\n" << std::noshowpos;

    std::cout << "(Eq11, BKTG)   │ " << std::setw(7) << std::setprecision(4) << meanVarBKTG
              << " ± " << std::setw(6) << std::setprecision(4) << seVarBKTG
              << " │ " << std::setw(6) << std::setprecision(4) << stdVarBKTG
              << "  │ " << std::setw(7) << std::setprecision(3) << std::showpos << biasBKTGVar << "%\n" << std::noshowpos;

    std::cout << "(Eq11, BKQE)   │ " << std::setw(7) << std::setprecision(4) << meanVarBKQE
              << " ± " << std::setw(6) << std::setprecision(4) << seVarBKQE
              << " │ " << std::setw(6) << std::setprecision(4) << stdVarBKQE
              << "  │ " << std::setw(7) << std::setprecision(3) << std::showpos << biasBKQEVar << "%\n" << std::noshowpos;

}


/**
 * Std Dev BK is the true volatility - methods should converge to this value as Δt -> 0
 * We would expect:
 *   Variance Bias:
 * - Euler:  Large negative bias (-25% to -30%)   Known issue
 * - TG:     Small positive bias (+5% to +10%)    Better
 * - QE:     Very small bias (±1%)                Best
 * - BK:     Reference (exact)                    Perfect
*/



/**
 * TEST: Discretization Convergence Study
 * Demonstrates that all schemes converge to BK as time step → 0
 */
TEST_F(PathSimulator2DTest, DiscretizationConvergence)
{
    auto seed = 1;
    const size_t numPaths = 1000;  // More paths for better MC accuracy
    
    // Test different time step sizes
    std::vector<int> numSteps = {4, 12, 24, 48};  // Quarterly, monthly, bi-weekly, weekly
    
    std::cout << "\n";
    std::cout << "Discretization Convergence Study: Bias vs. BK as Δt → 0\n";
    std::cout << "    Number of paths:  " << numPaths << "\n";
    std::cout << "    Maturity:         1.0 year\n";
    std::cout << "\n";
    std::cout << "    Steps │    Δt    │  Price Bias vs. BK (abs)  │  Variance Bias vs. BK (%)  \n";
    std::cout << "          │  (days)  │  Euler    TG      QE    │  Euler    TG      QE       \n";
    std::cout << "    ──────┼──────────┼─────────────────────────┼────────────────────────────\n";
    
    for (int n : numSteps) {
        // Create time steps: [0, Δt, 2Δt, ..., T]
        std::vector<double> steps;
        double T = 1.0;  // 1 year
        double dt = T / n;
        for (int i = 0; i <= n; ++i) {
            steps.push_back(i * dt);
        }
        
        // Create simulators with current time step
        EulerPathSimulator2D simEuler(steps, *hestonModel, seed);
        BKExactPathSimulator2D simBK(steps, *hestonModel, seed, NewtonMethod::Optimized);
        TGPathSimulator2D simTG(steps, *hestonModel, seed);
        QEPathSimulator2D simQE(steps, *hestonModel, seed);
        
        // Collect final values
        std::vector<double> pricesEuler, pricesBK, pricesTG, pricesQE;
        std::vector<double> varsEuler, varsBK, varsTG, varsQE;
        
        for (size_t i = 0; i < numPaths; ++i) {
            auto [pE, vE] = simEuler.paths();
            pricesEuler.push_back(pE.back());
            varsEuler.push_back(vE.back());
            
            auto [pB, vB] = simBK.paths();
            pricesBK.push_back(pB.back());
            varsBK.push_back(vB.back());
            
            auto [pT, vT] = simTG.paths();
            pricesTG.push_back(pT.back());
            varsTG.push_back(vT.back());
            
            auto [pQ, vQ] = simQE.paths();
            pricesQE.push_back(pQ.back());
            varsQE.push_back(vQ.back());
        }
        
        // Compute means
        auto mean = [](const std::vector<double>& v) {
            double sum = 0.0;
            for (double x : v) sum += x;
            return sum / v.size();
        };
        
        double meanPriceEuler = mean(pricesEuler);
        double meanPriceBK = mean(pricesBK);
        double meanPriceTG = mean(pricesTG);
        double meanPriceQE = mean(pricesQE);
        
        double meanVarEuler = mean(varsEuler);
        double meanVarBK = mean(varsBK);
        double meanVarTG = mean(varsTG);
        double meanVarQE = mean(varsQE);

        // Calculate biases in price and variance as absolutes
        double biasEulerPrice = (meanPriceEuler - meanPriceBK);
        double biasTGPrice = (meanPriceTG - meanPriceBK);
        double biasQEPrice = (meanPriceQE - meanPriceBK);

        // double biasEulerVar = (meanVarEuler - meanVarBK);
        // double biasTGVar = (meanVarTG - meanVarBK);
        // double biasQEVar = (meanVarQE - meanVarBK);

        // Calculate biases in price & variances as ratios
        // double biasEulerPrice = (meanPriceEuler - meanPriceBK) / meanPriceBK * 100.0;
        // double biasTGPrice = (meanPriceTG - meanPriceBK) / meanPriceBK * 100.0;
        // double biasQEPrice = (meanPriceQE - meanPriceBK) / meanPriceBK * 100.0;
        //
        double biasEulerVar = (meanVarEuler - meanVarBK) / meanVarBK * 100.0;
        double biasTGVar = (meanVarTG - meanVarBK) / meanVarBK * 100.0;
        double biasQEVar = (meanVarQE - meanVarBK) / meanVarBK * 100.0;
        
        // Print row
        double days = dt * 252;  // Convert to trading days
        std::cout << "    " << std::setw(5) << n 
                  << " │ " << std::setw(6) << std::setprecision(1) << std::fixed << days
                  << " │ " << std::setw(6) << std::setprecision(2) << std::showpos << biasEulerPrice
                  << "" << std::setw(6) << biasTGPrice
                  << "" << std::setw(6) << biasQEPrice << ""
                  << " │ " << std::setw(6) << biasEulerVar
                  << "" << std::setw(6) << biasTGVar
                  << "" << std::setw(6) << biasQEVar << "%\n" << std::noshowpos;
    }
    
    std::cout << "\n";
    std::cout << "    Expected: All biases tend to 0% as steps increase (Δt -> 0)\n";
    std::cout << std::endl;
}




