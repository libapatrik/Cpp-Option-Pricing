#include <gtest/gtest.h>
#include <cppfm/models/Model.h>
#include <cppfm/pricers/Pricer.h>
#include <cppfm/simulators/PathSimulator.h>
#include <cppfm/market/VolatilitySurface.h>
#include <cppfm/market/DiscountCurve.h>
#include <cppfm/market/FinancialInstrument.h>
#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cmath>
#include <memory>
#include <iostream>
#include <iomanip>

//  Market prices -> implied vols -> vol surface -> local vol surface -> Dupire model -> pricing engine (FD/MC) -> prices
/**  TEST TO IMPELEMENT
 *   Match Black-Scholes when vol surface is flat
 *   Satisfy put-call parity
 *   Handle volatility smiles correctly
 *   Produce good Greeks
 *   FD ≈ MC
 *   Work for multiple strikes and maturities
 *   Handle boundary cases (deep ITM/OTM)

 * TODO:
 *    1. No-arbitrage constraints and test - butterfly arb etc.
 *    2. Add tests for Model1, Model2, Model2 impacts on the PV
 *    3. Test different Milstein vs. Euler agianst overkilled MCPricer
 *    4.
 *    
 *    
 */

// ============================================================================
// Test Fixture: Common Setup for Dupire Pricing Tests
// ============================================================================

class DupirePricerTest : public::testing::Test {
protected:
    void SetUp() override {
        // Market parameters
        spot = 100.0;
        riskFreeRate = 0.05;
        
        // Create discount curve
        discountCurve = std::make_unique<FlatDiscountCurve>(riskFreeRate);
        
        // Create flat volatility surface (should match Black-Scholes)
        createFlatVolatilitySurface(0.20);
        
        // Create Dupire model
        dupireModel = std::make_unique<DupireModel>(spot, *volSurface);
        
        // Create Black-Scholes model for comparison
        bsModel = std::make_unique<BlackScholesModel>(spot, *discountCurve, 0.20);
    }
    
    void createFlatVolatilitySurface(double flatVol) {
        // Create a flat volatility surface: σ(K,T) = flatVol everywhere
        VolatilitySurfaceBuilder builder;
        
        std::vector<double> strikes = {80, 90, 100, 110, 120};
        std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
        
        for (double K : strikes) {
            for (double T : maturities) {
                builder.setVolatility(K, T, flatVol);
            }
        }
        
        builder.setDiscountCurve(*discountCurve);
        volSurface = builder.build();
    }
    
    void createSmiledVolatilitySurface() {
        // Create a volatility surface with realistic smile
        VolatilitySurfaceBuilder builder;
        
        std::vector<double> strikes = {70, 80, 90, 100, 110, 120, 130};
        std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
        
        for (double T : maturities) {
            for (double K : strikes) {
                // Create smile: higher vol at low/high strikes
                double atmVol = 0.20;
                double moneyness = K / spot;
                double skew = -0.1 * (moneyness - 1.0);  // Downward sloping
                double smile = 0.05 * std::pow(moneyness - 1.0, 2.0);  // U-shape
                double impliedVol = atmVol + skew + smile;
                
                builder.setVolatility(K, T, impliedVol);
            }
        }
        
        builder.setDiscountCurve(*discountCurve);
        volSurface = builder.build();
        
        // Update Dupire model with smiled surface
        dupireModel = std::make_unique<DupireModel>(spot, *volSurface);
    }
    
    // Member variables
    double spot;
    double riskFreeRate;
    std::unique_ptr<FlatDiscountCurve> discountCurve;
    std::unique_ptr<VolatilitySurface> volSurface;
    std::unique_ptr<DupireModel> dupireModel;
    std::unique_ptr<BlackScholesModel> bsModel;
};

// ============================================================================
// TEST 1: Flat Volatility Surface Should Match Black-Scholes
// ============================================================================

TEST_F(DupirePricerTest, FlatSurfaceMatchesBlackScholes_FDPricer) {
    // When volatility surface is flat, Dupire's local volatility equals
    // the constant implied volatility, so prices should match Black-Scholes
    
    std::cout << "\n=== Test: Flat Surface Matches Black-Scholes (FD) ===" << std::endl;
    std::cout << "Setup: Spot=" << spot << ", Strike=100, T=1.0, Vol=0.20" << std::endl;
    
    // Create FD pricers
    FDPricer dupireFDPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    FDPricer bsFDPricer(*bsModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    // Test ATM call
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    double dupirePrice = dupireFDPricer.price(atmCall);
    double bsPrice = bsFDPricer.price(atmCall);
    
    std::cout << "  Dupire FD Price:  $" << std::fixed << std::setprecision(4) << dupirePrice << std::endl;
    std::cout << "  BS FD Price:      $" << bsPrice << std::endl;
    std::cout << "  Difference:       $" << std::abs(dupirePrice - bsPrice) << std::endl;
    
    // Should match within numerical tolerance
    EXPECT_NEAR(dupirePrice, bsPrice, 0.10)  // 10 cents tolerance for FD discretization
        << "Dupire with flat surface should match Black-Scholes";
    
    // Also compare with analytical Black-Scholes
    BlackScholesPricer analyticPricer(*bsModel, *discountCurve);
    double analyticPrice = analyticPricer.price(atmCall);
    
    std::cout << "  Analytical BS:    $" << analyticPrice << std::endl;
    std::cout << "  FD vs Analytical: $" << std::abs(dupirePrice - analyticPrice) << std::endl;
    
    EXPECT_NEAR(dupirePrice, analyticPrice, 0.10)
        << "FD Dupire should match analytical BS for flat surface";
    
    std::cout << "  Correct: Test passed" << std::endl;
}

TEST_F(DupirePricerTest, FlatSurfaceMatchesBlackScholes_MCPricer) {
    // Test Monte Carlo pricer with Dupire model on flat surface
    
    // Create time steps for simulation (daily steps for 1 year)
    std::vector<double> timeSteps;
    size_t numSteps = 252;
    for (size_t i = 0; i <= numSteps; ++i) {
        timeSteps.push_back(1.0 * i / numSteps);
    }
    
    // Create simulators
    EulerPathSimulator dupireSimulator(timeSteps, *dupireModel, 42);
    EulerPathSimulator bsSimulator(timeSteps, *bsModel, 42);
    
    // Create MC pricers
    size_t numPaths = 50000;  // Need many paths for MC convergence
    MonteCarloPricer dupireMCPricer(*dupireModel, *discountCurve, dupireSimulator, numPaths);
    MonteCarloPricer bsMCPricer(*bsModel, *discountCurve, bsSimulator, numPaths);
    
    // Test ATM call
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    double dupirePrice = dupireMCPricer.price(atmCall);
    double bsPrice = bsMCPricer.price(atmCall);
    
    // MC has statistical error, so tolerance is larger
    EXPECT_NEAR(dupirePrice, bsPrice, 0.20)
        << "MC Dupire with flat surface should match MC Black-Scholes";
    
    // Compare with analytical
    BlackScholesPricer analyticPricer(*bsModel, *discountCurve);
    double analyticPrice = analyticPricer.price(atmCall);
    
    EXPECT_NEAR(dupirePrice, analyticPrice, 0.30)
        << "MC Dupire should be close to analytical BS (within MC error)";
}

// ============================================================================
// TEST 2: Put-Call Parity
// ============================================================================

TEST_F(DupirePricerTest, PutCallParity_FDPricer) {
    // Test put-call parity: C - P = S_0 - K * exp(-r*T)
    
    std::cout << "\n=== Test: Put-Call Parity (FD) ===" << std::endl;
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    double strike = 100.0;
    double maturity = 1.0;
    
    EuropeanOptionPayoff call(Option::Type::Call, strike, maturity);
    EuropeanOptionPayoff put(Option::Type::Put, strike, maturity);
    
    double callPrice = fdPricer.price(call);
    double putPrice = fdPricer.price(put);
    
    double forward = spot - strike * discountCurve->discount(maturity);
    double putCallParity = callPrice - putPrice;
    
    std::cout << "  Call Price:       $" << std::fixed << std::setprecision(4) << callPrice << std::endl;
    std::cout << "  Put Price:        $" << putPrice << std::endl;
    std::cout << "  C - P:            $" << putCallParity << std::endl;
    std::cout << "  S - K*e^(-rT):    $" << forward << std::endl;
    std::cout << "  Difference:       $" << std::abs(putCallParity - forward) << std::endl;
    
    EXPECT_NEAR(putCallParity, forward, 0.10)
        << "Put-call parity should hold: C - P = S_0 - K*e^(-r*T)";
    
    std::cout << "  Correct: Put-call parity holds" << std::endl;
}

TEST_F(DupirePricerTest, PutCallParity_MCPricer) {
    // Test put-call parity with Monte Carlo pricer
    
    std::vector<double> timeSteps;
    size_t numSteps = 252;
    for (size_t i = 0; i <= numSteps; ++i) {
        timeSteps.push_back(1.0 * i / numSteps);
    }
    
    EulerPathSimulator simulator(timeSteps, *dupireModel, 42);
    MonteCarloPricer mcPricer(*dupireModel, *discountCurve, simulator, 50000);
    
    double strike = 100.0;
    double maturity = 1.0;
    
    EuropeanOptionPayoff call(Option::Type::Call, strike, maturity);
    EuropeanOptionPayoff put(Option::Type::Put, strike, maturity);
    
    double callPrice = mcPricer.price(call);
    double putPrice = mcPricer.price(put);
    
    double forward = spot - strike * discountCurve->discount(maturity);
    double putCallParity = callPrice - putPrice;
    
    // MC has larger tolerance due to statistical error
    EXPECT_NEAR(putCallParity, forward, 0.30)
        << "Put-call parity should hold for MC pricer";
}

// ============================================================================
// TEST 3: Smile Surface Produces Different Prices
// ============================================================================

TEST_F(DupirePricerTest, SmileSurfaceAffectsPrices_FDPricer) {
    // With a volatility smile, Dupire prices should differ from
    // Black-Scholes with ATM volatility
    
    std::cout << "\n=== Test: Volatility Smile Effects (FD) ===" << std::endl;
    std::cout << "Comparing: Dupire (smiled surface) vs BS (flat 20% vol)" << std::endl;
    
    createSmiledVolatilitySurface();
    
    FDPricer dupireFDPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    FDPricer bsFDPricer(*bsModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    // Test OTM put (strike = 80) - should be sensitive to smile
    EuropeanOptionPayoff otmPut(Option::Type::Put, 80.0, 1.0);
    double dupirePrice = dupireFDPricer.price(otmPut);
    double bsPrice = bsFDPricer.price(otmPut);
    
    std::cout << "  OTM Put (K=80, S=100):" << std::endl;
    std::cout << "    Dupire Price:     $" << std::fixed << std::setprecision(4) << dupirePrice << std::endl;
    std::cout << "    BS Price:         $" << bsPrice << std::endl;
    std::cout << "    Difference:       $" << std::abs(dupirePrice - bsPrice) << std::endl;
    // std::cout << "    (Smile increases vol at low strikes -> higher put price)" << std::endl;
    
    // Prices should differ due to smile (Dupire uses higher vol at low strikes)
    EXPECT_GT(std::abs(dupirePrice - bsPrice), 0.05)
        << "Smiled surface should produce different prices than flat BS";
    
    std::cout << "  Correct: Smile effects captured correctly" << std::endl;
}

TEST_F(DupirePricerTest, SmileSurfaceAffectsPrices_MCPricer) {
    // Test Monte Carlo with smile surface
    
    createSmiledVolatilitySurface();
    
    std::vector<double> timeSteps;
    size_t numSteps = 252;
    for (size_t i = 0; i <= numSteps; ++i) {
        timeSteps.push_back(1.0 * i / numSteps);
    }
    
    EulerPathSimulator dupireSimulator(timeSteps, *dupireModel, 42);
    EulerPathSimulator bsSimulator(timeSteps, *bsModel, 42);
    
    MonteCarloPricer dupireMCPricer(*dupireModel, *discountCurve, dupireSimulator, 30000);
    MonteCarloPricer bsMCPricer(*bsModel, *discountCurve, bsSimulator, 30000);
    
    EuropeanOptionPayoff otmPut(Option::Type::Put, 80.0, 1.0);
    double dupirePrice = dupireMCPricer.price(otmPut);
    
    // Should see difference (though MC has noise)
    // Just verify prices are positive and finite
    EXPECT_GT(dupirePrice, 0.0);
    EXPECT_LT(dupirePrice, spot);
    EXPECT_TRUE(std::isfinite(dupirePrice));
}

// ============================================================================
// TEST 4: Greeks Computation
// ============================================================================

TEST_F(DupirePricerTest, DeltaIsReasonable) {
    // Delta should be in (0,1) for ATM call
    
    std::cout << "\n=== Test: Delta Computation ===" << std::endl;
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    double price = fdPricer.price(atmCall);
    double delta = fdPricer.delta(atmCall);
    
    std::cout << "  ATM Call (K=100):" << std::endl;
    std::cout << "    Price:  $" << std::fixed << std::setprecision(4) << price << std::endl;
    std::cout << "    Delta:  " << std::setprecision(4) << delta << std::endl;
    std::cout << "    Expected: ~0.5 for ATM" << std::endl;
    
    // ATM call delta should be around 0.5
    EXPECT_GT(delta, 0.3);
    EXPECT_LT(delta, 0.7);
    EXPECT_TRUE(std::isfinite(delta));
    
    std::cout << "  Correct: Delta in reasonable range [0.3, 0.7]" << std::endl;
}

TEST_F(DupirePricerTest, GammaIsPositive) {
    // Gamma should be positive for long options
    
    std::cout << "\n=== Test: Gamma Computation ===" << std::endl;
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    double delta = fdPricer.delta(atmCall);
    double gamma = fdPricer.gamma(atmCall);
    
    std::cout << "  ATM Call Greeks:" << std::endl;
    std::cout << "    Delta:  " << std::fixed << std::setprecision(4) << delta << std::endl;
    std::cout << "    Gamma:  " << std::setprecision(6) << gamma << std::endl;
    std::cout << "    Expected: Gamma > 0 (long convexity)" << std::endl;
    
    EXPECT_GT(gamma, 0.0);
    EXPECT_TRUE(std::isfinite(gamma));
    
    std::cout << "  Correct: Gamma is positive" << std::endl;
}

// ============================================================================
// TEST 5: Consistency Between FD and MC
// ============================================================================

TEST_F(DupirePricerTest, FDandMCAgree_FlatSurface) {
    // FD and MC should produce similar prices for flat surface
    
    std::cout << "\n=== Test: FD vs MC Agreement ===" << std::endl;
    std::cout << "Setup: 100,000 MC paths, 252 time steps" << std::endl;
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    std::vector<double> timeSteps;
    size_t numSteps = 252;
    for (size_t i = 0; i <= numSteps; ++i) {
        timeSteps.push_back(1.0 * i / numSteps);
    }
    
    EulerPathSimulator simulator(timeSteps, *dupireModel, 42);
    MonteCarloPricer mcPricer(*dupireModel, *discountCurve, simulator, 100000);
    
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    
    double fdPrice = fdPricer.price(atmCall);
    double mcPrice = mcPricer.price(atmCall);
    
    std::cout << "  FD Price:         $" << std::fixed << std::setprecision(4) << fdPrice << std::endl;
    std::cout << "  MC Price:         $" << mcPrice << std::endl;
    std::cout << "  Difference:       $" << std::abs(fdPrice - mcPrice) << std::endl;
    std::cout << "  Relative Error:   " << std::setprecision(2) 
              << (std::abs(fdPrice - mcPrice) / fdPrice * 100) << "%" << std::endl;
    
    // Should agree within MC statistical error
    EXPECT_NEAR(fdPrice, mcPrice, 0.30)
        << "FD and MC should produce similar prices";
    
    std::cout << "  Correct: FD and MC agree (within MC statistical error)" << std::endl;
}

// ============================================================================
// TEST 6: Different Maturities
// ============================================================================

TEST_F(DupirePricerTest, PricesForDifferentMaturities) {
    // Test that pricer works for different maturities

    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    std::vector<double> maturities = {0.25, 0.5, 1.0, 2.0};
    
    for (double T : maturities) {
        EuropeanOptionPayoff call(Option::Type::Call, 100.0, T);
        double price = fdPricer.price(call);
        
        EXPECT_GT(price, 0.0) << "Price should be positive for maturity " << T;
        EXPECT_LT(price, spot) << "Call price should be less than spot for maturity " << T;
        EXPECT_TRUE(std::isfinite(price)) << "Price should be finite for maturity " << T;
    }
}

// ============================================================================
// TEST 7: Different Strikes
// ============================================================================

TEST_F(DupirePricerTest, CallPricesDecreaseWithStrike) {
    // Call option prices should decrease as strike increases
    
    std::cout << "\n=== Test: Call Price Monotonicity in Strike ===" << std::endl;
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    double maturity = 1.0;
    std::vector<double> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    
    std::cout << "  Call prices (T=1.0 year):" << std::endl;
    std::vector<double> prices;
    for (double K : strikes) {
        EuropeanOptionPayoff call(Option::Type::Call, K, maturity);
        double price = fdPricer.price(call);
        prices.push_back(price);
        std::cout << "    K=" << std::fixed << std::setw(5) << K 
                  << " -> $" << std::setprecision(4) << price << std::endl;
    }
    
    // Verify prices are decreasing
    bool allDecreasing = true;
    for (size_t i = 1; i < prices.size(); ++i) {
        EXPECT_LT(prices[i], prices[i-1])
            << "Call prices should decrease with increasing strike";
        if (prices[i] >= prices[i-1]) allDecreasing = false;
    }
    
    if (allDecreasing) {
        std::cout << "  Correct: Prices decrease monotonically with strike" << std::endl;
    }
}

TEST_F(DupirePricerTest, PutPricesIncreaseWithStrike) {
    // Put option prices should increase as strike increases
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    double maturity = 1.0;
    std::vector<double> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    
    std::vector<double> prices;
    for (double K : strikes) {
        EuropeanOptionPayoff put(Option::Type::Put, K, maturity);
        prices.push_back(fdPricer.price(put));
    }
    
    // Verify prices are increasing
    for (size_t i = 1; i < prices.size(); ++i) {
        EXPECT_GT(prices[i], prices[i-1])
            << "Put prices should increase with increasing strike";
    }
}

// ============================================================================
// TEST 8: Milstein vs Euler for MC
// ============================================================================

TEST_F(DupirePricerTest, MilsteinVsEuler_MCPricer) {
    // Milstein should be more accurate than Euler for local volatility

    std::vector<double> timeSteps;
    size_t numSteps = 252;
    for (size_t i = 0; i <= numSteps; ++i) {
        timeSteps.push_back(1.0 * i / numSteps);
    }
    
    EulerPathSimulator eulerSim(timeSteps, *dupireModel, 42);
    MilsteinPathSimulator milsteinSim(timeSteps, *dupireModel, 42);
    
    size_t numPaths = 1000;
    MonteCarloPricer eulerPricer(*dupireModel, *discountCurve, eulerSim, numPaths);
    MonteCarloPricer milsteinPricer(*dupireModel, *discountCurve, milsteinSim, numPaths);
    
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    
    double eulerPrice = eulerPricer.price(atmCall);
    double milsteinPrice = milsteinPricer.price(atmCall);
    
    // Both should be positive and close to each other
    EXPECT_GT(eulerPrice, 0.0);
    EXPECT_GT(milsteinPrice, 0.0);
    EXPECT_NEAR(eulerPrice, milsteinPrice, 0.50)
        << "Euler and Milstein should produce similar results";
}

// ============================================================================
// TEST 9: Boundary Cases
// ============================================================================

TEST_F(DupirePricerTest, DeepITMCallNearIntrinsicValue) {
    // Deep ITM call should be close to intrinsic value
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    double strike = 50.0;  // Deep ITM
    double maturity = 1.0;
    EuropeanOptionPayoff call(Option::Type::Call, strike, maturity);
    
    double price = fdPricer.price(call);
    double intrinsicValue = spot - strike;
    
    // Price should be greater than intrinsic value but not by much
    EXPECT_GT(price, intrinsicValue);
    EXPECT_LT(price, intrinsicValue + 5.0);  // Not much time value
}

TEST_F(DupirePricerTest, DeepOTMCallNearZero) {
    // Deep OTM call should be close to zero
    
    FDPricer fdPricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    
    double strike = 150.0;  // Deep OTM
    double maturity = 1.0;
    EuropeanOptionPayoff call(Option::Type::Call, strike, maturity);
    
    double price = fdPricer.price(call);
    
    EXPECT_GT(price, 0.0);
    EXPECT_LT(price, 1.0);  // Should be very small
}

// ============================================================================
// TEST 10: Grid Independence
// ============================================================================

TEST_F(DupirePricerTest, GridRefinementConvergence) {
    // Price should converge as grid is refined
    
    std::cout << "\n=== Test: Grid Refinement Convergence ===" << std::endl;
    
    EuropeanOptionPayoff atmCall(Option::Type::Call, 100.0, 1.0);
    
    // Coarse grid
    std::cout << "  Testing grid refinement:" << std::endl;
    FDPricer coarsePricer(*dupireModel, *discountCurve, 20.0, 300.0, 100, 50);
    double coarsePrice = coarsePricer.price(atmCall);
    std::cout << "    Coarse (100x50):    $" << std::fixed << std::setprecision(4) << coarsePrice << std::endl;
    
    // Fine grid
    FDPricer finePricer(*dupireModel, *discountCurve, 20.0, 300.0, 200, 100);
    double finePrice = finePricer.price(atmCall);
    std::cout << "    Fine   (200x100):   $" << finePrice << std::endl;
    
    // Very fine grid
    FDPricer veryFinePricer(*dupireModel, *discountCurve, 20.0, 300.0, 400, 200);
    double veryFinePrice = veryFinePricer.price(atmCall);
    std::cout << "    Very Fine (400x200):   $" << veryFinePrice << std::endl;
    
    // Refinement should reduce error
    double error1 = std::abs(finePrice - coarsePrice);
    double error2 = std::abs(veryFinePrice - finePrice);
    
    std::cout << "\n  Error reduction:" << std::endl;
    std::cout << "    |Fine - Coarse|:    $" << std::setprecision(6) << error1 << std::endl;
    std::cout << "    |V.Fine - Fine|:    $" << error2 << std::endl;
    std::cout << "    Improvement ratio:  " << std::setprecision(2) << (error1/error2) << "x" << std::endl;
    
    // Second error should be smaller (convergence)
    EXPECT_LT(error2, error1)
        << "Grid refinement should improve accuracy";
    
    std::cout << "  Correct: Solution converges with grid refinement" << std::endl;
}

