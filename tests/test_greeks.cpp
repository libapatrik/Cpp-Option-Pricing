#include <gtest/gtest.h>
#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cppfm/market/FinancialInstrument.h>
#include <cppfm/market/DiscountCurve.h>
#include <cmath>

/*
═══════════════════════════════════════════════════════════════════════════════
                    GREEKS: The Option Sensitivities
═══════════════════════════════════════════════════════════════════════════════
First-Order (Price Sensitivities):
  Δ (Delta)   -> How option price changes with spot price (∂V/∂S)
  ν (Vega)    -> How option price changes with volatility (∂V/∂σ)
  Θ (Theta)   -> How option price changes with time (∂V/∂t)
  ρ (Rho)     -> How option price changes with interest rate (∂V/∂r)

Second-Order (Convexity & Cross-effects):
  Γ (Gamma)   -> Rate of change of delta (∂2V/∂S2)
  Vanna       -> How delta changes with volatility (∂2V/∂S∂σ)
  Volga       -> How vega changes with volatility (∂2V/∂σ2)

TEST COVERAGE:
  - Basic properties and value ranges
  - Behavior across moneyness (ITM, ATM, OTM)
  - Boundary conditions (at maturity)
  - Put-call parity relationships
  - Numerical verification (finite differences)
  - Extreme parameter regimes
  - Time-dependent interest rates (discount curves)

═══════════════════════════════════════════════════════════════════════════════
    NOTES:  
        
*/

// ═══════════════════════════════════════════════════════════════════════════════
//                          TEST CONFIGURATION
// ═══════════════════════════════════════════════════════════════════════════════

// Test tolerances - balancing accuracy with numerical stability
constexpr double PRICE_TOL = 1e-6;         // For price comparisons
constexpr double GREEK_TOL = 1e-4;         // For analytical Greeks
constexpr double NUMERICAL_TOL = 1e-3;     // Looser for finite differences

// Standard market parameters (realistic mid-market conditions)
constexpr double S0 = 100.0;      // Spot price: $100
constexpr double K = 100.0;       // Strike: $100 (at-the-money)
constexpr double r = 0.05;        // Risk-free rate: 5% per annum
constexpr double sigma = 0.2;     // Volatility: 20% per annum
constexpr double T = 1.0;         // Time to maturity: 1 year

// Helper function for numerical derivatives
double numericalDerivative(
    std::function<double(double)> f,
    double x,
    double h = 1e-4)
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double numericalSecondDerivative(std::function<double(double)> f, double x, double h = 1e-4)
{
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

// ═══════════════════════════════════════════════════════════════════════════════
//                          DELTA TESTS (∂V/∂S)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS DELTA?
// Delta measures how much an option's price changes when the underlying asset
// price moves by $1. It represents the hedge ratio - how many shares of stock
// you need to buy/sell to neutralize the option position.
//
// FINANCIAL INTERPRETATION:
// - Call Delta ∈ [0, 1]:  Probability of finishing in-the-money (risk-neutral)
// - Put Delta ∈ [-1, 0]:  Negative because puts lose value as spot increases
// - Delta = 0.5:          Roughly at-the-money (adjusted for drift)
// - Delta -> 1:            Deep in-the-money (acts like stock)
// - Delta -> 0:            Deep out-of-the-money (almost worthless)
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, DeltaCallBasic) {
    // GOAL: Verify that call delta is within valid bounds and reasonable for ATM
    double delta = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Call);
    
    // Delta must be between 0 and 1 for calls (can't change more than the stock)
    EXPECT_GE(delta, 0.0) << "Call delta cannot be negative";
    EXPECT_LE(delta, 1.0) << "Call delta cannot exceed 1 (max leverage is 1:1)";
    
    // ATM call delta ≈ 0.5-0.65 (above 0.5 due to positive drift from r > 0)
    // 50% probability that the option will expiry ITM
    // The forward price F = S*e^(rT) > S, so ATM by strike is ITM by forward
    // F = S*e^(rT) = 100*e^(0.05*1) = 105.12710963718821 > S = 100
    EXPECT_NEAR(delta, 0.6, 0.15) 
        << "ATM call delta should be ~0.6 (slightly ITM on forward basis)";
}

TEST(GreeksTest, DeltaPutBasic) {
    // GOAL: Verify put delta is negative and within valid range
    double delta = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Put);
    
    // Put delta must be between -1 and 0 (inverse relationship with spot)
    EXPECT_GE(delta, -1.0) << "Put delta cannot be more negative than -1";
    EXPECT_LE(delta, 0.0) << "Put delta must be negative (puts lose value as S↑)";
    
    // ATM put delta ≈ -0.35 to -0.5 (less negative due to forward > spot)
    EXPECT_NEAR(delta, -0.4, 0.15)
        << "ATM put delta should be ~-0.4 (less negative than -0.5 due to drift)";
}

TEST(GreeksTest, DeltaPutCallParity) {
    // GOAL: Verify put-call parity holds for delta
    //
    // THEORY: From put-call parity: C - P = S - K*e^(-rT)
    // Taking ∂/∂S of both sides: ∂C/∂S - ∂P/∂S = 1
    // Therefore: Δ_Call - Δ_Put = 1 
    //
    // PRACTICAL MEANING: A synthetic long stock position (long call + short put)
    // has delta = 1, behaving exactly like the underlying stock.
    
    double deltaCall = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Call);
    double deltaPut = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Put);
    
    double expectedDiff = 1.0;
    EXPECT_NEAR(deltaCall - deltaPut, expectedDiff, GREEK_TOL)
        << "Put-call parity for delta: Long call + Short put = Long stock (Δ = 1)";
}

TEST(GreeksTest, DeltaMoneyness) {
    // GOAL: Test delta behavior across different moneyness levels
    //
    // MONEYNESS SPECTRUM:
    // ITM (K=80):  High probability of exercise -> Delta near 1.0
    // ATM (K=100): 50/50 probability -> Delta around 0.5-0.6
    // OTM (K=120): Low probability of exercise -> Delta near 0.0
    
    double deltaITM = BlackScholesFormulas::delta(S0, 80.0, r, sigma, T, Option::Type::Call);
    double deltaATM = BlackScholesFormulas::delta(S0, 100.0, r, sigma, T, Option::Type::Call);
    double deltaOTM = BlackScholesFormulas::delta(S0, 120.0, r, sigma, T, Option::Type::Call);
    
    // Delta should monotonically decrease as strike increases
    EXPECT_GT(deltaITM, deltaATM) 
        << "In-the-money delta should exceed at-the-money delta";
    EXPECT_GT(deltaATM, deltaOTM)
        << "At-the-money delta should exceed out-of-the-money delta";
    
    // Deep ITM -> acts like stock (delta ≈ 1)
    EXPECT_NEAR(deltaITM, 1.0, 0.1)
        << "Deep ITM call should behave like the underlying stock";
    
    // OTM still has time value with 1 year to expiry
    EXPECT_LT(deltaOTM, 0.4)
        << "OTM call delta should be < 0.4 (but still positive due to time value)";
}

TEST(GreeksTest, DeltaAtMaturity) {
    // GOAL: Test boundary condition - delta becomes deterministic at expiry
    //
    // AT EXPIRY (T=0): No time value remains, only intrinsic value
    // - ITM call: Definitely exercise -> Delta = 1 (certainty)
    // - OTM call: Definitely expires worthless -> Delta = 0 (no sensitivity)
    // - ATM call: Dirac delta at K (discontinuity)
    
    double T_zero = 0.0;
    
    // Spot > Strike: ITM call has delta = 1 (will exercise with certainty)
    double deltaITM = BlackScholesFormulas::delta(110.0, K, r, sigma, T_zero, Option::Type::Call);
    EXPECT_NEAR(deltaITM, 1.0, GREEK_TOL)
        << "At expiry, ITM call has Δ=1 (acts exactly like stock)";
    
    // Spot < Strike: OTM call has delta = 0 (worthless, no sensitivity)
    double deltaOTM = BlackScholesFormulas::delta(90.0, K, r, sigma, T_zero, Option::Type::Call);
    EXPECT_NEAR(deltaOTM, 0.0, GREEK_TOL)
        << "At expiry, OTM call has Δ=0 (no value, no sensitivity)";
}

TEST(GreeksTest, DeltaNumericalVerification) {
    // GOAL: Verify analytical delta formula against numerical differentiation
    //
    // METHOD: Compute ∂V/∂S using finite differences: [V(S+h) - V(S-h)] / 2h
    // This confirms our closed-form formula is correct (no algebra errors)
    
    auto callPrice = [](double spot) {
        return BlackScholesFormulas::callPrice(spot, K, r, sigma, T);
    };
    
    double analyticDelta = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Call);
    double numericalDelta = numericalDerivative(callPrice, S0);
    
    EXPECT_NEAR(analyticDelta, numericalDelta, NUMERICAL_TOL)
        << "Analytical delta (from formula) should match numerical derivative";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                          GAMMA TESTS (∂2V/∂S2)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS GAMMA?
// Gamma measures the rate of change of delta - i.e., how quickly your hedge needs
// to be adjusted as the underlying moves. It represents the convexity of the
// option price with respect to the spot price.
//
// FINANCIAL INTERPRETATION:
// - Γ > 0 (always):     Options have positive convexity (good!)
// - High Γ:             Delta changes rapidly -> need frequent rehedging
// - Maximum at ATM:     Highest uncertainty about exercise
// - Low Γ (deep ITM/OTM): Delta is stable -> easy to hedge
//
// WHY IT MATTERS:
// - Delta hedging P&L ≈ 1/2 * Γ * (ΔS)2 (gamma profit/loss)
// - Market makers earn theta, pay for rehedging costs (gamma risk)
// - "Long gamma" = profit from volatility (scalping)
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, GammaBasic) {
    // GOAL: Verify gamma is positive and within reasonable bounds
    //
    // Gamma is ALWAYS positive for long options (both calls and puts)
    // This is why options have positive convexity!
    
    double gamma = BlackScholesFormulas::gamma(S0, K, r, sigma, T);
    
    EXPECT_GT(gamma, 0.0)
        << "Gamma must be positive (options have positive convexity)";
    
    // Typical gamma for ATM 1Y option is 0.01-0.02 (1-2% per $1 move)
    EXPECT_LT(gamma, 0.1)
        << "Gamma should be in reasonable range (<0.1 for annual maturity)";
}

TEST(GreeksTest, GammaMoneyness) {
    // GOAL: Verify that gamma peaks at-the-money
    //
    // FINANCIAL INTUITION: ATM options have maximum uncertainty about exercise
    // - Small moves in S can swing the option from ITM to OTM or vice versa
    // - Therefore delta changes most rapidly near K -> maximum gamma
    // - Deep ITM/OTM: Delta is stable (near 1 or 0) -> low gamma
    
    double gammaITM = BlackScholesFormulas::gamma(S0, 80.0, r, sigma, T);
    double gammaATM = BlackScholesFormulas::gamma(S0, 100.0, r, sigma, T);
    double gammaOTM = BlackScholesFormulas::gamma(S0, 120.0, r, sigma, T);
    
    EXPECT_GT(gammaATM, gammaITM)
        << "ATM gamma should exceed ITM gamma (max uncertainty at-the-money)";
    EXPECT_GT(gammaATM, gammaOTM)
        << "ATM gamma should exceed OTM gamma (delta changes fastest at K)";
}

TEST(GreeksTest, GammaAtMaturity) {
    // GOAL: Test boundary condition - gamma vanishes at expiry
    //
    // AT EXPIRY: Delta becomes step function (0 if OTM, 1 if ITM)
    // Derivative of step function = Dirac delta (infinite spike at K)
    // In practice, gamma -> 0 everywhere except exactly at K
    // This represents: no more uncertainty, no need to rehedge
    
    double T_zero = 0.0;
    double gamma = BlackScholesFormulas::gamma(S0, K, r, sigma, T_zero);
    
    EXPECT_NEAR(gamma, 0.0, GREEK_TOL)
        << "At expiry, gamma -> 0 (delta no longer changes, hedging frozen)";
}

TEST(GreeksTest, GammaNumericalVerification) {
    // GOAL: Verify analytical gamma formula using second-order finite differences
    //
    // METHOD: Compute ∂2V/∂S^2 using: [V(S+h) - 2V(S) + V(S-h)] / h^2
    // This validates the closed-form expression for convexity
    
    auto callPrice = [](double spot) {
        return BlackScholesFormulas::callPrice(spot, K, r, sigma, T);
    };
    
    double analyticGamma = BlackScholesFormulas::gamma(S0, K, r, sigma, T);
    double numericalGamma = numericalSecondDerivative(callPrice, S0);
    
    EXPECT_NEAR(analyticGamma, numericalGamma, NUMERICAL_TOL)
        << "Analytical gamma should match numerical second derivative";
}

TEST(GreeksTest, GammaPutCallIdentity) {
    // GOAL: Verify that calls and puts have identical gamma
    //
    // THEORY: From put-call parity C - P = S - K*e^(-rT)
    // Taking ∂2/∂S^2 of both sides: ∂2C/∂S^2 - ∂2P/∂S^2 = 0
    // Therefore: Γ_Call = Γ_Put (exactly the same!)
    //
    // INTUITION: Both have same convexity profile around strike
    
    auto callPrice = [](double spot) {
        return BlackScholesFormulas::callPrice(spot, K, r, sigma, T);
    };
    auto putPrice = [](double spot) {
        return BlackScholesFormulas::putPrice(spot, K, r, sigma, T);
    };
    
    double gammaCall = numericalSecondDerivative(callPrice, S0);
    double gammaPut = numericalSecondDerivative(putPrice, S0);
    
    EXPECT_NEAR(gammaCall, gammaPut, NUMERICAL_TOL)
        << "Put and call gamma must be identical (from put-call parity)";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                          VEGA TESTS (∂V/∂σ)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS VEGA? 
// Vega measures how option value changes with volatility.
//
// FINANCIAL INTERPRETATION:
// - ν > 0 (always):     Higher volatility -> higher option value (both C and P)
// - "Long vol":         Buying options = long vega (profit from vol increase)
// - "Short vol":        Selling options = short vega (profit from vol decrease)
// - Maximum at ATM:     Most sensitive to changes in implied volatility
//
// MARKET APPLICATIONS:
// - Volatility trading: Express views on implied vs realized volatility
// - Vega hedging:       Balance options portfolio against vol moves
// - Vol arbitrage:      Trade mispriced implied volatility
// - Vega = S2*Γ*σ*T:    Direct relationship with gamma
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, VegaBasic) {
    // GOAL: Verify vega is positive and in reasonable range
    //
    // Vega must ALWAYS be positive for long options
    // More uncertainty = more optionality = higher value
    
    double vega = BlackScholesFormulas::vega(S0, K, r, sigma, T);
    
    EXPECT_GT(vega, 0.0)
        << "Vega must be positive (more volatility always increases option value)";
    
    // Typical vega for $100 stock with 1Y maturity is 20-40 (cents per vol point)
    EXPECT_LT(vega, 100.0)
        << "Vega should be in reasonable range for standard parameters";
}

TEST(GreeksTest, VegaMoneyness) {
    // GOAL: Verify vega peaks at-the-money
    //
    // INTUITION: ATM options have maximum sensitivity to volatility
    // - Deep ITM: Already in-the-money, additional vol has limited impact
    // - Deep OTM: Too far away, extra vol helps but limited value
    // - ATM: Right at the inflection point, vol moves matter most
    //
    // TRADING: "Sell ATM straddles" = collect maximum vega premium
    
    double vegaITM = BlackScholesFormulas::vega(S0, 80.0, r, sigma, T);
    double vegaATM = BlackScholesFormulas::vega(S0, 100.0, r, sigma, T);
    double vegaOTM = BlackScholesFormulas::vega(S0, 120.0, r, sigma, T);
    
    EXPECT_GT(vegaATM, vegaITM)
        << "ATM vega should exceed ITM (maximum vol sensitivity at-the-money)";
    EXPECT_GT(vegaATM, vegaOTM)
        << "ATM vega should exceed OTM (volatility matters most near strike)";
}

TEST(GreeksTest, VegaAtMaturity) {
    // GOAL: Test boundary condition - vega vanishes at expiry
    //
    // AT EXPIRY: Option value = max(S-K, 0) for calls (deterministic)
    // No time left -> volatility doesn't matter -> vega = 0
    // This makes intuitive sense: you know exactly what you'll get
    
    double T_zero = 0.0;
    double vega = BlackScholesFormulas::vega(S0, K, r, sigma, T_zero);
    
    EXPECT_NEAR(vega, 0.0, GREEK_TOL)
        << "At expiry, vega = 0 (no time left, vol doesn't matter)";
}

TEST(GreeksTest, VegaTimeEffect) {
    // GOAL: Verify vega increases with time to maturity
    //
    // INTUITION: More time = more opportunity for volatility to matter
    // - LEAPS (long-dated options) have high vega
    // - Vega ~ √T approximately (from diffusion)
    //
    // TRADING: "Calendar spread" exploits different vega across maturities
    
    double vega_short = BlackScholesFormulas::vega(S0, K, r, sigma, 0.25);  // 3 months
    double vega_long = BlackScholesFormulas::vega(S0, K, r, sigma, 1.0);    // 1 year
    
    EXPECT_GT(vega_long, vega_short)
        << "Longer maturity -> higher vega (more time for vol to impact price)";
}

TEST(GreeksTest, VegaNumericalVerification) {
    // GOAL: Verify analytical vega formula using finite differences
    //
    // METHOD: Bump volatility by small amount, measure price change
    // Vega = [V(σ+h) - V(σ-h)] / 2h
    
    auto callPrice = [](double vol) {
        return BlackScholesFormulas::callPrice(S0, K, r, vol, T);
    };
    
    double analyticVega = BlackScholesFormulas::vega(S0, K, r, sigma, T);
    double numericalVega = numericalDerivative(callPrice, sigma);
    
    EXPECT_NEAR(analyticVega, numericalVega, NUMERICAL_TOL)
        << "Analytical vega should match numerical derivative w.r.t. volatility";
}

TEST(GreeksTest, VegaPutCallIdentity) {
    // GOAL: Verify calls and puts have identical vega
    //
    // THEORY: From put-call parity, taking ∂/∂σ:
    // ∂C/∂σ = ∂P/∂σ (volatility affects both equally)
    //
    // INTUITION: Both have same exposure to volatility changes
    
    auto callPrice = [](double vol) {
        return BlackScholesFormulas::callPrice(S0, K, r, vol, T);
    };
    auto putPrice = [](double vol) {
        return BlackScholesFormulas::putPrice(S0, K, r, vol, T);
    };
    
    double vegaCall = numericalDerivative(callPrice, sigma);
    double vegaPut = numericalDerivative(putPrice, sigma);
    
    EXPECT_NEAR(vegaCall, vegaPut, NUMERICAL_TOL)
        << "Call and put vega must be identical (both benefit from vol equally)";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                          THETA TESTS (∂V/∂t or -∂V/∂T)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS THETA?
// Theta measures time decay - how option value changes as time passes.
// Convention: Often quoted as -∂V/∂T (loss per day) so it's typically negative.
//
// FINANCIAL INTERPRETATION:
// - Θ < 0 (usually):    Options lose value as time passes (time decay)
// - "Theta decay":      Weekend/holiday effect - passage of time erodes value
// - Maximum at ATM:     Most time value at-the-money -> most decay
// - Deep ITM/OTM:       Mostly intrinsic/worthless -> less time decay
//
// MARKET APPLICATIONS:
// - Option selling:     Collect theta by selling options (short volatility)
// - Calendar spreads:   Exploit different decay rates across maturities
// - "Time is money":    Theta is the price of optionality
// - Gamma-Theta tradeoff: Θ ≈ -1/2*σ2*S2*Γ (from Black-Scholes PDE)
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, ThetaCallBasic) {
    // GOAL: Verify call theta is negative (time decay)
    //
    // For ATM calls with r > 0, theta is always negative
    // Time decay reflects the erosion of extrinsic (time) value
    
    double theta = BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Call);
    
    // Most calls have negative theta (you pay for time)
    EXPECT_LT(theta, 0.0)
        << "Call theta should be negative (options lose value over time)";
}

TEST(GreeksTest, ThetaPutBasic) {
    // GOAL: Verify put theta is well-defined (can be positive or negative!)
    //
    // SURPRISING FACT: Deep ITM puts can have positive theta!
    // - Why? Present value effect: K*e^(-rT) -> K as T->0
    // - The "interest carry" can dominate time decay for ITM puts
    // - Most ATM/OTM puts have negative theta like calls
    
    double theta = BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Put);
    
    // Put theta's sign depends on moneyness and r
    EXPECT_TRUE(std::isfinite(theta))
        << "Put theta should be finite (sign depends on moneyness and rates)";
}

TEST(GreeksTest, ThetaAtMaturity) {
    // GOAL: Test boundary condition - theta vanishes at expiry
    //
    // AT EXPIRY: No time value left -> no decay -> theta = 0
    // The option value is purely intrinsic (if any)
    
    double T_zero = 0.0;
    double thetaCall = BlackScholesFormulas::theta(S0, K, r, sigma, T_zero, Option::Type::Call);
    double thetaPut = BlackScholesFormulas::theta(S0, K, r, sigma, T_zero, Option::Type::Put);
    
    EXPECT_NEAR(thetaCall, 0.0, GREEK_TOL)
        << "At expiry, theta = 0 (no time value to decay)";
    EXPECT_NEAR(thetaPut, 0.0, GREEK_TOL)
        << "At expiry, theta = 0 (no time value to decay)";
}

TEST(GreeksTest, ThetaMoneyness) {
    // GOAL: Verify ATM options have maximum time decay
    //
    // INTUITION: ATM options are "pure time value"
    // - ITM: Mostly intrinsic value (less to decay)
    // - OTM: Little total value (less absolute decay)
    // - ATM: Maximum extrinsic value -> maximum absolute theta
    //
    // TRADING: "Sell ATM options" maximizes theta collection
    
    double thetaITM = BlackScholesFormulas::theta(S0, 80.0, r, sigma, T, Option::Type::Call);
    double thetaATM = BlackScholesFormulas::theta(S0, 100.0, r, sigma, T, Option::Type::Call);
    double thetaOTM = BlackScholesFormulas::theta(S0, 120.0, r, sigma, T, Option::Type::Call);
    
    // ATM should have most negative theta (maximum time decay in absolute terms)
    EXPECT_LT(thetaATM, thetaITM)
        << "ATM theta should be most negative (maximum time value to decay)";
    EXPECT_LT(thetaATM, thetaOTM)
        << "ATM theta decays fastest (peak extrinsic value at-the-money)";
}

TEST(GreeksTest, ThetaNumericalVerification) {
    // GOAL: Verify analytical theta formula using finite differences
    //
    // CONVENTION: Theta = -∂V/∂T (negative of time derivative)
    // This makes theta negative for most options (intuitive: value decreases)
    
    auto callPrice = [](double time) {
        return BlackScholesFormulas::callPrice(S0, K, r, sigma, time);
    };
    
    double analyticTheta = BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Call);
    double numericalTheta = -numericalDerivative(callPrice, T);
    
    EXPECT_NEAR(analyticTheta, numericalTheta, NUMERICAL_TOL * 10)
        << "Analytical theta should match -∂V/∂T (looser tolerance due to singularity)";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                          RHO TESTS (∂V/∂r)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS RHO?
// Rho measures sensitivity to interest rate changes. Often the least important
// Greek in practice, but critical in high-rate environments.
//
// FINANCIAL INTERPRETATION:
// - ρ_Call > 0:         Calls benefit from higher rates (forward drift effect)
// - ρ_Put < 0:          Puts hurt by higher rates (present value of strike)
// - Increases with T:   Longer maturity -> more rate exposure
// - Important for Long-dated options sensitive to rate changes 
//
// WHY RATES MATTER:
// - Forward price F = S*e^(rT) increases with r -> helps calls
// - Present value K*e^(-rT) decreases with r -> hurts puts (strikes worth less)
// - Low rate environment (2010s): Rho was negligible
// - High rate environment (2022-2024): Rho matters again!
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, RhoCallBasic) {
    // GOAL: Verify call rho is positive
    //
    // INTUITION: Higher rates -> higher forward price -> calls worth more
    // The stock "grows" faster in the risk-neutral measure
    
    double rho = BlackScholesFormulas::rho(S0, K, r, sigma, T, Option::Type::Call);
    
    EXPECT_GT(rho, 0.0)
        << "Call rho should be positive (calls benefit from higher interest rates)";
}

TEST(GreeksTest, RhoPutBasic) {
    // GOAL: Verify put rho is negative
    //
    // INTUITION: Higher rates -> strike's present value decreases
    // The payoff K-S is worth less in present value terms -> puts worth less
    
    double rho = BlackScholesFormulas::rho(S0, K, r, sigma, T, Option::Type::Put);
    
    EXPECT_LT(rho, 0.0)
        << "Put rho should be negative (puts hurt by higher interest rates)";
}

TEST(GreeksTest, RhoAtMaturity) {
    // GOAL: Test boundary condition - rho vanishes at expiry
    //
    // AT EXPIRY: No time for rates to matter -> rho = 0
    // The payoff is immediate, no discounting or drift effects
    
    double T_zero = 0.0;
    double rhoCall = BlackScholesFormulas::rho(S0, K, r, sigma, T_zero, Option::Type::Call);
    double rhoPut = BlackScholesFormulas::rho(S0, K, r, sigma, T_zero, Option::Type::Put);
    
    EXPECT_NEAR(rhoCall, 0.0, GREEK_TOL)
        << "At expiry, rho = 0 (no time for rates to affect value)";
    EXPECT_NEAR(rhoPut, 0.0, GREEK_TOL)
        << "At expiry, rho = 0 (immediate payoff, no discounting)";
}

TEST(GreeksTest, RhoTimeEffect) {
    // GOAL: Verify rho increases with time to maturity
    //
    // INTUITION: Longer maturity -> more time for rates to compound
    // - Short-dated options: rho ≈ 0 (rates don't matter much)
    // - LEAPS (2+ years): significant rho exposure
    // - Rho ~ T approximately
    
    double rho_short = BlackScholesFormulas::rho(S0, K, r, sigma, 0.25, Option::Type::Call);  // 3M
    double rho_long = BlackScholesFormulas::rho(S0, K, r, sigma, 1.0, Option::Type::Call);    // 1Y
    
    EXPECT_GT(std::abs(rho_long), std::abs(rho_short))
        << "Longer maturity -> higher |rho| (more time for interest rates to compound)";
}

TEST(GreeksTest, RhoNumericalVerification) {
    // GOAL: Verify analytical rho formula using finite differences
    //
    // METHOD: Bump interest rate, measure price change
    // Rho = [V(r+h) - V(r-h)] / 2h
    
    auto callPrice = [](double rate) {
        return BlackScholesFormulas::callPrice(S0, K, rate, sigma, T);
    };
    
    double analyticRho = BlackScholesFormulas::rho(S0, K, r, sigma, T, Option::Type::Call);
    double numericalRho = numericalDerivative(callPrice, r);
    
    EXPECT_NEAR(analyticRho, numericalRho, NUMERICAL_TOL)
        << "Analytical rho should match numerical derivative w.r.t. interest rate";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                     VANNA TESTS (∂2V/∂S∂σ)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS VANNA?
// Vanna is a second-order cross-Greek measuring how delta changes with volatility,
// or equivalently, how vega changes with spot price. It's a "cross-gamma".
//
// FINANCIAL INTERPRETATION:
// - Vanna = ∂Δ/∂σ = ∂ν/∂S    Cross-derivative (symmetric by Schwarz's theorem)
// - Can be positive or negative depending on moneyness
// - Important for exotic option hedging
// - "Volga and vanna" = convexity risks in volatility trading
//
// WHY IT MATTERS:
// - Delta hedging errors: When vol changes, delta changes too (vanna effect)
// - Vega hedging errors: When spot moves, vega exposure changes (vanna again)
// - Critical for volatility arbitrage and exotic derivatives
// - Together with volga, captures second-order vol risk
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, VannaBasic) {
    // GOAL: Verify vanna is well-defined
    //
    // Vanna can be positive or negative depending on strike/spot relationship
    // It measures how correlated your delta and vega risks are
    
    double vanna = BlackScholesFormulas::vanna(S0, K, r, sigma, T);
    
    // Vanna can have either sign - just check it's finite and reasonable
    EXPECT_TRUE(std::isfinite(vanna))
        << "Vanna should be finite (sign depends on moneyness)";
}

TEST(GreeksTest, VannaAtMaturity) {
    // GOAL: Test boundary condition - vanna vanishes at expiry
    //
    // AT EXPIRY: Both delta and vega freeze -> no cross-sensitivity
    // Delta becomes step function, vega = 0
    
    double T_zero = 0.0;
    double vanna = BlackScholesFormulas::vanna(S0, K, r, sigma, T_zero);
    
    EXPECT_NEAR(vanna, 0.0, GREEK_TOL)
        << "At expiry, vanna = 0 (no time for cross-effects between spot and vol)";
}

TEST(GreeksTest, VannaPutCallIdentity) {
    // GOAL: Verify calls and puts have identical vanna
    //
    // THEORY: From put-call parity, taking ∂2/∂S∂σ:
    // Vanna_Call = Vanna_Put (cross-derivatives commute and match)
    //
    // METHOD: Compute ∂Δ/∂σ numerically for both
    
    auto callPrice = [](double spot, double vol) {
        return BlackScholesFormulas::callPrice(spot, K, r, vol, T);
    };
    auto putPrice = [](double spot, double vol) {
        return BlackScholesFormulas::putPrice(spot, K, r, vol, T);
    };
    
    // Numerical cross-derivative for call: ∂(∂C/∂S)/∂σ
    auto callDeltaFunc = [&](double vol) {
        auto priceFunc = [&](double spot) { return callPrice(spot, vol); };
        return numericalDerivative(priceFunc, S0);
    };
    double vannaCall = numericalDerivative(callDeltaFunc, sigma);
    
    // Numerical cross-derivative for put: ∂(∂P/∂S)/∂σ
    auto putDeltaFunc = [&](double vol) {
        auto priceFunc = [&](double spot) { return putPrice(spot, vol); };
        return numericalDerivative(priceFunc, S0);
    };
    double vannaPut = numericalDerivative(putDeltaFunc, sigma);
    
    EXPECT_NEAR(vannaCall, vannaPut, NUMERICAL_TOL * 10)
        << "Call and put vanna must be identical (from put-call parity)";
}

// ═══════════════════════════════════════════════════════════════════════════════
//                     VOLGA TESTS (∂2V/∂σ2)
// ═══════════════════════════════════════════════════════════════════════════════
//
// WHAT IS VOLGA?
// Volga (also called "vomma" or "volga") measures the convexity of option value
// with respect to volatility - i.e., how vega changes as vol changes.
//
// FINANCIAL INTERPRETATION:
// - Volga = ∂2V/∂σ2 = ∂ν/∂σ    Second derivative w.r.t. volatility
// - Always positive (like gamma but for vol)
// - Maximum for ATM options (typically)
// - "Convexity in volatility space"
//
// WHY IT MATTERS:
// - Volatility trading: Profit from changes in implied vol
// - Vega hedging errors: When vol moves, vega changes (volga effect)
// - Variance swaps: Pure volga exposure
// - VIX options: High volga means explosive moves
// - Together with vanna, captures all second-order vol risk
//
// MARKET APPLICATIONS:
// - Long straddles have positive volga (convex in vol)
// - Benefit from volatility of volatility (vol-of-vol)
// - Key risk in structured products and exotic options
//
// ═══════════════════════════════════════════════════════════════════════════════

TEST(GreeksTest, VolgaBasic) {
    // GOAL: Verify volga is positive
    //
    // Volga must ALWAYS be positive for long options
    // Higher vol -> higher vega -> positive second derivative
    
    double volga = BlackScholesFormulas::volga(S0, K, r, sigma, T);
    
    EXPECT_GT(volga, 0.0)
        << "Volga must be positive (options have positive convexity in volatility)";
}

TEST(GreeksTest, VolgaMoneyness) {
    // GOAL: Verify volga is positive across all moneyness levels
    //
    // VOLGA PROFILE: More complex than vega
    // - ATM: High volga (maximum vega sensitivity)
    // - ITM/OTM: Can have secondary peaks depending on parameters
    // - Always positive but magnitude varies
    //
    // TRADING: "Long ATM straddles" = positive volga exposure
    
    double volgaITM = BlackScholesFormulas::volga(S0, 80.0, r, sigma, T);
    double volgaATM = BlackScholesFormulas::volga(S0, 100.0, r, sigma, T);
    double volgaOTM = BlackScholesFormulas::volga(S0, 120.0, r, sigma, T);
    
    EXPECT_GT(volgaITM, 0.0) << "ITM volga should be positive";
    EXPECT_GT(volgaATM, 0.0) << "ATM volga should be positive (typically maximum)";
    EXPECT_GT(volgaOTM, 0.0) << "OTM volga should be positive";
}

TEST(GreeksTest, VolgaAtMaturity) {
    // GOAL: Test boundary condition - volga vanishes at expiry
    //
    // AT EXPIRY: Vega = 0 -> no sensitivity to changes in vega -> volga = 0
    // The option value is deterministic, volatility irrelevant
    
    double T_zero = 0.0;
    double volga = BlackScholesFormulas::volga(S0, K, r, sigma, T_zero);
    
    EXPECT_NEAR(volga, 0.0, GREEK_TOL)
        << "At expiry, volga = 0 (vega is zero, so its derivative is zero too)";
}

TEST(GreeksTest, VolgaNumericalVerification) {
    // GOAL: Verify analytical volga formula using second-order finite differences
    //
    // METHOD: Compute ∂2V/∂σ2 using: [V(σ+h) - 2V(σ) + V(σ-h)] / h2
    // This validates the closed-form expression for vol convexity
    
    auto callPrice = [](double vol) {
        return BlackScholesFormulas::callPrice(S0, K, r, vol, T);
    };
    
    double analyticVolga = BlackScholesFormulas::volga(S0, K, r, sigma, T);
    double numericalVolga = numericalSecondDerivative(callPrice, sigma);
    
    EXPECT_NEAR(analyticVolga, numericalVolga, NUMERICAL_TOL * 10)
        << "Analytical volga should match numerical second derivative w.r.t. volatility";
}

// ============================================================================
// DISCOUNT CURVE VERSIONS
// ============================================================================

TEST(GreeksTest, DeltaWithDiscountCurve) {
    FlatDiscountCurve curve(r);
    
    // Compare flat rate version with discount curve version
    double deltaFlat = BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Call);
    double deltaCurve = BlackScholesFormulas::delta(S0, K, curve, sigma, T, Option::Type::Call);
    
    EXPECT_NEAR(deltaFlat, deltaCurve, GREEK_TOL);
}

TEST(GreeksTest, VegaWithDiscountCurve) {
    FlatDiscountCurve curve(r);
    
    double vegaFlat = BlackScholesFormulas::vega(S0, K, r, sigma, T);
    double vegaCurve = BlackScholesFormulas::vega(S0, K, curve, sigma, T);
    
    EXPECT_NEAR(vegaFlat, vegaCurve, GREEK_TOL);
}

TEST(GreeksTest, ThetaWithDiscountCurve) {
    FlatDiscountCurve curve(r);
    
    double thetaFlat = BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Call);
    double thetaCurve = BlackScholesFormulas::theta(S0, K, curve, sigma, T, Option::Type::Call);
    
    EXPECT_NEAR(thetaFlat, thetaCurve, GREEK_TOL * 10);
}

TEST(GreeksTest, GammaWithDiscountCurve) {
    FlatDiscountCurve curve(r);
    
    double gammaFlat = BlackScholesFormulas::gamma(S0, K, r, sigma, T);
    double gammaCurve = BlackScholesFormulas::gamma(S0, K, curve, sigma, T);
    
    EXPECT_NEAR(gammaFlat, gammaCurve, GREEK_TOL);
}

// ============================================================================
// CROSS-VALIDATION AND RELATIONSHIPS
// ============================================================================

TEST(GreeksTest, BlackScholesGammaVegaRelationship) {
    // For Black-Scholes, there's an exact relationship between Gamma and Vega:
    // This comes from the relationship between the standard normal PDF evaluations
    double gamma = BlackScholesFormulas::gamma(S0, K, r, sigma, T);
    double vega = BlackScholesFormulas::vega(S0, K, r, sigma, T);
    
    // Vega = S * φ(d1) * sqrt(T)
    // Gamma = φ(d1) / (S * σ * sqrt(T))
    // Rearranging: φ(d1) = Gamma * S * σ * sqrt(T)
    // Substituting: Vega = S * [Gamma * S * σ * sqrt(T)] * sqrt(T) = S^2 * Gamma * σ * T
    double expectedVega = S0 * S0 * gamma * sigma * T;
    
    // Tight tolerance - this is an exact relationship
    EXPECT_NEAR(vega, expectedVega, 0.01);
}

TEST(GreeksTest, ParameterSensitivityConsistency) {
    // All Greeks should change smoothly with parameters
    std::vector<double> spots = {90.0, 95.0, 100.0, 105.0, 110.0};
    std::vector<double> deltas;
    
    for (double spot : spots) {
        deltas.push_back(BlackScholesFormulas::delta(spot, K, r, sigma, T, Option::Type::Call));
    }
    
    // Delta should be monotonically increasing with spot
    for (size_t i = 1; i < deltas.size(); ++i) {
        EXPECT_GT(deltas[i], deltas[i-1]);
    }
}

// ============================================================================
// EXTREME PARAMETER TESTS
// ============================================================================

TEST(GreeksTest, HighVolatilityRegime) {
    // Greeks should still be well-defined for high volatility
    double high_vol = 1.0;  // 100% volatility
    
    double delta = BlackScholesFormulas::delta(S0, K, r, high_vol, T, Option::Type::Call);
    double gamma = BlackScholesFormulas::gamma(S0, K, r, high_vol, T);
    double vega = BlackScholesFormulas::vega(S0, K, r, high_vol, T);
    
    EXPECT_TRUE(std::isfinite(delta));
    EXPECT_TRUE(std::isfinite(gamma));
    EXPECT_TRUE(std::isfinite(vega));
    EXPECT_GT(vega, 0.0);
}

TEST(GreeksTest, LongMaturity) {
    // Greeks should be well-defined for long maturities
    double long_T = 5.0;  // 5 years
    
    double delta = BlackScholesFormulas::delta(S0, K, r, sigma, long_T, Option::Type::Call);
    double vega = BlackScholesFormulas::vega(S0, K, r, sigma, long_T);
    
    EXPECT_TRUE(std::isfinite(delta));
    EXPECT_TRUE(std::isfinite(vega));
    
    // Long-dated ATM call delta will be higher due to long-term drift r*T
    // For r=5%, T=5, the forward is much higher than spot
    EXPECT_GT(delta, 0.5);
    EXPECT_LT(delta, 1.0);
}

TEST(GreeksTest, ShortMaturity) {
    // Greeks should handle short maturities
    double short_T = 0.01;  // ~3.65 days
    
    double gamma = BlackScholesFormulas::gamma(S0, K, r, sigma, short_T);
    double theta = BlackScholesFormulas::theta(S0, K, r, sigma, short_T, Option::Type::Call);
    
    EXPECT_TRUE(std::isfinite(gamma));
    EXPECT_TRUE(std::isfinite(theta));
    
    // Short-dated ATM options have large gamma (pin risk)
    EXPECT_GT(gamma, BlackScholesFormulas::gamma(S0, K, r, sigma, T));
}
