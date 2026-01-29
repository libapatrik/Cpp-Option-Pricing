//
// Created by Patrik  Liba on 24/10/2025.
//

#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cppfm/utils/Utils.h>

#include <cmath>
#include <algorithm>

// ============================================================================
// Helper Functions
// ============================================================================
namespace {
    /**
     * Standard normal cumulative distribution function
     */
    double normCdf(double x) {
        return Utils::stdNormCdf(x);
    }
    
    /**
     * Standard normal probability density function
     * φ(x) = (1/sqrt(2*π)) * e^(-x²/2)
     */
    double normPdf(double x) {
        return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * x * x);
    }
}

// ============================================================================
// d1/d2 Helpers - Flat Rate Version
// ============================================================================

/**
 * Black-Scholes d1 parameter with flat rate
 * d₁ = [ln(S/K) + (r + σ²/2)T] / (σ√T)
 */
double BlackScholesFormulas::d1(double spot, double strike, double rate, double volatility, double maturity)
{
    // Check edge cases
    if (maturity <= 0.0 || volatility <= 0.0)
    {
        return 0.0;
    }
    const double sqrtTime = std::sqrt(maturity);
    const double numerator = std::log(spot / strike) + (rate + 0.5 * volatility * volatility) * maturity;
    const double denominator = volatility * sqrtTime;  

    return numerator / denominator;
}

/**
 * Black-Scholes d2 parameter with flat rate
 * d2 = d1 - σ√T
 * 
 * Alternative formula: d2 = [ln(S/K) + (r - σ^2/2)T] / (σ√T)
 */
double BlackScholesFormulas::d2(double spot, double strike, double rate, double volatility, double maturity)
{
    // Check edge case
    if (maturity <= 0.0 || volatility <= 0.0)
    {
        return 0.0;
    }
    // Faster: just subtract σ√T from d1
    return d1(spot, strike, rate, volatility, maturity) - volatility * std::sqrt(maturity);
}

// ============================================================================
// d1/d2 Helpers - Discount Curve Version (Time-Dependent Rates)
// ============================================================================

/**
 * Black-Scholes d1 parameter with DiscountCurve
 * 
 * KEY DIFFERENCE from flat rate:
 * - Flat rate: uses r·T in the numerator
 * - Discount curve: uses -ln(B(T)) which equals int_0^T r(s)ds
 * 
 * Formula: d₁ = [ln(S/K) - ln(B(T)) + σ^2T/2] / (σ√T)
 * 
 * Where B(T) = discountCurve.discount(T) is the discount factor at time T
 * 
 * Financial Interpretation:
 * -ln(B(T)) = integral of instantaneous forward rates from 0 to T
 * This handles time-varying interest rates properly!
 */
double BlackScholesFormulas::d1(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity)
{
    // Check edge cases
    if (maturity <= 0.0 || volatility <= 0.0)
    {
        return 0.0;
    }
    
    // Get discount factor B(T) from the curve
    const double discountFactor = discountCurve.discount(maturity);
    const double sqrtTime = std::sqrt(maturity);

    // Use -ln(B(T)) instead of r·T
    // So equivalent to int_0^T r(s)ds for time-dependent rates
    const double numerator = std::log(spot / strike) - std::log(discountFactor) + 0.5 * volatility * volatility * maturity;
    const double denominator = volatility * sqrtTime;

    return numerator / denominator;
}

/**
 * Black-Scholes d2 parameter with DiscountCurve
 * d2 = d1 - σ√T (same relationship as flat rate)
 */
double BlackScholesFormulas::d2(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity)
{
    if (maturity <= 0.0 || volatility <= 0.0)
    {
        return 0.0;
    }
    // Faster: reuse d1 calculation
    return d1(spot, strike, discountCurve, volatility, maturity) - volatility * std::sqrt(maturity);
}

// ============================================================================
// Pricing Functions - Flat Rate Version
// ============================================================================

/**
 * Generic option pricing with option type as parameter
 * This is the CLEAN way - no if/else in caller code!
 */
double BlackScholesFormulas::price(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType)
{
    return (optionType == Option::Type::Call)
            ? callPrice(spot, strike, rate, volatility, maturity)
            : putPrice(spot, strike, rate, volatility, maturity);
}

/**
 * Black-Scholes Call Option Formula
 * C = S*N(d1) - K*e^(-rT)*N(d2)
 */
double BlackScholesFormulas::callPrice(double spot, double strike, double rate, double volatility, double maturity)
{
    // At maturity: return intrinsic value
    if (maturity <= 0.0)
    {
        return std::max(spot - strike, 0.0);  // max(S - K, 0)
    }
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double discountFactor = std::exp(-rate * maturity);

    return spot * normCdf(d1Val) - strike * discountFactor * normCdf(d2Val);
}

/**
 * Black-Scholes Put Option Formula
 * P = K*e^(-rT)*N(-d2) - S*N(-d1)
 */
double BlackScholesFormulas::putPrice(double spot, double strike, double rate, double volatility, double maturity)
{
    // At maturity: return intrinsic value
    if (maturity <= 0.0)
    {
        return std::max(strike - spot, 0.0);  // max(K - S, 0)
    }
    
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double discountFactor = std::exp(-rate * maturity);

    return strike * discountFactor * normCdf(-d2Val) - spot * normCdf(-d1Val);
}

// ============================================================================
// Pricing Functions - Discount Curve Version 
// ============================================================================

/**
 * Generic option pricing with DiscountCurve
 */
double BlackScholesFormulas::price(double spot, double strike, const DiscountCurve& discountCurve, 
                                   double volatility, double maturity, Option::Type optionType)
{
    return (optionType == Option::Type::Call)
            ? callPrice(spot, strike, discountCurve, volatility, maturity)
            : putPrice(spot, strike, discountCurve, volatility, maturity);
}

/**
 * Black-Scholes Call with DiscountCurve
 * 
 * Use discountCurve.discount(T) instead of exp(-r*T)
 */
double BlackScholesFormulas::callPrice(double spot, double strike, const DiscountCurve& discountCurve, 
                                       double volatility, double maturity)
{
    if (maturity <= 0.0)
    {
        return std::max(spot - strike, 0.0);  
    }
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    const double d2Val = d2(spot, strike, discountCurve, volatility, maturity);

    const double discountFactor = discountCurve.discount(maturity);  

    return spot * normCdf(d1Val) - strike * discountFactor * normCdf(d2Val);
}

/**
 * Black-Scholes Put with DiscountCurve
 */
double BlackScholesFormulas::putPrice(double spot, double strike, const DiscountCurve& discountCurve, 
                                      double volatility, double maturity)
{
    if (maturity <= 0.0)
    {
        return std::max(strike - spot, 0.0);  
    }
    
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    const double d2Val = d2(spot, strike, discountCurve, volatility, maturity);
    
    const double discountFactor = discountCurve.discount(maturity);

    return strike * discountFactor * normCdf(-d2Val) - spot * normCdf(-d1Val);
}

/**
 * Greeks: first order derivatives
 */
/**
 * Delta: ∂V/∂S (rate of change with respect to spot)
 * Call Delta: Δ = N(d1) in [0, 1]
 * Put Delta: Δ = N(d1) - 1 in [-1, 0]
 */
double BlackScholesFormulas::delta(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType)
{
    if (maturity <= 0.0)
    { 
        // At maturity: Delta = 1 if ITM, 0 if OTM
        if (optionType == Option::Type::Call)
        {
            return (spot > strike) ? 1.0 : 0.0;
        } else
        {
            return (spot < strike) ? -1.0 : 0.0;
        }
    }
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    return (optionType == Option::Type::Call)
            ? normCdf(d1Val)
            : - normCdf(-d1Val); // or normCdf(d1Val) -1
}

double BlackScholesFormulas::delta(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity, Option::Type optionType)
{
    if (maturity <= 0.0) {
        if (optionType == Option::Type::Call) {
            return (spot > strike) ? 1.0 : 0.0;
        } else {
            return (spot < strike) ? -1.0 : 0.0;
        }
    }
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    return (optionType == Option::Type::Call)
            ? normCdf(d1Val)
            : - normCdf(-d1Val);
}

double BlackScholesFormulas::vega(double spot, double strike, double rate, double volatility, double maturity)
{ // with given rate
    if (maturity <= 0.0 || volatility <= 0.0) {
        return 0.0;
    }
    // Vega: ν = S·φ(d1)·sqrt(T)
    //
    // Same for calls and puts (by put-call parity)
    //
    // Financial interpretation:
    // - Measures P&L impact of 1% volatility change
    // - Maximum for ATM options
    // - Decreases as option moves ITM or OTM
    // - Increases with time to maturity
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);    // vega/100?

    return spot * normPdf(d1Val) * sqrtTime;
}

double BlackScholesFormulas::vega(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity)
{ // with instantaneous curve
    if (maturity <= 0.0 || volatility <= 0.0) {
        return 0.0;
    }
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);

    return spot * normPdf(d1Val) * sqrtTime;  // vega/100?
}

double BlackScholesFormulas::theta(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType)
{ // with given rate
    if (maturity <= 0.0) {
        return 0.0;
    }
    // Theta measures time decay: Θ = ∂V/∂t (negative for long positions)
    //
    // Call Theta: Θ_call = -S·φ(d₁)·σ/(2√T) - rK·e^(-rT)·N(d₂)
    // Put Theta:  Θ_put = -S·φ(d₁)·σ/(2√T) + rK·e^(-rT)·N(-d₂)
    //
    // Financial interpretation:
    // - Options lose value as time passes (time decay)
    // - Theta is most negative for ATM options near expiry
    // - Long options have negative Theta (lose money daily)
    // - Short options have positive Theta (make money daily)
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);
    const double discountFactor = std::exp(-rate * maturity);

    // Common term for both calls and puts
    const double term1 = -spot * normPdf(d1Val) * volatility / (2.0 * sqrtTime);

    if (optionType == Option::Type::Call) {
        const double term2 = -rate * strike * discountFactor * normCdf(d2Val);
        return term1 + term2;
    } else {
        const double term2 = rate * strike * discountFactor * normCdf(-d2Val);
        return term1 + term2;
    }
}

double BlackScholesFormulas::theta(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity, Option::Type optionType)
{
    if (maturity <= 0.0) {
        return 0.0;
    }
    // For time-dependent rates, we need the instantaneous rate r(t)
    const double rate = discountCurve.instantaneousRate(maturity);
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    const double d2Val = d2(spot, strike, discountCurve, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);
    const double discountFactor = discountCurve.discount(maturity);

    const double term1 = -spot * normPdf(d1Val) * volatility / (2.0 * sqrtTime);

    if (optionType == Option::Type::Call) {
        const double term2 = -rate * strike * discountFactor * normCdf(d2Val);
        return term1 + term2;
    } else {
        const double term2 = rate * strike * discountFactor * normCdf(-d2Val);
        return term1 + term2;    // θ/365?
    }
}

double BlackScholesFormulas::rho(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType)
{
    if (maturity <= 0.0) {
        return 0.0;
    }

    // Rho: ρ = ∂V/∂r.      
    //
    // Call Rho: ρ_call = K·T·e^(-rT)·N(d2)  (positive)
    // Put Rho:  ρ_put = -K·T·e^(-rT)·N(-d₂) (negative)
    //
    // Financial interpretation:
    // - Measures sensitivity to interest rate changes
    // - Call options benefit from higher rates (positive Rho)
    // - Put options suffer from higher rates (negative Rho)

    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double discountFactor = std::exp(-rate * maturity);

    if (optionType == Option::Type::Call) {
        return strike * maturity * discountFactor * normCdf(d2Val);        // ρ/100?
    } else {
        return -strike * maturity * discountFactor * normCdf(-d2Val);    // ρ/100?
    }
    // TODO: Rho with dividend yield
}

// ============================================================================
// Greeks - Second Order (Convexity)
// ============================================================================

/**
 * Gamma: ∂2V/∂S^2 (convexity - how fast Delta changes)
 * Same for calls and puts
 * Formula: Γ = φ(d1) / (S*σ*√T)
 */
double BlackScholesFormulas::gamma(double spot, double strike, double rate, double volatility, double maturity)
{
    if (maturity <= 0.0 || volatility <= 0.0 || spot <= 0.0) {
        return 0.0;
    }
    // Formula: Γ = φ(d1) / (S*σ*√T)
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);
    return normPdf(d1Val) / (spot * volatility * sqrtTime);   // γ/100?
}

/**
 * Gamma with DiscountCurve
 */
double BlackScholesFormulas::gamma(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity)
{
    if (maturity <= 0.0 || volatility <= 0.0 || spot <= 0.0) {
        return 0.0;
    }
    const double d1Val = d1(spot, strike, discountCurve, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);

    return normPdf(d1Val) / (spot * volatility * sqrtTime);   // γ/100?
}

/** NOTE:
 * Greeks like vega, vanna and Volga/vomma that involve partial differentials with respect to  are in some sense ‘invalid’
 * in the context of Black-Scholes, since in its derivation we assume that  is constant.
 * We might interpret them as applying to a model in which  was slightly variable but otherwise was close to constant for all,
 *  etc. Vega, for example, would then measure the sensitivity to changes in the mean level of.
 *  For some types of derivatives, e.g. binary puts and calls,
 *  it can be difficult to interpret how these particular sensitivities should be understood.
 *
 * Vanna: ∂2V/∂S∂σ (cross-Greek)
 * Same for calls and puts!
 */
double BlackScholesFormulas::vanna(double spot, double strike, double rate, double volatility, double maturity)
{
    if (maturity <= 0.0 || volatility <= 0) {
        return 0.0;
    }
    // Formula: Vanna = -φ(d1) * d2/(σ·√T)
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);

    return -normPdf(d1Val) * d2Val / (volatility * sqrtTime);
}

/**
 * Volga: ∂2V/∂σ^2 (vol convexity)
 * Same for calls and puts
 */
double BlackScholesFormulas::volga(double spot, double strike, double rate, double volatility, double maturity)
{
    if (maturity <= 0 || volatility <= 0) {
        return 0.0;
    }

    // Formula: Volga = S*φ(d1) * √T * d1*d2/σ
    // ∂^2V/∂σ^2 = vega(d1*d2/σ) 
    // Measures convexity of option value with respect to volatility
    //
    // Financial interpretation:
    // - How Vega changes with volatility
    // - Positive for vanilla options (Vega increases with vol)
    // - Important for volatility smile dynamics
    // - Used in vol-of-vol hedging strategies
    const double d1Val = d1(spot, strike, rate, volatility, maturity);
    const double d2Val = d2(spot, strike, rate, volatility, maturity);
    const double sqrtTime = std::sqrt(maturity);

    // return (spot * normPdf(d1Val) * sqrtTime * d1Val * d2Val / volatility);  // volga/100? Or, no! Δ has no division.
    return vega(spot, strike, rate, volatility, maturity) * (d1Val * d2Val / volatility); // gives the same result
}

// ============================================================================
// BS Implied Volatility - Newton-Raphson
// ============================================================================

/** GRZ Chapter 4
 * Implied Volatility with flat rate
 * Uses Newton-Raphson: σ_{n+1} = σ_n - (BS(σ_n) - Price) / Vega(σ_n)
 */
double BlackScholesFormulas::impliedVolatility(double spot, double strike, double rate, double maturity,
                                               double marketPrice, Option::Type optionType,
                                               double initialGuess, size_t maxIter, double tol)
{
    // Edge case: at or past maturity
    if (maturity <= 0.0) {
        return 0.0;
    }

    double sigma = initialGuess;
    const double minVol = 1e-6;
    const double maxVol = 5.0;

    for (size_t i = 0; i < maxIter; ++i) {
        double bsPrice = price(spot, strike, rate, sigma, maturity, optionType);
        double vegaVal = vega(spot, strike, rate, sigma, maturity);

        double diff = bsPrice - marketPrice;

        // Check convergence
        if (std::abs(diff) < tol) {
            return sigma;
        }

        // Safeguard against zero vega
        if (std::abs(vegaVal) < 1e-12) {
            // Fallback to bisection step
            if (diff > 0) {
                sigma *= 0.5;
            } else {
                sigma *= 1.5;
            }
        } else {
            // Newton-Raphson update
            sigma -= diff / vegaVal;
        }

        // Keep within bounds
        sigma = std::max(minVol, std::min(sigma, maxVol));
    }

    return sigma; // Return best guess if not converged
}

/**
 * Implied Volatility with DiscountCurve
 */
double BlackScholesFormulas::impliedVolatility(double spot, double strike, const DiscountCurve& discountCurve,
                                               double maturity, double marketPrice, Option::Type optionType,
                                               double initialGuess, size_t maxIter, double tol)
{
    if (maturity <= 0.0) {
        return 0.0;
    }

    double sigma = initialGuess;
    const double minVol = 1e-6;
    const double maxVol = 5.0;

    for (size_t i = 0; i < maxIter; ++i) {
        double bsPrice = price(spot, strike, discountCurve, sigma, maturity, optionType);
        double vegaVal = vega(spot, strike, discountCurve, sigma, maturity);

        double diff = bsPrice - marketPrice;

        if (std::abs(diff) < tol) {
            return sigma;
        }

        if (std::abs(vegaVal) < 1e-12) {
            if (diff > 0) {
                sigma *= 0.5;
            } else {
                sigma *= 1.5;
            }
        } else {
            sigma -= diff / vegaVal;
        }

        sigma = std::max(minVol, std::min(sigma, maxVol));
    }

    return sigma;
}