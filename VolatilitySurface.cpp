//
// Created by Patrik  Liba on 05/10/2025.
//

#include "VolatilitySurface.h"
#include "BlackScholesFormulas.h"
#include "Utils.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>

/**
 * VolatilitySurface Implementation
 */
VolatilitySurface::VolatilitySurface(const std::vector<double>& strikes,                        // [strike]; 1d vector
                                     const std::vector<double>& maturities,                     // [maturity]; 1d vector
                                     const std::vector<std::vector<double>>& volatilities,      // [maturity][strike]; 2d matrix
                                     const DiscountCurve& discountCurve,                        // discount curve
                                     InterpolationType interpolationType)                       // interpolation type
    : _strikes(strikes), _maturities(maturities), _volatilities(volatilities),
      _interpolationType(interpolationType), _discountCurve(discountCurve.clone())
{
    // Validate input data
    if (_strikes.empty() || _maturities.empty() || _volatilities.empty()) {
        throw std::invalid_argument("VolatilitySurface: Empty data provided");
    }
    
    if (_volatilities.size() != _maturities.size()) {
        throw std::invalid_argument("VolatilitySurface: Mismatch between maturities and volatility matrix rows");
    }
    
    for (const auto& row : _volatilities) {
        if (row.size() != _strikes.size()) {
            throw std::invalid_argument("VolatilitySurface: Mismatch between strikes and volatility matrix columns");
        }
    }
    
    // Check for sorted data
    if (!std::is_sorted(_strikes.begin(), _strikes.end())) {
        throw std::invalid_argument("VolatilitySurface: Strikes must be sorted");
    }
    
    if (!std::is_sorted(_maturities.begin(), _maturities.end())) {
        throw std::invalid_argument("VolatilitySurface: Maturities must be sorted");
    }
    
    // Check for positive volatilities
    for (const auto& row : _volatilities) {
        for (double vol : row) {
            if (vol <= 0.0) {
                throw std::invalid_argument("VolatilitySurface: All volatilities must be positive");
            }
        }
    }
    
    initializeInterpolators(); // proceeds to next step - call private method to create interpolators for strikes and maturities
}
/**
 *   Smile Interpolation: For each maturity, interpolate volatility across strikes (the volatility smile)
 *   Term Structure Interpolation: For each strike, interpolate volatility across maturities (how volatility changes over time)
 */
void VolatilitySurface::initializeInterpolators()
{
    // Initialize 1D interpolators for each maturity (smile interpolation)
    // cubic spline interpolation along strikes
    _smileInterpolators.reserve(_maturities.size());                   // reserve space for maturities in _smileInterpolators vector
    for (size_t i = 0; i < _maturities.size(); ++i) {                    // for each maturity
        if (_interpolationType == InterpolationType::CubicSpline) {      // check type: do cubic spline interpolation
            _smileInterpolators.push_back(
                std::make_unique<CubicSplineInterpolation>(            // creates a unique pointer to a new CubicSplineInterpolation object
                    _strikes, _volatilities[i], 
                    CubicSplineInterpolation::BoundaryType::Natural));   // boundary type: natural
        } else {
            // Linear interpolation for smiles
            _smileInterpolators.push_back(
                std::make_unique<LinearInterpolation>(_strikes, _volatilities[i]));
        }
    }
    
    // Initialize 1D interpolators for each strike (term structure interpolation)
    _termStructureInterpolators.reserve(_strikes.size());          // reserve space for strikes
    for (size_t j = 0; j < _strikes.size(); ++j) {                   // for each strike
        std::vector<double> termVols;                                // volatilities across maturities for strike j
        termVols.reserve(_maturities.size());                      // reserve space for maturities
        for (size_t i = 0; i < _maturities.size(); ++i) {            // for each maturity
            termVols.push_back(_volatilities[i][j]);                 // collect volatilities for strike j
        }
        
        if (_interpolationType == InterpolationType::CubicSpline) {
            _termStructureInterpolators.push_back(
                std::make_unique<CubicSplineInterpolation>(
                    _maturities, termVols,
                    CubicSplineInterpolation::BoundaryType::Natural));
        } else {
            _termStructureInterpolators.push_back(
                std::make_unique<LinearInterpolation>(_maturities, termVols));
        }
    }
}

double VolatilitySurface::impliedVolatility(double strike, double maturity) const
{
    // Check bounds
    if (strike < _strikes.front() || strike > _strikes.back() ||
        maturity < _maturities.front() || maturity > _maturities.back()) {
        throw std::out_of_range("VolatilitySurface: Point outside interpolation range");
    }
    
    
    // For other interpolation types, use 1D interpolation along maturity first
    // Find the closest maturity index
    auto maturityIt = std::lower_bound(_maturities.begin(), _maturities.end(), maturity); // first element >=
    size_t maturityIdx = std::distance(_maturities.begin(), maturityIt);
    
    if (maturityIdx == 0) {
        return _smileInterpolators[0]->interpolate(strike);
    }
    
    if (maturityIdx >= _maturities.size()) {
        return _smileInterpolators.back()->interpolate(strike);
    }
    
    // Interpolate between adjacent maturities
    double t1 = _maturities[maturityIdx - 1]; // Check careful here!
    double t2 = _maturities[maturityIdx];
    double vol1 = _smileInterpolators[maturityIdx - 1]->interpolate(strike);
    double vol2 = _smileInterpolators[maturityIdx]->interpolate(strike);
    
    // Linear interpolation in time
    double weight = (maturity - t1) / (t2 - t1);
    return vol1 + weight * (vol2 - vol1);
}

double VolatilitySurface::impliedVolatilityForwardMoneyness(double strike, double maturity, double spot) const
{
    // Section 1.3 - Interpolation along maturities
    // Implement equation (1.25) with forward moneyness k_F_T = K/F_T
    
    // LECTURE NOTES: Use DiscountCurve
    // Forward price: F_T = S_0 * e^{\int_0^T r(s)ds} = S_0 / B(T)
    // where B(T) = e^{-\int_0^T r(s)ds} is the discount factor
    double forward = spot / _discountCurve->discount(maturity);
    double forwardMoneyness = strike / forward;
    
    // Step 2: Find maturity interval
    auto maturityIt = std::lower_bound(_maturities.begin(), _maturities.end(), maturity); // careful then decide idx
    size_t maturityIdx = std::distance(_maturities.begin(), maturityIt);
    
    if (maturityIdx == 0) {
        // Extrapolation: T < T_1
        double t1 = _maturities[0];
        double forward1 = spot / _discountCurve->discount(t1);
        double strike1 = forwardMoneyness * forward1;
        double vol1 = _smileInterpolators[0]->interpolate(strike1);
        return vol1;
    }
    
    if (maturityIdx >= _maturities.size()) {
        // Extrapolation: T > T_M
        double tM = _maturities.back();
        double forwardM = spot / _discountCurve->discount(tM);
        double strikeM = forwardMoneyness * forwardM;
        double volM = _smileInterpolators.back()->interpolate(strikeM);
        return volM;
    }
    
    // Step 3: Interpolation between T_i and T_{i+1}
    double t1 = _maturities[maturityIdx - 1];
    double t2 = _maturities[maturityIdx];
    
    // Step 4: Extract strikes corresponding to forward moneyness
    double forward1 = spot / _discountCurve->discount(t1);
    double forward2 = spot / _discountCurve->discount(t2);
    double strike1 = forwardMoneyness * forward1;
    double strike2 = forwardMoneyness * forward2;
    
    // Step 5: Get volatilities at these strikes
    double vol1 = _smileInterpolators[maturityIdx - 1]->interpolate(strike1);
    double vol2 = _smileInterpolators[maturityIdx]->interpolate(strike2);
    
    // Step 6: Linear interpolation in variance
    // Equation v(T,k) = (σ*)^2(T, k * S0 e^{\int_0^T r(s)ds}) * T
    double variance1 = vol1 * vol1 * t1;
    double variance2 = vol2 * vol2 * t2;
    
    // Linear interpolation in time for variance
    double weight = (maturity - t1) / (t2 - t1);
    double interpolatedVariance = variance1 + weight * (variance2 - variance1);
    
    // Step 7: Convert back to volatility
    // Equation (1.25) σ*(T,K) = sqrt(v(T,k_F_T) / T)
    // Ensure we don't divide by zero
    if (maturity <= 1e-12) {
        return vol1; // Return the closer volatility for very small times
    }
    return std::sqrt(std::max(interpolatedVariance / maturity, 1e-12));
}

double VolatilitySurface::localVolatility(double spot, double time) const
{
    return computeDupireLocalVolatility(spot, time);
}

/**
 * Dupire Local Volatility Formula - MAIN CALCULATION HERE
 */
double VolatilitySurface::computeDupireLocalVolatility(double spot, double time) const
{
    /*
     * Dupire Local Volatility Formula
     * - What local volatility function σ(S,t) would reproduce the observed market prices?
     *
     * Mathematical Formula:
     * sigma_D(T,K)^2 = sigma*(T,K)^2 * [numerator] / [denominator]
     * 
     * Where:
     * - sigma*(T,K) is the implied volatility at strike K and time T
     * - numerator = 1 + 2T/sigma* * (dSigma/dT + r*dSigma/dK)
     * - denominator = 1 + 2d1*(K*dSigma/dK*sqrtT) + d1*d2*(K*dSigma/dK*sqrtT)^2 + (K*d2Sigma/dK2*sqrtT)*(K*sigma*sqrtT)
     * 
     * Formula:
     * 1. Accounts for volatility smile curvature (d2Sigma/dK2)
     * 2. Considers time decay of the smile (dSigma/dT)  
     * 3. Incorporates strike sensitivity (dSigma/dK)
     * 4. Uses Black-Scholes d1 and d2 parameters for proper weighting
     */
    
    double strike = spot; // S = K; for local volatility at current spot level
    
    // LECTURE NOTES: Use DiscountCurve properly for time-dependent rates
    // For time-dependent rates, we need the instantaneous rate r(t) at time t
    // This is computed as: r(t) = -d/dt[log(B(t))] where B(t) is the discount factor
    // For numerical approximation: r(t) = -[log(B(t+ε)) - log(B(t))]/ε

    // const double eps = 1e-6;
    // double discount_t = _discountCurve->discount(time);
    // double discount_t_plus = _discountCurve->discount(time + eps);

    // No mkt data for short-rate
    // double riskFreeRate = -std::log(discount_t_plus / discount_t) / eps;
    double riskFreeRate = _discountCurve->instantaneousRate(time);

    // Get implied volatility and its derivatives from the interpolated surface
    double impliedVol = impliedVolatility(strike, time);
    double dSigma_dT = impliedVolatilityDerivativeTime(strike, time);
    double dSigma_dK = impliedVolatilityDerivativeStrike(strike, time);
    double d2Sigma_dK2 = impliedVolatilitySecondDerivativeStrike(strike, time);
    
    // Compute Black-Scholes d1 and d2 parameters
    // for denominator
    double sqrtTime = impliedVol * std::sqrt(time);
    double d1 = (std::log(spot / strike) + (riskFreeRate + 0.5 * impliedVol * impliedVol) * time) / sqrtTime;
    double d2 = d1 - sqrtTime;
    
    // Compute the numerator: 1 + 2T/σ* × (∂σ*/∂T + r×∂σ*/∂K)
    // This captures how the implied volatility surface changes over time and across strikes
    double numerator = 1.0 + (2.0 * time / impliedVol) * (dSigma_dT + riskFreeRate * dSigma_dK);
    
    // Compute the denominator components
    double K_dSigma_dK_sqrtT = strike * dSigma_dK * std::sqrt(time);
    double K_d2Sigma_dK2_sqrtT = strike * d2Sigma_dK2 * std::sqrt(time);
    double K_sigma_sqrtT = strike * impliedVol * std::sqrt(time);
    
    // The denominator captures: 
    // - Linear smile effects (2d1 term)
    // - Quadratic smile effects (d1d2 term)
    // - Curvature effects (∂2σ*/∂K2 term)
    double denominator = 1.0 + 2.0 * d1 * K_dSigma_dK_sqrtT + 
                        d1 * d2 * K_dSigma_dK_sqrtT * K_dSigma_dK_sqrtT +
                        K_d2Sigma_dK2_sqrtT * K_sigma_sqrtT;
    
    // Safety checks: ensure we don't get invalid results
    if (denominator <= 0.0 || !std::isfinite(numerator) || !std::isfinite(denominator)) {
        return impliedVol; // Fallback to implied volatility if computation fails
    }
    
    // Apply the Dupire formula
    double localVolSquared = impliedVol * impliedVol * numerator / denominator;
    
    // Ensure we don't take square root of negative number
    // If localVol < 0, floor it at very small value. 0.03%
    if (localVolSquared < 0.0) {
        return impliedVol; // Fallback to implied volatility
    }

    double localVol = std::sqrt(localVolSquared);
    
    return localVol;
}


double VolatilitySurface::blackScholesCall(double spot, double strike, double time, double volatility) const
{
    return BlackScholesFormulas::callPrice(spot, strike, *_discountCurve, volatility, time);
}

double VolatilitySurface::blackScholesPut(double spot, double strike, double time, double volatility) const
{
    return BlackScholesFormulas::putPrice(spot, strike, *_discountCurve, volatility, time);
}

double VolatilitySurface::blackScholesVega(double spot, double strike, double time, double volatility) const
{
    return BlackScholesFormulas::vega(spot, strike, *_discountCurve, volatility, time);
}

double VolatilitySurface::blackScholesGamma(double spot, double strike, double time, double volatility) const
{
    return BlackScholesFormulas::gamma(spot, strike, *_discountCurve, volatility, time);
}

double VolatilitySurface::blackScholesTheta(double spot, double strike, double time, double volatility) const
{
    // Theta for Call - Dupire formula typically uses call theta
    return BlackScholesFormulas::theta(spot, strike, *_discountCurve, volatility, time, Option::Type::Call);
}



/**
 * Implied Volatility Time Derivative (dSigma/dT)
 */
double VolatilitySurface::impliedVolatilityDerivativeTime(double strike, double maturity) const
{
    /* Implied Volatility Time Derivative (dSigma/dT)
     * 
     * This computes how the implied volatility changes as we move along the time dimension.
     * Dupire formula uses it to account for time decay of the volatility smile.
     * 
     * Uses finite difference approximation:
     * dSigma/dT = [sigma*(K, T+eps) - sigma*(K, T)] / eps
     * 
     * A positive value means implied volatility increases with time (smile steepens).
     * A negative value means implied volatility decreases with time (smile flattens).
     */
    
    const double eps = 1e-4; // Small time step for finite difference
    
    double vol_center = impliedVolatility(strike, maturity);
    double vol_time_plus = impliedVolatility(strike, maturity + eps);
    
    return (vol_time_plus - vol_center) / eps;
}

/**
 * Implied Volatility Strike Derivative (dSigma/dK)
 */
double VolatilitySurface::impliedVolatilityDerivativeStrike(double strike, double maturity) const
{
    /*
     * Implied Volatility Strike Derivative (dSigma/dK)
     * 
     * This computes how the implied volatility changes as we move across the strike dimension.
     * It measures the slope of the volatility smile at a given point.
     * 
     * Uses centered finite difference approximation:
     * dSigma/dK = [sigma*(K+eps, T) - sigma*(K-eps, T)] / (2*eps)
     * 
     * A positive value means implied volatility increases with strike (upward sloping smile).
     * A negative value means implied volatility decreases with strike (downward sloping smile).
     * Near zero indicates a flat smile at that point.
     */
    
    const double eps = 1e-4; // Small strike step for finite difference
    
    double vol_strike_plus = impliedVolatility(strike + eps, maturity);
    double vol_strike_minus = impliedVolatility(strike - eps, maturity);
    
    return (vol_strike_plus - vol_strike_minus) / (2.0 * eps);
}

/**
 * Implied Volatility Second Strike Derivative (d2Sigma/dK2)
 */
double VolatilitySurface::impliedVolatilitySecondDerivativeStrike(double strike, double maturity) const
{
    /*
     * Implied Volatility Second Strike Derivative (d2Sigma/dK2)
     * 
     * This computes the curvature of the volatility smile at a given point.
     * It measures how the slope of the smile changes across strikes.
     * 
     * Uses finite difference approximation:
     * d2Sigma/dK2 ≈ [sigma*(K+eps, T) - 2*sigma*(K, T) + sigma*(K-eps, T)] / eps^2
     * 
     * A positive value indicates convexity (smile curves upward like a parabola).
     * A negative value indicates concavity (smile curves downward).
     * Near zero indicates the smile is approximately linear at that point.
     * 
     * In Dupire's formula this accounts for smile curvature effects.
     */
    
    const double eps = 1e-4; // Small strike step for finite difference
    
    double vol_center = impliedVolatility(strike, maturity);
    double vol_strike_plus = impliedVolatility(strike + eps, maturity);
    double vol_strike_minus = impliedVolatility(strike - eps, maturity);
    
    return (vol_strike_plus - 2.0 * vol_center + vol_strike_minus) / (eps * eps);
}


double VolatilitySurface::impliedVolatilityTimeDerivative(double strike, double maturity) const
{
    return impliedVolatilityDerivativeTime(strike, maturity);
}

double VolatilitySurface::impliedVolatilityStrikeDerivative(double strike, double maturity) const
{
    return impliedVolatilityDerivativeStrike(strike, maturity);
}

double VolatilitySurface::impliedVolatilitySecondStrikeDerivative(double strike, double maturity) const
{
    return impliedVolatilitySecondDerivativeStrike(strike, maturity);
}


std::vector<double> VolatilitySurface::volatilitySmile(double maturity, const std::vector<double>& strikes) const
{
    std::vector<double> result;
    result.reserve(strikes.size());
    
    for (double strike : strikes) {
        result.push_back(impliedVolatility(strike, maturity));
    }
    
    return result;
}

std::vector<double> VolatilitySurface::volatilityTermStructure(double strike, const std::vector<double>& maturities) const
{
    std::vector<double> result;
    result.reserve(maturities.size());
    
    for (double maturity : maturities) {
        result.push_back(impliedVolatility(strike, maturity));
    }
    
    return result;
}


std::pair<std::pair<double, double>, std::pair<double, double>> VolatilitySurface::getBounds() const
{
    return {{_strikes.front(), _strikes.back()}, 
            {_maturities.front(), _maturities.back()}};
}

VolatilitySurface* VolatilitySurface::clone() const
{
    return new VolatilitySurface(_strikes, _maturities, _volatilities, *_discountCurve, _interpolationType);
}

bool VolatilitySurface::operator==(const VolatilitySurface& other) const
{
    return _strikes == other._strikes &&
           _maturities == other._maturities &&
           _volatilities == other._volatilities &&
           _interpolationType == other._interpolationType;
}


/**
 * VolatilitySurfaceBuilding Implementation
 */
VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::addStrike(double strike)
{
    _strikes.push_back(strike);
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::addMaturity(double maturity)
{
    _maturities.push_back(maturity);
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::setVolatility(double strike, double maturity, double volatility)
{
    // Find or create entry for this (strike-maturity pair)
    auto strikeIt = std::find(_strikes.begin(), _strikes.end(), strike);
    auto maturityIt = std::find(_maturities.begin(), _maturities.end(), maturity);
    
    if (strikeIt == _strikes.end()) {
        _strikes.push_back(strike);
    }
    if (maturityIt == _maturities.end()) {
        _maturities.push_back(maturity);
    }
    
    // Could be improved by maintaining a map of (strike, maturity) -> volatility ?
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::setInterpolationType(VolatilitySurface::InterpolationType type)
{
    _interpolationType = type;
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::setDiscountCurve(const DiscountCurve& discountCurve)
{ // Set the discount curve for forward moneyness interpolation
    _discountCurve = std::unique_ptr<DiscountCurve>(discountCurve.clone()); 
    // Clone the discount curve on the heap to avoid memory leaks
    // unique_ptr is a smart pointer that automatically deletes the object when it goes out of scope
    // if something goes wrong, the object will be deleted automatically - so no need to delete manually
    return *this;
}

std::unique_ptr<VolatilitySurface> VolatilitySurfaceBuilder::build()
{
    sortAndDeduplicate();       // Sort strikes/maturities and remove duplicates
    buildVolatilityMatrix();    // Build the volatility matrix
    
    if (!_discountCurve) { // Check: if no discount curve was set, use a default flat discount curve
        _discountCurve = std::make_unique<FlatDiscountCurve>(0.05);  // Default 5% rate
        // make_unique creates object FlatDiscountCurve on the heap and sets the unique_ptr to it
        // If FlatDiscountCurve throws an exception, the object will be deleted automatically - so no need to delete manually
    }

    // Create the volatility surface
    // make_unique creates object VolatilitySurface on the heap and sets the unique_ptr to it
    // If VolatilitySurface throws an exception, the object will be deleted automatically - so no need to delete manually
    return std::make_unique<VolatilitySurface>(_strikes, _maturities, _volatilities, *_discountCurve, _interpolationType);
    // here we dereference the unique_ptr to get the actual DiscountCurve object
    // returns unique_ptr to the VolatilitySurface object - the caller owns the surface
}

void VolatilitySurfaceBuilder::sortAndDeduplicate()
{
    // Sort and remove duplicates
    std::sort(_strikes.begin(), _strikes.end());
    _strikes.erase(std::unique(_strikes.begin(), _strikes.end()), _strikes.end());
    
    std::sort(_maturities.begin(), _maturities.end());
    _maturities.erase(std::unique(_maturities.begin(), _maturities.end()), _maturities.end());
}

void VolatilitySurfaceBuilder::buildVolatilityMatrix()
{
    // Initialize volatility matrix with default values
    _volatilities.resize(_maturities.size());
    for (auto& row : _volatilities) {
        row.resize(_strikes.size(), 0.2); // Default 20% volatility
    }
    
    // In practice, we buy this matrix from data provider
    // So, the actual volatility data that was set via setVolatility()
}


/**
 * DON'T DO THIS
 * DiscountCurve* _discountCurve;  // Raw pointer
 * Manual memory management
 * Easy to forget delete
 * Exception unsafe
 * Memory leaks possible
 */