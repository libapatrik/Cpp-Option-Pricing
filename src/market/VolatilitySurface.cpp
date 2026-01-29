//
// Created by Patrik  Liba on 05/10/2025.
//

#include <cppfm/market/VolatilitySurface.h>
#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cppfm/utils/Utils.h>
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
                                     SmileInterpolationType smileInterpolationType,             // smile interpolation type
                                     MaturityInterpolationType maturityInterpolationType)       // maturity interpolation type
    : _strikes(strikes), _maturities(maturities), _volatilities(volatilities),
      _smileInterpolationType(smileInterpolationType), 
      _maturityInterpolationType(maturityInterpolationType),
      _discountCurve(discountCurve.clone())
{
    validateInputData();      // validate input data
    initializeInterpolators(); // proceeds to next step - call private method to create interpolators for strikes and maturities
}

void VolatilitySurface::validateInputData() const
{
    // Check for empty data
    if (_strikes.empty() || _maturities.empty() || _volatilities.empty()) {
        throw std::invalid_argument("VolatilitySurface: Empty data provided");
    }
    
    // Check matrix dimensions
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
}
/**
 *   Smile Interpolation: For each maturity, interpolate volatility across strikes (the volatility smile)
 *   Term Structure Interpolation: For each strike, interpolate volatility across maturities (how volatility changes over time)
 */
void VolatilitySurface::initializeInterpolators()
{
    // Initialize 1D interpolators for each maturity (smile interpolation)
    // cubic spline interpolation along strikes
    _smileInterpolators.reserve(_maturities.size());                     // reserve space for maturities in _smileInterpolators vector
    for (size_t i = 0; i < _maturities.size(); ++i) {                    // for each maturity
        if (_smileInterpolationType == SmileInterpolationType::CubicSpline) {      // check type: do cubic spline interpolation
            _smileInterpolators.push_back(
                std::make_unique<CubicSplineInterpolation>(              // creates a unique pointer to a new CubicSplineInterpolation object
                    _strikes, _volatilities[i], 
                    CubicSplineInterpolation::BoundaryType::Natural));   // boundary type: natural
        } else {
            // Linear interpolation for smiles
            _smileInterpolators.push_back(
                std::make_unique<LinearInterpolation>(_strikes, _volatilities[i]));
        }
    }
    
    // Initialize 1D interpolators for each strike (term structure interpolation)
    _termStructureInterpolators.reserve(_strikes.size());            // reserve space for strikes
    for (size_t j = 0; j < _strikes.size(); ++j) {                   // for each strike
        std::vector<double> termVols;                                // volatilities across maturities for strike j
        termVols.reserve(_maturities.size());                        // reserve space for maturities
        for (size_t i = 0; i < _maturities.size(); ++i) {            // for each maturity
            termVols.push_back(_volatilities[i][j]);                 // collect volatilities for strike j
        }
        
        if (_smileInterpolationType == SmileInterpolationType::CubicSpline) {
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
    // Dispatch to appropriate interpolation method based on _maturityInterpolationType
    switch (_maturityInterpolationType) {
        case MaturityInterpolationType::Bilinear:
            return impliedVolatilityBilinear(strike, maturity);
        case MaturityInterpolationType::ForwardMoneyness:
            return impliedVolatilityForwardMoneyness(strike, maturity);
        default:
            throw std::runtime_error("Unknown maturity interpolation type");
    }
}

double VolatilitySurface::impliedVolatilityBilinear(double strike, double maturity) const
{
    // Simple bilinear interpolation: interpolate strikes → interpolate time
    // Check bounds
    if (strike < _strikes.front() || strike > _strikes.back() ||
        maturity < _maturities.front() || maturity > _maturities.back()) {
        throw std::out_of_range("VolatilitySurface: Point outside interpolation range");
    }
    
    // Find the closest maturity index
    auto maturityIt = std::lower_bound(_maturities.begin(), _maturities.end(), maturity); // first element >=
    size_t maturityIdx = std::distance(_maturities.begin(), maturityIt);
    
    if (maturityIdx == 0) {
        return (*_smileInterpolators[0])(strike); // Use operator() for automatic extrapolation handling
    }
    
    if (maturityIdx >= _maturities.size()) {
        return (*_smileInterpolators.back())(strike); // Use operator() for automatic extrapolation handling
    }
    
    // Interpolate between adjacent maturities
    double t1 = _maturities[maturityIdx - 1];
    double t2 = _maturities[maturityIdx];
    double vol1 = (*_smileInterpolators[maturityIdx - 1])(strike); // Use operator() for automatic extrapolation handling
    double vol2 = (*_smileInterpolators[maturityIdx])(strike); // Use operator() for automatic extrapolation handling
    
    // Linear interpolation in volatility space
    double weight = (maturity - t1) / (t2 - t1);
    return vol1 + weight * (vol2 - vol1);
}

double VolatilitySurface::impliedVolatilityForwardMoneyness(double strike, double maturity) const
{
    // Section 1.3 - Interpolation along maturities using constant forward moneyness
    // Implements equation (1.27) from lecture notes
    //
    // Algorithm: For a given maturity T ∈ [T_i, T_{i+1}] and strike K:
    // 1. Compute forward moneyness: k_F_T = K/F_T = (K/S_0) * e^{-∫_0^T r(s)ds}
    // 2. Extract strikes K^(i) at pillar maturities with same forward moneyness:
    //    K^(i) = k_F_T × F_{T_i} = K × B(T)/B(T_i)
    //    where B(t) = e^{-∫_0^t r(s)ds} is the discount factor
    // 3. Interpolate in variance space: v(T,k) = (σ*)^2(T, K^(i))T_i
    //
    // Key insight: S_0 cancels out in K^(i) = K × B(T)/B(T_i), so spot not needed!
    
    // Discount factor at target maturity T
    double B_T = _discountCurve->discount(maturity);
    
    // Find maturity interval
    auto maturityIt = std::lower_bound(_maturities.begin(), _maturities.end(), maturity);
    size_t maturityIdx = std::distance(_maturities.begin(), maturityIt);

    if (maturityIdx == 0) {
        // Extrapolation: T < T_1
        double t1 = _maturities[0];
        double B_t1 = _discountCurve->discount(t1);
        // K^(1) = K × B(T)/B(T_1) - strike at T_1 with same forward moneyness
        double strike1 = strike * B_T / B_t1;
        double vol1 = (*_smileInterpolators[0])(strike1); // Use operator() for automatic extrapolation handling
        return vol1;
    }
    
    if (maturityIdx >= _maturities.size()) {
        // Extrapolation: T > T_M
        double tM = _maturities.back();
        double B_tM = _discountCurve->discount(tM);
        // K^(M) = K × B(T)/B(T_M)
        double strikeM = strike * B_T / B_tM;
        double volM = (*_smileInterpolators.back())(strikeM); // Use operator() for automatic extrapolation handling
        return volM;
    }

    // Interpolation between T_i and T_{i+1}
    double t1 = _maturities[maturityIdx - 1];
    double t2 = _maturities[maturityIdx];
    
    // Discount factors at pillar maturities
    double B_t1 = _discountCurve->discount(t1);
    double B_t2 = _discountCurve->discount(t2);
    
    // Extract strikes corresponding to constant forward moneyness
    // K^(i) = K × B(T)/B(T_i) and K^(i+1) = K × B(T)/B(T_{i+1})
    double strike1 = strike * B_T / B_t1;
    double strike2 = strike * B_T / B_t2;
    
    // Get implied volatilities at these strikes from smile interpolators
    // Use operator() for automatic extrapolation handling if strikes are outside range
    double vol1 = (*_smileInterpolators[maturityIdx - 1])(strike1);
    double vol2 = (*_smileInterpolators[maturityIdx])(strike2);
    
    // Linear interpolation in variance space
    // Variance quantity: v(T,k) = (σ*)^2(T, K^(i)) × T
    double variance1 = vol1 * vol1 * t1;
    double variance2 = vol2 * vol2 * t2;
    
    // Linearly interpolate variance
    double weight = (maturity - t1) / (t2 - t1);
    double interpolatedVariance = variance1 + weight * (variance2 - variance1);
    
    // Convert back to volatility: σ*(T,K) = sqrt(v(T,k_F_T) / T)
    // Equation (1.27) from lecture notes
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
    // Ensure no division by zero
    constexpr double MIN_TIME = 1e-12;
    if (time <= MIN_TIME) {
        return impliedVol; // for very small times, return implied volatility
    }

    // For times before the first maturity pillar, Dupire derivatives are unreliable
    // Use flat extrapolation from the first reliable time point
    double firstMaturity = _maturities.front();
    double minReliableTime = firstMaturity * 0.5;  // Require at least half of first pillar
    if (time < minReliableTime && firstMaturity > MIN_TIME) {
        // Recursively compute local vol at the first reliable time
        // and use flat extrapolation backwards
        return computeDupireLocalVolatility(spot, minReliableTime);
    }
    
    double volSqrtTime = impliedVol * std::sqrt(time);
    double d1 = (std::log(spot / strike) + (riskFreeRate + 0.5 * impliedVol * impliedVol) * time) / volSqrtTime;
    double d2 = d1 - volSqrtTime;
    
    // Compute the numerator: 1 + 2T/σ* × (∂σ*/∂T + r×K×∂σ*/∂K)
    // This captures how the implied volatility surface changes over time and across strikes
    // NOTE: The rate term requires multiplication by strike K (from Dupire formula derivation)
    double numerator = 1.0 + (2.0 * time / impliedVol) * (dSigma_dT + riskFreeRate * strike * dSigma_dK);
    
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
    
    // Check for non-finite values
    if (!std::isfinite(numerator) || !std::isfinite(denominator)) {
        throw std::runtime_error("Dupire: Non-finite values detected at K=" + 
                                std::to_string(strike) + ", T=" + std::to_string(time));
    }
    // TODO: Add arbitrage checks for the surface + tests
    // Check for negative denominator (potential arbitrage)
    if (denominator <= 0.0) {
        // WARNING: This indicates potential arbitrage in the implied volatility surface
        // Apply correction: use absolute value and continue with warning
        denominator = std::abs(denominator);
    }

    // Check for negative numerator
    if (numerator < 0.0) {
        // Apply correction: use absolute value
        numerator = std::abs(numerator);
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
     * Uses centered finite difference approximation via NumericalDerivatives utility.
     *
     * A positive value means implied volatility increases with time (smile steepens).
     * A negative value means implied volatility decreases with time (smile flattens).
     */

    // Create lambda: fix strike K, differentiate with respect to time T
    auto volAsTimeFunction = [this, strike](double t) {
        return this->impliedVolatility(strike, t);
    };

    // Step size: 1% relative to maturity, but at least 0.01 years (about 3-4 trading days)
    double h = std::max(0.01 * maturity, 0.01);
    return NumericalDerivatives::firstDerivative(volAsTimeFunction, maturity, h);
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
     * Uses centered finite difference approximation via NumericalDerivatives utility.
     *
     * A positive value means implied volatility increases with strike (upward sloping smile).
     * A negative value means implied volatility decreases with strike (downward sloping smile).
     * Near zero indicates a flat smile at that point.
     */

    // Create lambda: fix maturity T, differentiate with respect to strike K
    auto volAsStrikeFunction = [this, maturity](double k) {
        return this->impliedVolatility(k, maturity);
    };

    // Vol-adaptive step size: scale with implied vol and sqrt(T)
    // This is approximately 2% of the ATM straddle width
    double vol = impliedVolatility(strike, maturity);
    double h = 0.02 * strike * vol * std::sqrt(maturity);
    h = std::max(h, 0.005 * strike);  // Floor at 0.5% of strike for very low vol/short T
    return NumericalDerivatives::firstDerivative(volAsStrikeFunction, strike, h);
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
     * Uses centered second difference via NumericalDerivatives utility.
     *
     * A positive value indicates convexity (smile curves upward like a parabola).
     * A negative value indicates concavity (smile curves downward).
     * Near zero indicates the smile is approximately linear at that point.
     *
     * In Dupire's formula this accounts for smile curvature effects.
     *
     * Note: Second derivative is numerically sensitive (h^2 in denominator).
     * Step size scales with implied vol and sqrt(T) for stability across vol regimes.
     */

    // Create lambda: fix maturity T, differentiate twice with respect to strike K
    auto volAsStrikeFunction = [this, maturity](double k) {
        return this->impliedVolatility(k, maturity);
    };

    // Vol-adaptive step size: ~5% of ATM straddle width
    // Higher vol -> smoother prices -> need larger step to capture curvature
    double vol = impliedVolatility(strike, maturity);
    double h = 0.05 * strike * vol * std::sqrt(maturity);
    h = std::max(h, 0.01 * strike);  // Floor at 1% of strike for very low vol/short T
    return NumericalDerivatives::secondDerivative(volAsStrikeFunction, strike, h);
}

std::pair<std::pair<double, double>, std::pair<double, double>> VolatilitySurface::getBounds() const
{
    return {{_strikes.front(), _strikes.back()}, 
            {_maturities.front(), _maturities.back()}};
}

std::unique_ptr<VolatilitySurface> VolatilitySurface::clone() const
{
    return std::make_unique<VolatilitySurface>(_strikes, _maturities, _volatilities, 
                                 *_discountCurve, 
                                 _smileInterpolationType, 
                                 _maturityInterpolationType);
}

bool VolatilitySurface::operator==(const VolatilitySurface& other) const
{
    return _strikes == other._strikes &&
           _maturities == other._maturities &&
           _volatilities == other._volatilities &&
           _smileInterpolationType == other._smileInterpolationType &&
           _maturityInterpolationType == other._maturityInterpolationType;
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
    // Store strike and maturity if not already present
    auto strikeIt = std::find(_strikes.begin(), _strikes.end(), strike);
    auto maturityIt = std::find(_maturities.begin(), _maturities.end(), maturity);
    
    if (strikeIt == _strikes.end()) {
        _strikes.push_back(strike);
    }
    if (maturityIt == _maturities.end()) {
        _maturities.push_back(maturity);
    }
    
    // Store volatility in map: (strike, maturity) -> volatility
    _volatilityMap[{strike, maturity}] = volatility;
    
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::setSmileInterpolationType(VolatilitySurface::SmileInterpolationType type)
{
    _smileInterpolationType = type;
    return *this;
}

VolatilitySurfaceBuilder& VolatilitySurfaceBuilder::setMaturityInterpolationType(VolatilitySurface::MaturityInterpolationType type)
{
    _maturityInterpolationType = type;
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
    return std::make_unique<VolatilitySurface>(_strikes, _maturities, _volatilities, *_discountCurve, 
                                               _smileInterpolationType, _maturityInterpolationType);
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
    
    // Fill matrix from stored volatility map
    // Matrix indexing: _volatilities[maturity_index][strike_index]
    for (size_t i = 0; i < _maturities.size(); ++i) {
        for (size_t j = 0; j < _strikes.size(); ++j) {
            auto key = std::make_pair(_strikes[j], _maturities[i]);
            auto it = _volatilityMap.find(key);
            if (it != _volatilityMap.end()) {
                _volatilities[i][j] = it->second; // Use stored volatility
            }
            // Otherwise keep default value (0.2)
        }
    }
}