//
// Created by Patrik  Liba on 05/10/2025.
//

#ifndef CPPFM_VOLATILITYSURFACE_H
#define CPPFM_VOLATILITYSURFACE_H

#include "InterpolationSchemes.h"
#include "DiscountCurve.h"
#include <vector>
#include <memory>

/**
 * VolatilitySurface class for market data interpolation
 * Section 1.1 - Implied Volatility Surface
 * Implements equation (1.4): σ*(T_i, K_j) = C^{-1}_{T_i, K_j}(C^{Mkt}(T_i, K_j))
 */
class VolatilitySurface
{
public:
    /**
     * Constructor with market data
     * @param strikes Vector of strike prices
     * @param maturities Vector of maturities (in years)
     * @param volatilities Matrix of implied volatilities [maturity][strike]
     * @param interpolationType Type of interpolation to use
     */
    enum class InterpolationType { Linear, CubicSpline };
    
    VolatilitySurface(const std::vector<double>& strikes,
                      const std::vector<double>& maturities,
                      const std::vector<std::vector<double>>& volatilities,
                      const DiscountCurve& discountCurve,
                      InterpolationType interpolationType = InterpolationType::CubicSpline);
    
    /**
     * Get implied volatility at given strike and maturity
     * @param strike Strike price
     * @param maturity Time to maturity (in years)
     * @return Implied volatility
     */
    double getImpliedVolatility(double strike, double maturity) const;
    
    /**
     * Get implied volatility using forward moneyness interpolation
     * @param strike Strike price
     * @param maturity Time to maturity (in years)
     * @param spot Current spot price
     * @return Implied volatility
     */
    double getImpliedVolatilityForwardMoneyness(double strike, double maturity, double spot) const;
    
    /**
     * Get local volatility using Dupire formula
     * @param spot Current spot price
     * @param time Current time
     * @return Local volatility
     */
    double getLocalVolatility(double spot, double time) const;
    
    /**
     * Get volatility smile for given maturity
     * @param maturity Time to maturity
     * @param strikes Vector of strikes to evaluate
     * @return Vector of implied volatilities
     */
    std::vector<double> getVolatilitySmile(double maturity, const std::vector<double>& strikes) const;
    
    /**
     * Get volatility term structure for given strike
     * @param strike Strike price
     * @param maturities Vector of maturities to evaluate
     * @return Vector of implied volatilities
     */
    std::vector<double> getVolatilityTermStructure(double strike, const std::vector<double>& maturities) const;
    
    
    /**
     * Get surface bounds
     * @return Pair of (minStrike, maxStrike) and (minMaturity, maxMaturity)
     */
    std::pair<std::pair<double, double>, std::pair<double, double>> getBounds() const;

    // Black-Sholes
    double getBlackScholesCall(double spot, double strike, double time, double volatility) const;
    double getBlackScholesPut(double spot, double strike, double time, double volatility) const;

    // Greeks
    double getBlackScholesVega(double spot, double strike, double time, double volatility) const;
    double getBlackScholesGamma(double spot, double strike, double time, double volatility) const;
    double getBlackScholesTheta(double spot, double strike, double time, double volatility) const;

     // dSigma/dT - how implied volatility changes with time
    double getImpliedVolatilityTimeDerivative(double strike, double maturity) const;

    // dSigma/dK - how implied volatility changes with strike (smile slope)
    double getImpliedVolatilityStrikeDerivative(double strike, double maturity) const;

     // d2Sigma/dK2 - smile curvature (convexity)
    double getImpliedVolatilitySecondStrikeDerivative(double strike, double maturity) const;

    // Pointer to new VolatilitySurface instance
    VolatilitySurface* clone() const;        // Clone method for polymorphic copying

    // True if surfaces are equal
    bool operator==(const VolatilitySurface& other) const;
    
    // Get raw data access
    const std::vector<double>& getStrikes() const { return _strikes; }
    const std::vector<double>& getMaturities() const { return _maturities; }
    const std::vector<std::vector<double>>& getVolatilities() const { return _volatilities; }
    const DiscountCurve& getDiscountCurve() const { return *_discountCurve; }

private:
    std::vector<double> _strikes;
    std::vector<double> _maturities;
    std::vector<std::vector<double>> _volatilities;
    InterpolationType _interpolationType;
    std::unique_ptr<DiscountCurve> _discountCurve;
    
    // Interpolation objects for different dimensions
    std::vector<std::unique_ptr<InterpolationSchemes>> _smileInterpolators;
    std::vector<std::unique_ptr<InterpolationSchemes>> _termStructureInterpolators;
    
    void initializeInterpolators(); // to create interpolators for strikes and matiruties
    double computeDupireLocalVolatility(double spot, double time) const;
    
    // Black-Scholes pricing for Dupire formula (using DiscountCurve)
    double blackScholesCall(double spot, double strike, double time, double volatility) const;
    double blackScholesPut(double spot, double strike, double time, double volatility) const;
    
    // Analytical Black-Scholes Greeks (using DiscountCurve)
    double blackScholesVega(double spot, double strike, double time, double volatility) const;
    double blackScholesGamma(double spot, double strike, double time, double volatility) const;
    double blackScholesTheta(double spot, double strike, double time, double volatility) const;
    
    
    
    // Implied volatility surface derivatives (used in Dupire formula)
    double getImpliedVolatilityDerivativeTime(double strike, double maturity) const;
    double getImpliedVolatilityDerivativeStrike(double strike, double maturity) const;
    double getImpliedVolatilitySecondDerivativeStrike(double strike, double maturity) const;
    
    
};

/**
 * Class to build the volatility surface from market data   - VolatilitySurfaceBuilder
 */
class VolatilitySurfaceBuilder
{
/**
 * TODO: We need to add the arbitrage checks here
 * As discussed Calendar Spread Arbitrage, Butterfly Arbitrage
 */
public:
    VolatilitySurfaceBuilder& addStrike(double strike);
    VolatilitySurfaceBuilder& addMaturity(double maturity);
    VolatilitySurfaceBuilder& setVolatility(double strike, double maturity, double volatility);
    VolatilitySurfaceBuilder& setInterpolationType(VolatilitySurface::InterpolationType type);
    VolatilitySurfaceBuilder& setDiscountCurve(const DiscountCurve& discountCurve);
    
    std::unique_ptr<VolatilitySurface> build();
    
private:
    std::vector<double> _strikes;
    std::vector<double> _maturities;
    std::vector<std::vector<double>> _volatilities;
    VolatilitySurface::InterpolationType _interpolationType = VolatilitySurface::InterpolationType::CubicSpline;
    std::unique_ptr<DiscountCurve> _discountCurve;
    
    void sortAndDeduplicate(); // Sort strikes/maturities and remove duplicates
    void buildVolatilityMatrix(); // Build the volatility matrix
};

#endif //CPPFM_VOLATILITYSURFACE_H