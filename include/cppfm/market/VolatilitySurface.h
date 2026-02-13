//
// Created by Patrik  Liba on 05/10/2025.
//

#ifndef CPPFM_VOLATILITYSURFACE_H
#define CPPFM_VOLATILITYSURFACE_H

#include <cppfm/utils/InterpolationSchemes.h>
#include <cppfm/market/DiscountCurve.h>
#include <vector>
#include <memory>
#include <map>
#include <utility>
#include <string>

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
     * @param discountCurve Discount curve for forward price calculations
     * @param smileInterpolationType Smile interpolation (across strikes): Linear or CubicSpline
     * @param maturityInterpolationType Maturity interpolation: Bilinear or ForwardMoneyness
     */
    enum class SmileInterpolationType { Linear, CubicSpline };
    enum class MaturityInterpolationType { Bilinear,          /// TODO: Fix naming: Simple linear interpolation in volatility space
                                           ForwardMoneyness   // Variance interpolation with constant forward moneyness (recommended)
    };
    
    VolatilitySurface(const std::vector<double>& strikes,
                      const std::vector<double>& maturities,
                      const std::vector<std::vector<double>>& volatilities,
                      const DiscountCurve& discountCurve,
                      SmileInterpolationType smileInterpolationType = SmileInterpolationType::CubicSpline,
                      MaturityInterpolationType maturityInterpolationType = MaturityInterpolationType::ForwardMoneyness);
    
    /**
     * Get implied volatility at given strike and maturity
     * Uses the interpolation method specified in constructor
     * @param strike Strike price K
     * @param maturity Time to maturity T (in years)
     * @return Implied volatility σ*(T,K)
     */
    double impliedVolatility(double strike, double maturity) const;
    
    /**
     * Get local volatility using Dupire formula
     * Delegates to computeDupireLocalVolatility (private implementation)
     * @param spot Current spot price
     * @param time Current time
     * @return Local volatility σ_local(S,t)
     */
    double localVolatility(double spot, double time) const;
    
    /**
     * Get surface bounds
     * @return Pair of (minStrike, maxStrike) and (minMaturity, maxMaturity)
     */
    std::pair<std::pair<double, double>, std::pair<double, double>> getBounds() const;

    // Pointer to new VolatilitySurface instance
    std::unique_ptr<VolatilitySurface> clone() const;        // Clone method for polymorphic copying

    // True if surfaces are equal
    bool operator==(const VolatilitySurface& other) const;
    
    // Raw data accessors
    const std::vector<double>& strikes() const { return _strikes; } 
    const std::vector<double>& maturities() const { return _maturities; }
    const std::vector<std::vector<double>>& volatilities() const { return _volatilities; }
    const DiscountCurve& discountCurve() const { return *_discountCurve; }

    friend class VolatilitySurfaceBuilder;

private:
    std::vector<double> _strikes;
    std::vector<double> _maturities;
    std::vector<std::vector<double>> _volatilities;
    SmileInterpolationType _smileInterpolationType;
    MaturityInterpolationType _maturityInterpolationType;
    std::unique_ptr<DiscountCurve> _discountCurve;
    double _initialSpot = 0.0;  // S_0 for Dupire d1; 0 = use strike as spot (backward compat)

    // Interpolation objects for different dimensions
    std::vector<std::unique_ptr<InterpolationScheme>> _smileInterpolators;            // 1 per maturity interpolate across strikes
    std::vector<std::unique_ptr<InterpolationScheme>> _termStructureInterpolators;    // 1 per strike interpolate acrsso maturities
    // An "is-a" relation does not apply to unique_ptr; only for class-to-class

    void validateInputData() const; // validate strikes, maturities, volatilities matrix
    void initializeInterpolators(); // to create interpolators for strikes and maturities
    double computeDupireLocalVolatility(double spot, double time) const;
    
    // Private interpolation methods - called by impliedVolatility based on _maturityInterpolationType
    double impliedVolatilityBilinear(double strike, double maturity) const;
    double impliedVolatilityForwardMoneyness(double strike, double maturity) const;
    
    // Implied volatility surface derivatives (used in Dupire formula) - private helpers
    double impliedVolatilityDerivativeTime(double strike, double maturity) const;
    double impliedVolatilityDerivativeStrike(double strike, double maturity) const;
    double impliedVolatilitySecondDerivativeStrike(double strike, double maturity) const;
    
    
};

// diagnostics from arbitrage checks in the builder
struct ArbitrageDiagnostic {
    int maturityIdx;       // which maturity slice
    int strikeIdx;         // which strike (-1 for calendar = all strikes at that pair)
    double value;          // offending value (density for butterfly, variance diff for calendar)
    std::string type;      // "butterfly" or "calendar"
};

/**
 * Class to build the volatility surface from market data   - VolatilitySurfaceBuilder
 */
class VolatilitySurfaceBuilder
{
public:
    VolatilitySurfaceBuilder& addStrike(double strike);
    VolatilitySurfaceBuilder& addMaturity(double maturity);
    VolatilitySurfaceBuilder& setVolatility(double strike, double maturity, double volatility);    // store volatility in a map
    VolatilitySurfaceBuilder& setSmileInterpolationType(VolatilitySurface::SmileInterpolationType type);
    VolatilitySurfaceBuilder& setMaturityInterpolationType(VolatilitySurface::MaturityInterpolationType type);
    VolatilitySurfaceBuilder& setDiscountCurve(const DiscountCurve& discountCurve);
    VolatilitySurfaceBuilder& setInitialSpot(double spot);

    std::unique_ptr<VolatilitySurface> build();

    const std::vector<ArbitrageDiagnostic>& diagnostics() const { return _diagnostics; }

private:
    std::vector<ArbitrageDiagnostic> _diagnostics;
    bool checkButterflyArbitrage();   // returns true if clean
    bool checkCalendarArbitrage();    // returns true if clean

    std::vector<double> _strikes;
    std::vector<double> _maturities;
    std::vector<std::vector<double>> _volatilities;
    std::map<std::pair<double, double>, double> _volatilityMap; // member to store (strike, maturity) -> volatility pairs
    VolatilitySurface::SmileInterpolationType _smileInterpolationType = VolatilitySurface::SmileInterpolationType::CubicSpline;
    VolatilitySurface::MaturityInterpolationType _maturityInterpolationType = VolatilitySurface::MaturityInterpolationType::ForwardMoneyness;
    std::unique_ptr<DiscountCurve> _discountCurve;
    double _initialSpot = 0.0;

    void sortAndDeduplicate(); // Sort strikes/maturities and remove duplicates
    void buildVolatilityMatrix(); // Build the volatility matrix
};

#endif //CPPFM_VOLATILITYSURFACE_H