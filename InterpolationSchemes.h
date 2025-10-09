//
// Created by Patrik  Liba on 05/10/2025.
//

#ifndef CPPFM_INTERPOLATIONSCHEMES_H
#define CPPFM_INTERPOLATIONSCHEMES_H

#include <vector>
#include <memory>

/**
 * Abstract base class for interpolation schemes
 */
class InterpolationSchemes
{
public:
    virtual ~InterpolationSchemes() = default;

    virtual double interpolate(double x) const = 0;

    virtual InterpolationSchemes* clone() const = 0; // Pointer to new interpolator instance

    // Pair of (min, max) x values
    virtual std::pair<double, double> getRange() const = 0;

    // Extrapolated value
    virtual double extrapolate(double x) const = 0;
};

/**
 * Linear interpolation implementation
 * Between two points (x0, y0) and (x1, y1) the interpolation is a linear function
 * y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
 */
class LinearInterpolation : public InterpolationSchemes
{
public:
    /**
     * Constructor with data points
     * @param xData X coordinates (must be sorted)
     * @param yData Y coordinates
     */
    LinearInterpolation(const std::vector<double>& xData, const std::vector<double>& yData);
    
    double interpolate(double x) const override;
    InterpolationSchemes* clone() const override;
    std::pair<double, double> getRange() const override;
    double extrapolate(double x) const override;
    
private:
    std::vector<double> _xData, _yData;
    double _xMin, _xMax;
    
    size_t findInterval(double x) const;
};

/**
 * Cubic spline interpolation with Thomas algorithm
 * For each interval [x0, x1], the interpolation is a cubic polynomial
 * S(x) = a0 + a1 * (x - x0) + a2 * (x - x0)^2 + a3 * (x - x0)^3
 * coeffs found by solving the tridiagonal system of equations
 * the system is solved using the Thomas algorithm
 */
class CubicSplineInterpolation : public InterpolationSchemes
{
public:
    /**
     * Constructor with data points and boundary conditions
     * @param xData X coordinates (must be sorted)
     * @param yData Y coordinates
     * @param boundaryType Boundary condition type
     */
    enum class BoundaryType { Natural, Clamped, Periodic };  // for future direction
    
    CubicSplineInterpolation(const std::vector<double>& xData, 
                            const std::vector<double>& yData,
                            BoundaryType boundaryType = BoundaryType::Natural);
    
    double interpolate(double x) const override;
    InterpolationSchemes* clone() const override;
    std::pair<double, double> getRange() const override;
    double extrapolate(double x) const override;

    double derivative(double x) const;

    double secondDerivative(double x) const;
    
private:
    std::vector<double> _xData, _yData;
    double _xMin, _xMax;
    BoundaryType _boundaryType;
    
    // LECTURE NOTES: Store pre-computed coefficients as per equations (1.26)-(1.27)
    std::vector<double> _alpha;  // α coefficients
    std::vector<double> _beta;   // β coefficients  
    std::vector<double> _gamma;  // γ coefficients
    std::vector<double> _delta; // δ coefficients
    
    void computeSplineCoefficients();
    void solveThomasAlgorithm();
    size_t findInterval(double x) const;

};


#endif //CPPFM_INTERPOLATIONSCHEMES_H