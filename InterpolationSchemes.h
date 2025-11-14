//
// Created by Patrik  Liba on 08/11/2025.
//

#ifndef CPPFM_INTERPOLATIONSCHEMES_H
#define CPPFM_INTERPOLATIONSCHEMES_H

#include<vector>


/**
 * Abstract base class for interpolation schemes
 *
 * DESIGN:
 * - Interpolation provides the analytical form within the data range
 * - Extrapolation is generic and uses derivatives at boundaries
 * - All schemes must provide analytical derivatives
 */

class ExtrapolationScheme;

class InterpolationSchemes
{
public:
    // constructor
    InterpolationSchemes(const std::vector<double>& xData, const std::vector<double>& yData, double xMin, double xMax, const ExtrapolationScheme& extrapolationScheme);

    virtual ~InterpolationSchemes() = default;
    virtual InterpolationSchemes* clone() const = 0;

    virtual double interpolate(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual double secondDerivative(double x) const = 0;
    // Pair of (min, max) x values
    virtual std::pair<double, double> getRange() const = 0;


    /**
     * Generic extrapolation using linear continuation based on boundary derivatives
     * This method is generic and doesn't depend on the specific interpolation scheme
     * It uses the analytical derivative at the boundary for linear extrapolation
     * @param x Point to extrapolate at (outside range)
     * @return Extrapolated value
     */

    // Keep .h methods as only declared. Define them in cpp file.
    virtual double extrapolate(double x) const
    {
        auto [xMin, xMax] = getRange();
        if (x < xMin) {
            // Left extrapolation: y(x) = y(xMin) + y'(xMin) * (x - xMin)
            double yMin = interpolate(xMin);
            double dy = derivative(xMin);
            double d2y = secondDerivative(xMin);
            double dx = x - xMin;
            // Taylor series expansion
            return yMin + dy * dx + 0.5 * d2y * dx * dx;      // Quadratic extrapolation
        } else {
            // Right extrapolation: y(x) = y(xMax) + y'(xMax) * (x - xMax)
            double yMax = interpolate(xMax);
            double dy = derivative(xMax);
            double d2y = secondDerivative(xMax);
            double dx = x - xMax;
            // Taylor series expansion
            return yMax + dy * dx + 0.5 * d2y * dx * dx;
        }
    }
    /**
     * Unified interface: interpolate or extrapolate as needed
     * @param x Point to evaluate at
     * @return Value at x
     */
    double operator()(double x) const // to handle routing to extrapolate() when needed
    {
        auto [xMin, xMax] = getRange();
        if (x < xMin || x > xMax) {
            return _extrapolationScheme->extrapolate(x);
        } else {
            return interpolate(x);
        }
    }
    // Class instance(param1, param2,...) // Constructor
    // LinearInterpolationScheme interp(params); // Constructor
    // interp(x) -> interp.operator()(x);

    // operator[](size_t row, size_t column)
    // Matrix A;
    // A[2,3]

    // operator=(double x)
    // double a;
    // double x;
    // a = x -> a.operator=(x)

protected:
    std::vector<double> _xData, _yData;
    double _xMin, _xMax;
    std::unique_ptr<ExtrapolationScheme> _extrapolationScheme;



}; // end of base class InterpolationSchemes1

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
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    InterpolationSchemes* clone() const override;
    std::pair<double, double> getRange() const override;

    // Using default generic extrapolation from base class
    // No need to override extrapolate()
private:
    std::vector<double> _xData, _yData;
    double _xMin, _xMax;

    size_t findInterval(double x) const;



}; // end of class LinearInterpolation


/**
 * Cubic spline interpolation with Thomas algorithm
 * Analytical form: S(x) = α(x-x0)³ + β(x-x0)² + γ(x-x0) + δ
 * First derivative: S'(x) = 3α(x-x0)² + 2β(x-x0) + γ
 * Second derivative: S''(x) = 6α(x-x0) + 2β
 */
class CubicSplineInterpolation : public InterpolationSchemes
{
public:
    enum class BoundaryType { Natural, Clamped, Periodic }; // for future direction
    CubicSplineInterpolation(const std::vector<double>& xData,
                             const std::vector<double>& yData,
                             BoundaryType boundaryType = BoundaryType::Natural);

    double interpolate(double x) const override;
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    InterpolationSchemes* clone() const override;
    std::pair<double, double> getRange() const override;

    // Using default generic extrapolation from base class
    // No need to override extrapolate()

private:
    std::vector<double> _xData, _yData;
    double _xMin, _xMax;
    BoundaryType _boundaryType;

    // Pre-computed coefficients from lecture notes (1.26-7)
    std::vector<double> _alpha;
    std::vector<double> _beta;
    std::vector<double> _gamma;
    std::vector<double> _delta;

    void computeSplineCoefficients();
    void solveThomasAlgorithm();
    size_t findInterval(double x) const;
};  // end of class CubicSplineInterpolation


// Abstract class
class ExtrapolationScheme
{
public:
    // Constructor
    virtual double extrapolate(double x) const = 0;

protected:
    double _xMin, _xMax;
};

class FlatExtrapolationScheme : public ExtrapolationScheme
{
public:
    double extrapolate(double x) const override;
};

class LinearExtrapolationScheme : public ExtrapolationScheme
{
public:
    double extrapolate(double x) const override;
private:
    double _dxMin, _dxMax;
};

class QuadraticExtrapolationScheme : public ExtrapolationScheme
{
public:
    double extrapolate(double x) const override;
private:
    double _dxMin, _dxMax;
    double _d2xMin, _d2xMax;
};

#endif //CPPFM_INTERPOLATIONSCHEMES_H
