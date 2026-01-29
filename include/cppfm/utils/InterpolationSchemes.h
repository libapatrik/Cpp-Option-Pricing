//
// Created by Patrik  Liba on 08/11/2025.
//

#ifndef CPPFM_INTERPOLATIONSCHEMES_H
#define CPPFM_INTERPOLATIONSCHEMES_H

#include <vector>
#include <memory>
#include <utility>


/**
 *  InterpolationScheme has a constructor that only require which extrapolation you want [ExtrapolationType]
 *  When calling the extrapolation() method, the data member _extrapolationType will determine how to extrapolation
 *  You need one protected method per extrapolation type
 */

// Extrapolation type enum
enum class ExtrapolationType { Flat, Linear, Quadratic };

// Forward declaration
class ExtrapolationScheme;

/**
 * Abstract base class for interpolation schemes
 *
 * DESIGN:
 * - Interpolation provides the analytical form within the data range
 * - Extrapolation is generic and uses derivatives at boundaries
 */

// ============================================================================
// BASE CLASS: InterpolationScheme
// ============================================================================
class InterpolationScheme
{
public:
    // Constructor with extrapolation type 
    InterpolationScheme(const std::vector<double>& xData, 
                        const std::vector<double>& yData,
                        ExtrapolationType extraType = ExtrapolationType::Quadratic);
    virtual ~InterpolationScheme() = default;
    
    // Core interface
    virtual double interpolate(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual double secondDerivative(double x) const = 0;
    virtual std::unique_ptr<InterpolationScheme> clone() const = 0;
    
    double operator()(double x) const; // Unified interface: routes to interpolate or extrapolate automatically
    std::pair<double, double> getRange() const; // Common implementation for all derived classes

protected:
    std::vector<double> _xData;
    std::vector<double> _yData;
    std::unique_ptr<ExtrapolationScheme> _extrapolationScheme;
    
    ExtrapolationType _extrapolationType;
    
    // Common validation helper
    void validateData() const;
    
    // Common interval finding using binary search
    size_t findInterval(double x) const;
};



/**
 * Linear interpolation: y = y0 + (y1-y0) * (x-x0) / (x1-x0)
 */
class LinearInterpolation : public InterpolationScheme
{
public:
    // Constructor with extrapolation type (default: Quadratic for backward compatibility)
    LinearInterpolation(const std::vector<double>& xData, 
                        const std::vector<double>& yData,
                        ExtrapolationType extraType = ExtrapolationType::Quadratic);

    double interpolate(double x) const override;
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    std::unique_ptr<InterpolationScheme> clone() const override;
};

// ============================================================================
// CUBIC SPLINE INTERPOLATION
// ============================================================================

/**
 * Cubic spline interpolation with Thomas algorithm
 * Analytical form: S(x) = α(x-x0)³ + β(x-x0)² + γ(x-x0) + δ
 * First derivative: S'(x) = 3α(x-x0)² + 2β(x-x0) + γ
 * Second derivative: S''(x) = 6α(x-x0) + 2β
 */
class CubicSplineInterpolation : public InterpolationScheme
{
public:
    enum class BoundaryType { Natural, Clamped, Periodic }; // for future direction
    
    // Constructor with extrapolation type (default: Quadratic for backward compatibility)
    CubicSplineInterpolation(const std::vector<double>& xData,
                             const std::vector<double>& yData,
                             BoundaryType boundaryType = BoundaryType::Natural,
                             ExtrapolationType extraType = ExtrapolationType::Quadratic);

    double interpolate(double x) const override;
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    std::unique_ptr<InterpolationScheme> clone() const override;

private:
    BoundaryType _boundaryType;

    // Pre-computed coefficients from lecture notes (1.26-7)
    std::vector<double> _alpha;
    std::vector<double> _beta;
    std::vector<double> _gamma;
    std::vector<double> _delta;

    void computeSplineCoefficients();
    void solveThomasAlgorithm();
};


// ============================================================================
// BASE CLASS: ExtrapolationScheme
// ============================================================================
class ExtrapolationScheme
{
public:
    virtual ~ExtrapolationScheme() = default;
    virtual std::unique_ptr<ExtrapolationScheme> clone() const = 0;
    virtual void initialize(const InterpolationScheme& interp) = 0;
    virtual double extrapolate(double x, const InterpolationScheme& interp) const = 0;
};

// ============================================================================
// EXTRAPOLATION SCHEMES
// ============================================================================

// Flat extrapolation: y = y(boundary)
class FlatExtrapolation : public ExtrapolationScheme
{
public:
    std::unique_ptr<ExtrapolationScheme> clone() const override;
    void initialize(const InterpolationScheme& interp) override; //extract the boundary information from the interpolation scheme
    double extrapolate(double x, const InterpolationScheme& interp) const override;

private:
    double _yMin, _yMax;
};

// Linear extrapolation: y = y(boundary) + y'(boundary) * dx
class LinearExtrapolation : public ExtrapolationScheme
{
public:
    std::unique_ptr<ExtrapolationScheme> clone() const override;
    void initialize(const InterpolationScheme& interp) override;
    double extrapolate(double x, const InterpolationScheme& interp) const override;  // avoid circular reference!

private:
    double _yMin, _yMax;
    double _dyMin, _dyMax;
};

// Quadratic extrapolation: Taylor series to 2nd order (default)
class QuadraticExtrapolation : public ExtrapolationScheme
{
public:
    std::unique_ptr<ExtrapolationScheme> clone() const override;
    void initialize(const InterpolationScheme& interp) override;
    double extrapolate(double x, const InterpolationScheme& interp) const override;

private:
    double _yMin, _yMax;
    double _dyMin, _dyMax;
    double _d2yMin, _d2yMax;
};

#endif //CPPFM_INTERPOLATIONSCHEMES_H
