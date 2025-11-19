//
// Created by Patrik  Liba on 08/11/2025.
//

#ifndef CPPFM_INTERPOLATIONSCHEMES_H
#define CPPFM_INTERPOLATIONSCHEMES_H

#include <vector>
#include <memory>
#include <utility>



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
    virtual ~InterpolationScheme() = default;
    
    // Core interface
    virtual double interpolate(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual double secondDerivative(double x) const = 0;
    // Pair of (min, max) x values
    virtual std::pair<double, double> getRange() const = 0;
    virtual std::unique_ptr<InterpolationScheme> clone() const = 0;
    
    // Extrapolation strategy, we use setter to set the choice of the extrapolation scheme
    // Default: QuadraticExtrapolation
    /**
     *  auto interp = std::make_unique<LinearInterpolation>(xData, yData, std::move(scheme)); // user specifies the scheme
     *  transfers the ownership of the unique_ptr from the caller to _extrapolationScheme without copying. 
     */
    void setExtrapolationScheme(std::unique_ptr<ExtrapolationScheme> scheme);  
    

    double operator()(double x) const; // Unified interface: routes to interpolate or extrapolate automatically

    
protected:
    std::unique_ptr<ExtrapolationScheme> _extrapolationScheme;
};



/**
 * Linear interpolation: y = y0 + (y1-y0) * (x-x0) / (x1-x0)
 */
class LinearInterpolation : public InterpolationScheme
{
public:
    LinearInterpolation(const std::vector<double>& xData, const std::vector<double>& yData);
    
    // Constructor with extrapolation choice
    LinearInterpolation(const std::vector<double>& xData, 
                       const std::vector<double>& yData,
                       std::unique_ptr<ExtrapolationScheme> extrapolationScheme);

    double interpolate(double x) const override;
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    std::unique_ptr<InterpolationScheme> clone() const override;
    std::pair<double, double> getRange() const override;

private:
    std::vector<double> _xData;
    std::vector<double> _yData;
    double _xMin;
    double _xMax;
    
    size_t findInterval(double x) const;
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
    CubicSplineInterpolation(const std::vector<double>& xData,
                             const std::vector<double>& yData,
                             BoundaryType boundaryType = BoundaryType::Natural);
    
    // Constructor with extrapolation choice
    CubicSplineInterpolation(const std::vector<double>& xData,
                             const std::vector<double>& yData,
                             BoundaryType boundaryType,
                             std::unique_ptr<ExtrapolationScheme> extrapolationScheme);

    double interpolate(double x) const override;
    double derivative(double x) const override;
    double secondDerivative(double x) const override;
    std::unique_ptr<InterpolationScheme> clone() const override;
    std::pair<double, double> getRange() const override;



private:
    std::vector<double> _xData;
    std::vector<double> _yData;
    double _xMin;
    double _xMax;
    BoundaryType _boundaryType;

    // Pre-computed coefficients from lecture notes (1.26-7)
    std::vector<double> _alpha;
    std::vector<double> _beta;
    std::vector<double> _gamma;
    std::vector<double> _delta;

    void computeSplineCoefficients();
    void solveThomasAlgorithm();
    size_t findInterval(double x) const;
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
    double extrapolate(double x, const InterpolationScheme& interp) const override;

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
