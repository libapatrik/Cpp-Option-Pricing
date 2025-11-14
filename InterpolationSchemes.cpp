//
// Created by Patrik  Liba on 08/11/2025.
//

#include "InterpolationSchemes.h"
#include "Utils.h"
#include<algorithm>
#include<stdexcept>
#include<cmath>


// ============================================================================
// LinearInterpolation Implementation
// ============================================================================

LinearInterpolation::LinearInterpolation(const std::vector<double>& xData, const std::vector<double>& yData)
    : _xData(xData), _yData(yData)
{   /// TODO: validate input data as function
    if (_xData.size() != _yData.size() || _xData.size() < 2) {
        throw std::invalid_argument("LinearInterpolation: Invalid data size");
    }

    // Check if xData is sorted
    if (!std::is_sorted(_xData.begin(), _xData.end())) {
        throw std::invalid_argument("LinearInterpolation: xData must be sorted");
    }

    _xMin = _xData.front();
    _xMax = _xData.back();
}

double LinearInterpolation::interpolate(double x) const
{
    // This method now assumes x is within range [_xMin, _xMax]
    // Bounds checking is done at the higher level (operator())

    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        return _yData.back();
    }

    // Linear interpolation: y = y0 + (y1-y0) * (x-x0)/(x1-x0)
    double x0 = _xData[idx];
    double x1 = _xData[idx + 1];
    double y0 = _yData[idx];
    double y1 = _yData[idx + 1];

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double LinearInterpolation::derivative(double x) const
{
    // For linear interpolation, the derivative is the slope of the interval
    // Analytical form: y = y0 + m(x - x0), so dy/dx = m

    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        idx = _xData.size() - 2; // Use last interval
    }

    // Slope m = (y1 - y0) / (x1 - x0)
    double slope = (_yData[idx + 1] - _yData[idx]) / (_xData[idx + 1] - _xData[idx]);

    return slope;
}

double LinearInterpolation::secondDerivative(double x) const
{
    // For linear interpolation, the second derivative is always 0
    // Since y = mx + b, dy/dx = m (constant), d²y/dx² = 0
    return 0.0;
}

InterpolationSchemes* LinearInterpolation::clone() const
{ 
    return new LinearInterpolation(*this);
}

std::pair<double, double> LinearInterpolation::getRange() const
{
    return {_xMin, _xMax};
}

size_t LinearInterpolation::findInterval(double x) const
{
    // Binary search O(log n) instead of linear search O(n)
    auto it = std::upper_bound(_xData.begin(), _xData.end(), x);

    if (it == _xData.begin()) {
        return 0;
    }

    return std::distance(_xData.begin(), it) - 1;
}

// ============================================================================
// CubicSplineInterpolation Implementation
// ============================================================================

CubicSplineInterpolation::CubicSplineInterpolation(const std::vector<double>& xData, const std::vector<double>& yData, BoundaryType boundaryType)
    : _xData(xData), _yData(yData), _boundaryType(boundaryType)
{   /// TODO: validate input data as function
    if (_xData.size() != _yData.size() || _xData.size() < 2) {
        throw std::invalid_argument("CubicSplineInterpolation: Invalid data size");
    }
    
    // Check if xData is sorted
    if (!std::is_sorted(_xData.begin(), _xData.end())) {
        throw std::invalid_argument("CubicSplineInterpolation: xData must be sorted");
    }
    
    _xMin = _xData.front();
    _xMax = _xData.back();
    
    computeSplineCoefficients();
}

void CubicSplineInterpolation::computeSplineCoefficients()
{ 
    size_t n = _xData.size();

    // Initialize coefficient vectors
    _alpha.resize(n-1);
    _beta.resize(n-1);
    _gamma.resize(n-1);
    _delta.resize(n-1);

    if (n == 2) {
        // Linear case - all coefficients except gamma and delta are zero
        _alpha[0] = 0.0;
        _beta[0] = 0.0;
        _gamma[0] = (_yData[1] - _yData[0]) / (_xData[1] - _xData[0]);
        _delta[0] = _yData[0];
        return;
    }

    solveThomasAlgorithm();
}

void CubicSplineInterpolation::solveThomasAlgorithm()
{
    size_t n = _xData.size();
    
    if (n < 3) {
        return; // Already handled in computeSplineCoefficients
    }
    
    // LECTURE NOTES IMPLEMENTATION: Following equations (1.20) and (1.21)
    // For natural splines, we need to solve for β_2, β_3, ..., β_{N-1}
    // β_1 = 0 and β_N = 0 (natural boundary conditions)
    
    size_t num_unknowns = n - 2;  // β_2, β_3, ..., β_{N-1}
    std::vector<double> lower(num_unknowns, 0.0);
    std::vector<double> diag(num_unknowns, 0.0);
    std::vector<double> upper(num_unknowns, 0.0);
    std::vector<double> rhs(num_unknowns, 0.0);
    
    // Build the tridiagonal system for β coefficients
    // Equation (1.20): β_{j+2}Δx_{j+1} + 2β_{j+1}(Δx_{j+1} + Δx_j) + β_j Δx_j = 3(...)
    for (size_t i = 0; i < num_unknowns; ++i) {
        size_t j = i + 1;  // j = 1, 2, ..., N-2
        
        double dx_j = _xData[j+1] - _xData[j];        // Δx_j
        double dx_j_prev = _xData[j] - _xData[j-1];   // Δx_{j-1}
        
        // Matrix coefficients for β_{j}, β_{j+1}, β_{j+2} following equation (1.20)
        lower[i] = dx_j_prev;                // coefficient of β_j
        diag[i] = 2.0 * (dx_j_prev + dx_j);  // coefficient of β_{j+1}
        upper[i] = dx_j;                     // coefficient of β_{j+2}
        
        // Right-hand side from equation (1.20)
        double slope_j = (_yData[j+1] - _yData[j]) / dx_j;
        double slope_j_prev = (_yData[j] - _yData[j-1]) / dx_j_prev;
        rhs[i] = 3.0 * (slope_j - slope_j_prev);
    }
    
    // Natural boundary conditions are already enforced:
    // - β[0] = 0 (not included in the system, already enforced)
    // - β[n-1] = 0 (not included in the system, will be kept at 0)
    
    // Solve using general TridiagonalSolver (replaces inline Thomas algorithm)
    std::vector<double> beta_interior = ThomasAlgorithm::solve(lower, diag, upper, rhs);
    
    // Construct full β vector with boundary conditions
    // Natural splines: β[0] = 0 and β[n-1] = 0 (boundary conditions)
    std::vector<double> beta(n, 0.0);  // Initialize all to 0
    
    // Fill in interior β values from solver result
    for (size_t i = 0; i < num_unknowns; ++i) {
        beta[i+1] = beta_interior[i];  // β[1], β[2], ..., β[n-2]
    }
    // Note: beta[0] = 0 and beta[n-1] = 0 remain unchanged (natural boundary conditions)
    
    // LECTURE NOTES: Compute all coefficients using equations (1.26)-(1.27)
    for (size_t j = 0; j < n-1; ++j) {
        double dx_j = _xData[j+1] - _xData[j];
        
        // Equation (1.26): α_j = (β_{j+1} - β_j) / (3Δx_j)
        _alpha[j] = (beta[j+1] - beta[j]) / (3.0 * dx_j);
        
        // Equation (1.27): γ_j = (y_{j+1} - y_j) / Δx_j - α_j * Δx_j^2 - β_j * Δx_j
        _gamma[j] = (_yData[j+1] - _yData[j]) / dx_j - _alpha[j] * dx_j * dx_j - beta[j] * dx_j;
        
        // δ_j = y_j (from equation 1.7)
        _delta[j] = _yData[j];
        
        // β_j coefficients (second derivatives at left endpoint of interval j)
        _beta[j] = beta[j];
    }
}

double CubicSplineInterpolation::interpolate(double x) const
{
    // This method now assumes x is within range [_xMin, _xMax]
    // Bounds checking is done at the higher level (operator())
    
    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        return _yData.back();
    }
    
    // LECTURE NOTES IMPLEMENTATION: Following equation (1.6)
    // S_j(x) = α_j(x-x_j)^3 + β_j(x-x_j)^2 + γ_j(x-x_j) + δ_j
    double x0 = _xData[idx];
    double dx = x - x0;    // (x - x_j)
    double dx2 = dx * dx;  // (x - x_j)^2
    double dx3 = dx2 * dx; // (x - x_j)^3
    
    // Use pre-computed coefficients
    double alpha_j = _alpha[idx];
    double beta_j = _beta[idx];
    double gamma_j = _gamma[idx];
    double delta_j = _delta[idx];
    
    // LECTURE NOTES: Apply equation (1.6) exactly
    // S_j(x) = α_j(x-x_j)^3 + β_j(x-x_j)^2 + γ_j(x-x_j) + δ_j
    double S_interp = alpha_j * dx3 + beta_j * dx2 + gamma_j * dx + delta_j;
    
    return S_interp;
}

double CubicSplineInterpolation::derivative(double x) const
{
    // ANALYTICAL DERIVATIVE: S'(x) = γ_j + 2β_j(x-x_j) + 3α_j(x-x_j)²
    // This works for both interpolation and extrapolation at boundaries
    
    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        // At the right boundary, use the derivative from the last interval
        idx = _xData.size() - 2;
    }
    
    double x0 = _xData[idx];
    double dx = x - x0;  // (x - x_j)
    double dx2 = dx * dx;
    
    // Use pre-computed coefficients
    double alpha_j = _alpha[idx];
    double beta_j = _beta[idx];
    double gamma_j = _gamma[idx];
    
    // LECTURE NOTES: S'(x) = γ_j + 2β_j(x-x_j) + 3α_j(x-x_j)²
    double dS_dx = gamma_j + 2.0 * beta_j * dx + 3.0 * alpha_j * dx2;
    
    return dS_dx;
}

double CubicSplineInterpolation::secondDerivative(double x) const
{
    // ANALYTICAL SECOND DERIVATIVE: S''(x) = 2β_j + 6α_j(x-x_j)
    
    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        // At the right boundary, use the second derivative from the last interval
        idx = _xData.size() - 2;
    }
    
    double x0 = _xData[idx];
    double dx = x - x0;  // (x - x_j)
    
    // Use pre-computed coefficients
    double alpha_j = _alpha[idx];
    double beta_j = _beta[idx];
    
    // LECTURE NOTES: S''(x) = 2β_j + 6α_j(x-x_j)
    return 2.0 * beta_j + 6.0 * alpha_j * dx;
}

InterpolationSchemes* CubicSplineInterpolation::clone() const
{
    return new CubicSplineInterpolation(*this);
}

std::pair<double, double> CubicSplineInterpolation::getRange() const
{
    return {_xMin, _xMax};
}

size_t CubicSplineInterpolation::findInterval(double x) const
{
    // upper_bound returns the first element greater than x
    // For x in [xData[i], xData[i+1]), we want to return i
    auto it = std::upper_bound(_xData.begin(), _xData.end(), x);
    
    // Handle edge case: if x <= xData[0], return interval 0
    if (it == _xData.begin()) {
        return 0;
    }
    
    return std::distance(_xData.begin(), it) - 1;
}

double QuadraticExtrapolationScheme::extrapolate(double x) const
{
    if (x <= _xMin)
    {
        double delta = x - _xMin;
        return _xMin + _dxMin * delta + 0.5 _d2xMin * delta * delta;
    }
    else if (x >= _xMax)
    {
        double delta = x - _xMax;
        return _xMax + _dxMax * delta + 0.5 _d2xMax * delta * delta;
    }
    else
        throw "x is not correct";
}
