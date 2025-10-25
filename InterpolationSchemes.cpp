//
// Created by Patrik  Liba on 05/10/2025.
//

#include "InterpolationSchemes.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>

// ============================================================================
// LinearInterpolation Implementation
// ============================================================================

LinearInterpolation::LinearInterpolation(const std::vector<double>& xData, const std::vector<double>& yData)
    : _xData(xData), _yData(yData)
{
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
    if (x < _xMin || x > _xMax) {
        return extrapolate(x);
    }
    
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
    return std::distance(_xData.begin(), it) - 1;
}

double LinearInterpolation::extrapolate(double x) const
{
    if (x < _xMin) {
        // Left extrapolation: linear continuation from first two points
        double slope = (_yData[1] - _yData[0]) / (_xData[1] - _xData[0]);
        return _yData[0] + slope * (x - _xData[0]);
    } else {
        // Right extrapolation: linear continuation from last two points
        size_t n = _xData.size();
        double slope = (_yData[n-1] - _yData[n-2]) / (_xData[n-1] - _xData[n-2]);
        return _yData[n-1] + slope * (x - _xData[n-1]);
    }
}

// ============================================================================
// CubicSplineInterpolation Implementation
// ============================================================================

CubicSplineInterpolation::CubicSplineInterpolation(const std::vector<double>& xData, const std::vector<double>& yData, BoundaryType boundaryType)
    : _xData(xData), _yData(yData), _boundaryType(boundaryType)
{
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
        // Linear case - all coefficients except delta are zero
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
    std::vector<double> a(num_unknowns), b(num_unknowns), c(num_unknowns), d(num_unknowns);
    
    // Build the tridiagonal system for β coefficients
    // Equation (1.20): β_{j+2}Δx_{j+1} + 2β_{j+1}(Δx_{j+1} + Δx_j) + β_j Δx_j = 3(...)
    for (size_t i = 0; i < num_unknowns; ++i) {
        size_t j = i + 1;  // j = 1, 2, ..., N-2
        
        double dx_j = _xData[j+1] - _xData[j];        // Δx_j
        double dx_j_prev = _xData[j] - _xData[j-1];   // Δx_{j-1}
        
        // Matrix coefficients for β_{j}, β_{j+1}, β_{j+2} following equation (1.20)
        if (i > 0) {
            a[i] = dx_j_prev;  // coefficient of β_j
        }
        b[i] = 2.0 * (dx_j_prev + dx_j);  // coefficient of β_{j+1}
        if (i < num_unknowns - 1) {
            c[i] = dx_j;  // coefficient of β_{j+2}
        }
        
        // Right-hand side from equation (1.20)
        double slope_j = (_yData[j+1] - _yData[j]) / dx_j;
        double slope_j_prev = (_yData[j] - _yData[j-1]) / dx_j_prev;
        d[i] = 3.0 * (slope_j - slope_j_prev);
    }
    
    // Natural boundary conditions are already enforced:
    // - β[0] = 0 (not included in the system, already enforced)
    // - β[n-1] = 0 (not included in the system, will be kept at 0)
    // The tridiagonal system equations from the loop above are correct for natural splines
    
    // LECTURE NOTES THOMAS ALGORITHM: Following equations (1.22)-(1.24)
    std::vector<double> c_prime(num_unknowns), r_prime(num_unknowns);
    
    // Forward elimination - Equations (1.22)-(1.23)
    c_prime[0] = c[0] / b[0];
    r_prime[0] = d[0] / b[0];
    
    for (size_t i = 1; i < num_unknowns; ++i) {
        double denom = b[i] - a[i] * c_prime[i-1];
        c_prime[i] = c[i] / denom;
        r_prime[i] = (d[i] - a[i] * r_prime[i-1]) / denom;
    }
    
    // Back substitution - Equation (1.24)
    // Natural splines: β[0] = 0 and β[n-1] = 0 (boundary conditions)
    // Solve for β[1], β[2], ..., β[n-2] only
    std::vector<double> beta(n, 0.0);  // Initialize all to 0
    
    // Back-substitute to find interior β values
    if (num_unknowns > 0) {
        // Start from the last unknown and work backwards
        beta[n-2] = r_prime[num_unknowns-1];  // β[n-2] (last interior point)
        
        for (int i = static_cast<int>(num_unknowns) - 2; i >= 0; --i) {
            beta[i+1] = r_prime[i] - c_prime[i] * beta[i+2];
        }
    }
    // Note: beta[0] = 0 and beta[n-1] = 0 remain unchanged (natural boundary conditions)
    
    // LECTURE NOTES: Compute all coefficients using equations (1.26)-(1.27)
    for (size_t j = 0; j < n-1; ++j) {
        double dx_j = _xData[j+1] - _xData[j];
        
        // Equation (1.26): α_j = (β_{j+1} - β_j) / (3Δx_j)
        // For interval j, we use second derivatives at nodes j and j+1
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
    if (x < _xMin || x > _xMax) { // check if out of bounds
        return extrapolate(x);
    }
    
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
    if (x < _xMin || x > _xMax) { // check if out of bounds
        throw std::out_of_range("CubicSplineInterpolation: x is outside interpolation range");
    }
    
    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        // At the right boundary, use the derivative from the last interval
        idx = _xData.size() - 2;
    }
    
    // LECTURE NOTES IMPLEMENTATION: Using pre-computed coefficients
    // S'(x) = γ_j + 2β_j(x-x_j) + 3α_j(x-x_j)²
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
    if (x < _xMin || x > _xMax) {
        throw std::out_of_range("CubicSplineInterpolation: x is outside interpolation range");
    }
    
    size_t idx = findInterval(x);
    if (idx >= _xData.size() - 1) {
        // At the right boundary, use the second derivative from the last interval
        idx = _xData.size() - 2;
    }
    
    // LECTURE NOTES IMPLEMENTATION: Using pre-computed coefficients
    // S''(x) = 2β_j + 6α_j(x-x_j)
    double x0 = _xData[idx];
    double dx = x - x0;  // (x - x_j)
    
    // Use pre-computed coefficients
    double alpha_j = _alpha[idx];
    double beta_j = _beta[idx];
    
    // LECTURE NOTES: S''(x) = 2β_j + 6α_j(x-x_j)
    double d2S_dx2 = 2.0 * beta_j + 6.0 * alpha_j * dx;
    
    return d2S_dx2;
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

double CubicSplineInterpolation::extrapolate(double x) const
{
    // LECTURE NOTES IMPLEMENTATION: Following equations (1.32), (1.26), and (1.27)
    // σ*(T_i, K) = σ*(T_i, K_1) + D_L × (K - K_1) if K < K_1
    // σ*(T_i, K) = σ*(T_i, K_N) + D_R × (K - K_N) if K > K_N
    
    if (x < _xMin) {
        // LECTURE NOTES: Left extrapolation using equations (1.26) and (1.27)
        // D_L = S'_1(x_1) = γ_1
        // Use pre-computed γ_0 coefficient
        double gamma0 = _gamma[0];
        
        return _yData[0] + gamma0 * (x - _xData[0]);
    } else {
        // LECTURE NOTES: Right extrapolation using equations (1.26) and (1.27)
        // D_R = S'_{N-1}(x_N) = 3α_{N-1}Δx_{N-1}^2 + 2β_{N-1}Δx_{N-1} + γ_{N-1}
        size_t n = _xData.size();
        size_t last_idx = n - 2;  // Index of last interval
        
        // Use pre-computed coefficients
        double alpha_n = _alpha[last_idx];
        double beta_n = _beta[last_idx];
        double gamma_n = _gamma[last_idx];
        
        double hn = _xData[n-1] - _xData[n-2];
        
        // D_R = 3α_{N-1}Δx_{N-1}^2 + 2β_{N-1}Δx_{N-1} + γ_{N-1}
        double D_R = 3.0 * alpha_n * hn * hn + 2.0 * beta_n * hn + gamma_n;
        
        return _yData[n-1] + D_R * (x - _xData[n-1]);
    }
}

// ============================================================================