#include "Grid.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================
Grid::Grid(double S_min, double S_max, size_t N_S,
                 double T_max, size_t N_t,
                 SpacingType spacingType)
    : _S_min(S_min), _S_max(S_max), _N_S(N_S),
      _T_max(T_max), _N_t(N_t), _spacingType(spacingType)
{
    // Validate parameters
    validateParameters();
    
    // Compute time step (constant for all time steps)
    _dt = _T_max / static_cast<double>(_N_t);
    
    // Construct temporal grid (always uniform)
    constructTimeGrid();
    
    // Construct spatial grid (uniform or log-spaced)
    if (_spacingType == SpacingType::Uniform) {
        constructUniformGrid();
    } else {
        constructLogSpacedGrid();
    }
}

// Clone
Grid* Grid::clone() const
{
    return new Grid(*this);
}

// Utility Methods
std::pair<std::pair<double, double>, std::pair<double, double>>Grid::getBounds() const
{
    return {{_S_min, _S_max}, {0.0, _T_max}};
}

// ============================================================================
// Grid Construction Methods
// ============================================================================

void Grid::validateParameters() const
{
    // Validate spatial domain
    if (_S_min < 0.0) {
        throw std::invalid_argument("PDEGrid: S_min must be non-negative");
    }
    if (_S_max <= _S_min) {
        throw std::invalid_argument("PDEGrid: S_max must be greater than S_min");
    }
    
    // Validate grid dimensions
    if (_N_S < 3) {
        throw std::invalid_argument("PDEGrid: Need at least 3 spatial points for finite differences");
    }
    if (_N_t < 1) {
        throw std::invalid_argument("PDEGrid: Need at least 1 time step");
    }
    
    // Validate temporal domain
    if (_T_max <= 0.0) {
        throw std::invalid_argument("PDEGrid: T_max must be positive");
    }
}

void Grid::constructTimeGrid()
{
    // Reserve space for efficiency
    _timeGrid.reserve(_N_t + 1);
    
    // Create uniform time grid: t_j = j * Δt for j = 0, 1, ..., N_t
    for (size_t j = 0; j <= _N_t; ++j) {
        _timeGrid.push_back(static_cast<double>(j) * _dt);
    }
}

void Grid::constructUniformGrid()
{
    // Reserve space for efficiency
    _spotGrid.reserve(_N_S);
    
    // Compute uniform spacing
    double dS = (_S_max - _S_min) / static_cast<double>(_N_S - 1);
    
    // Create uniform spatial grid: S_i = S_min + i * ΔS for i = 0, 1, ..., N_S-1
    for (size_t i = 0; i < _N_S; ++i) {
        _spotGrid.push_back(_S_min + static_cast<double>(i) * dS);
    }
}

void Grid::constructLogSpacedGrid()
{
    // Reserve space for efficiency
    _spotGrid.reserve(_N_S);
    
    // Log-spaced grid: S_i = exp(ξ_i) where ξ_i is uniformly spaced in log-space
    // This concentrates more grid points near S_min (useful for low strikes)
    // Handle S_min = 0 case (add small epsilon to avoid log(0))
    const double epsilon = 1e-10;
    double S_min_safe = std::max(_S_min, epsilon);
    
    if (_S_min < 0.0) {
        throw std::invalid_argument("Grid: S_min must be non-negative");
    }
    
    double log_S_min = std::log(S_min_safe);
    double log_S_max = std::log(_S_max);
    double d_log_S = (log_S_max - log_S_min) / static_cast<double>(_N_S - 1);
    
    // Create log-spaced grid: S_i = exp(log_S_min + i * Δlog_S)
    for (size_t i = 0; i < _N_S; ++i) {
        double log_S_i = log_S_min + static_cast<double>(i) * d_log_S;
        _spotGrid.push_back(std::exp(log_S_i));
    }
    
    // Ensure exact boundary values (avoid numerical errors)
    _spotGrid[0] = S_min_safe;
    _spotGrid[_N_S - 1] = _S_max;
}

// ============================================================================
// Grid Spacing Methods
// ============================================================================

double Grid::dS(size_t i) const
{
    // For finite difference schemes, we need the effective grid spacing at point i
    // Use centered difference where possible for better accuracy
    
    if (i == 0) {
        // At lower boundary, use forward difference
        return _spotGrid[1] - _spotGrid[0];
    } else if (i == _N_S - 1) {
        // At upper boundary, use backward difference
        return _spotGrid[_N_S - 1] - _spotGrid[_N_S - 2];
    } else {
        // Interior points: use centered difference
        // dS(i) = (S_{i+1} - S_{i-1}) / 2
        return (_spotGrid[i + 1] - _spotGrid[i - 1]) / 2.0;
    }
}

double Grid::dS_forward(size_t i) const
{
    if (i >= _N_S - 1) {
        throw std::out_of_range("PDEGrid::dS_forward: index out of range");
    }
    return _spotGrid[i + 1] - _spotGrid[i];
}

double Grid::dS_backward(size_t i) const
{
    if (i == 0) {
        throw std::out_of_range("PDEGrid::dS_backward: index out of range");
    }
    return _spotGrid[i] - _spotGrid[i - 1];
}


