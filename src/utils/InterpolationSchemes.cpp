#include <cppfm/utils/InterpolationSchemes.h>
#include <cppfm/utils/Utils.h>
#include<algorithm>
#include<stdexcept>
#include<cmath>


// ============================================================================
// InterpolationScheme Base Class Implementation
// ============================================================================
InterpolationScheme::InterpolationScheme(const std::vector<double>& xData, const std::vector<double>& yData, ExtrapolationType extraType)
    : _xData(xData), _yData(yData), _extrapolationType(extraType)
{
    // Validate input data
    validateData();
    
    // Create extrapolation scheme based on enum (but don't initialize yet)
    // Initialization happens in derived class constructors after their setup is complete
    switch (extraType) {
        case ExtrapolationType::Flat:
            _extrapolationScheme = std::make_unique<FlatExtrapolation>();
            break;
        case ExtrapolationType::Linear:
            _extrapolationScheme = std::make_unique<LinearExtrapolation>();
            break;
        case ExtrapolationType::Quadratic:
            _extrapolationScheme = std::make_unique<QuadraticExtrapolation>();
            break;
    }
}

void InterpolationScheme::validateData() const
{
    if (_xData.size() != _yData.size()) {
        throw std::invalid_argument("InterpolationScheme: xData and yData must have same size");
    }
    
    if (_xData.size() < 2) {
        throw std::invalid_argument("InterpolationScheme: At least 2 data points required");
    }

    if (!std::is_sorted(_xData.begin(), _xData.end())) {
        throw std::invalid_argument("InterpolationScheme: xData must be sorted in ascending order");
    }
}

std::pair<double, double> InterpolationScheme::getRange() const
{
    return {_xData.front(), _xData.back()};
}

size_t InterpolationScheme::findInterval(double x) const
{
    // Returns index i such that x is in [xData[i], xData[i+1])
    // Uses binary search: O(log n) complexity
    auto it = std::upper_bound(_xData.begin(), _xData.end(), x);
    
    if (it == _xData.begin()) {
        return 0;
    }

    size_t idx = std::distance(_xData.begin(), it) -1;
    // Clamp to valid range
    if (idx >= _xData.size() - 1) {
        idx = _xData.size() - 2;
    }
    return idx;
}

double InterpolationScheme::operator()(double x) const
{   // Main interface to interpolate or extrapolate automatically
    auto [xMin, xMax] = getRange();
        // return _extrapolationScheme->extrapolate(x);
    // Route to appropriate method based on x location
    if (x < xMin || x > xMax) {
        // Always use extrapolation strategy
        // If no custom scheme is set, default QuadraticExtrapolation is used
        return _extrapolationScheme->extrapolate(x, *this);   // *this passes the interpolation object
        // Why? Extrapolation needs to query ('ask for information') the interpolation object to get the boundary values and derivatives
    } else {
        return interpolate(x);
    }
}

// ============================================================================
// LinearInterpolation Implementation
// ============================================================================

LinearInterpolation::LinearInterpolation(const std::vector<double>& xData, 
                                         const std::vector<double>& yData,
                                         ExtrapolationType extraType)
    : InterpolationScheme(xData, yData, extraType)
{
    // Initialize extrapolation scheme now that setup is complete
    _extrapolationScheme->initialize(*this);
}

double LinearInterpolation::interpolate(double x) const
{
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

std::unique_ptr<InterpolationScheme> LinearInterpolation::clone() const
{
    auto cloned = std::make_unique<LinearInterpolation>(_xData, _yData, _extrapolationType);
    // _extrapolationScheme is already initialized in the constructor
    // But we need to clone it to preserve the state
    cloned->_extrapolationScheme = _extrapolationScheme->clone();
    cloned->_extrapolationScheme->initialize(*cloned);
    return cloned;
}



// ============================================================================
// CubicSplineInterpolation Implementation
// ============================================================================

CubicSplineInterpolation::CubicSplineInterpolation(const std::vector<double>& xData, 
                                                   const std::vector<double>& yData, 
                                                   BoundaryType boundaryType,
                                                   ExtrapolationType extraType)
    : InterpolationScheme(xData, yData, extraType)
    , _boundaryType(boundaryType)
{
    computeSplineCoefficients();
    
    // Initialize extrapolation scheme now that spline setup is complete
    _extrapolationScheme->initialize(*this);
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

    double dx = x - _xData[idx];  // (x - x_j)
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

    double dx = x - _xData[idx];  // (x - x_j)
    
    // Use pre-computed coefficients
    double alpha_j = _alpha[idx];
    double beta_j = _beta[idx];

    // LECTURE NOTES: S''(x) = 2β_j + 6α_j(x-x_j)
    return 2.0 * beta_j + 6.0 * alpha_j * dx;
}

std::unique_ptr<InterpolationScheme> CubicSplineInterpolation::clone() const
{
    auto cloned = std::make_unique<CubicSplineInterpolation>(_xData, _yData, _boundaryType, _extrapolationType);
    // _extrapolationScheme is already initialized in the constructor
    // But we need to clone it to preserve the state
    cloned->_extrapolationScheme = _extrapolationScheme->clone();
    cloned->_extrapolationScheme->initialize(*cloned);
    return cloned;
}


// ============================================================================
// ExtrapolationScheme Implementations
// ============================================================================

// FlatExtrapolation
std::unique_ptr<ExtrapolationScheme> FlatExtrapolation::clone() const
{
    auto cloned = std::make_unique<FlatExtrapolation>();
    cloned->_yMin = _yMin;
    cloned->_yMax = _yMax;
    return cloned;
}

void FlatExtrapolation::initialize(const InterpolationScheme& interp)
{
    auto [xMin, xMax] = interp.getRange();
    _yMin = interp.interpolate(xMin);
    _yMax = interp.interpolate(xMax);
}

double FlatExtrapolation::extrapolate(double x, const InterpolationScheme& interp) const
{
    auto [xMin, xMax] = interp.getRange();
    return (x < xMin) ? _yMin : _yMax;        // use the boundary values 
}

// LinearExtrapolation
std::unique_ptr<ExtrapolationScheme> LinearExtrapolation::clone() const
{
    auto cloned = std::make_unique<LinearExtrapolation>();
    cloned->_yMin = _yMin;
    cloned->_yMax = _yMax;
    cloned->_dyMin = _dyMin;
    cloned->_dyMax = _dyMax;
    return cloned;
}

void LinearExtrapolation::initialize(const InterpolationScheme& interp)
{
    auto [xMin, xMax] = interp.getRange();
    _yMin = interp.interpolate(xMin);
    _yMax = interp.interpolate(xMax);
    _dyMin = interp.derivative(xMin);
    _dyMax = interp.derivative(xMax);
}

double LinearExtrapolation::extrapolate(double x, const InterpolationScheme& interp) const
{
    auto [xMin, xMax] = interp.getRange();
    
    if (x < xMin) {
        return _yMin + _dyMin * (x - xMin);
    } else {
        return _yMax + _dyMax * (x - xMax);
    }
}

// QuadraticExtrapolation
std::unique_ptr<ExtrapolationScheme> QuadraticExtrapolation::clone() const
{   // Deep copy of member variables
    // 1. Create a new QuadraticExtrapolation object
    auto cloned = std::make_unique<QuadraticExtrapolation>();
    // 2. Copy all the member variables from 'this' to the cloned object
    cloned->_yMin = _yMin;         // copy the boundary values
    cloned->_yMax = _yMax;         // copy the boundary values
    cloned->_dyMin = _dyMin;       // copy the derivative at the boundary
    cloned->_dyMax = _dyMax;       // copy the derivative at the boundary
    cloned->_d2yMin = _d2yMin;     // copy the second derivative at the boundary
    cloned->_d2yMax = _d2yMax;     // copy the second derivative at the boundary
    
    return cloned;                 // return the cloned object
}
// Lazy/eager evaluation and cache
void QuadraticExtrapolation::initialize(const InterpolationScheme& interp)
{
    // Compute expensives once during initialization
    // and cache (store values - so don't need to compute again) the results for later use
    auto [xMin, xMax] = interp.getRange();
    _yMin = interp.interpolate(xMin);          // cache the boundary values
    _yMax = interp.interpolate(xMax);
    _dyMin = interp.derivative(xMin);          // cache the derivative at the boundary
    _dyMax = interp.derivative(xMax);
    _d2yMin = interp.secondDerivative(xMin);   // cache the second derivative at the boundary
    _d2yMax = interp.secondDerivative(xMax);
}

double QuadraticExtrapolation::extrapolate(double x, const InterpolationScheme& interp) const
{
    auto [xMin, xMax] = interp.getRange();
    double result;

    if (x < xMin) {
        double dx = x - xMin;
        result = _yMin + _dyMin * dx + 0.5 * _d2yMin * dx * dx;
    } else {
        double dx = x - xMax;
        result = _yMax + _dyMax * dx + 0.5 * _d2yMax * dx * dx;
    }

    /// TODO: Add fallbacks for min and max vols?
    const double MIN_VOL = 0.01;   // 1%
    const double MAX_VOL = 2.0;    // 200%

    return std::clamp(result, MIN_VOL, MAX_VOL);   // return the clamped value
}