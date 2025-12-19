#include "Utils.h"
#include <cmath>
#include <stdexcept>
#include <complex>
#include <numbers>
#include <boost/math/distributions/non_central_chi_squared.hpp>  // for sampling from non-central chi-squared distribution for CIRSampler
#include <boost/math/special_functions/bessel.hpp>               // for modified Bessel functions I_ν(z)
#include <boost/random/uniform_01.hpp>

// ============================================================================
// * Statistics
// ============================================================================

constexpr double PI = std::numbers::pi; // C++ 20 pi

// Standard normal CDF
double Utils::stdNormCdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

double Utils::stdNormPdf(double x)
{
    // PDF of standard normal distribution: φ(x) = (1/√(2π)) * exp(-x²/2)
    return (1.0 / std::sqrt(2.0 * PI)) * std::exp(-0.5 * x * x);
}

std::complex<double> Utils::stdNormChF(double omega)
{
    // Characteristic function of N(0,1): φ(ω) = exp(-ω²/2)
    return std::exp(-0.5 * omega * omega);
}

double Utils::inverseNormalCDF(double u)
{
    /**
     * Moro's Inverse Normal CDF Algorithm (1995)
     * 
     * Approximates Φ⁻¹(u) where Φ is the standard normal CDF
     * 
     * Accuracy:
     * - Central region (0.08 < u < 0.92): |error| < 3×10⁻⁹
     * - Tail regions (u ≤ 0.08 or u ≥ 0.92): |error| < 10⁻⁷
     * 
     * Reference: Moro, B. (1995), "The Full Monte", RISK, Vol. 8, No. 2
     */
    
    // Clamp to safe range to prevent log(0) or log(negative)
    const double eps = 1e-15;
    u = std::max(eps, std::min(u, 1.0 - eps));
    
    // Transform to standard normal variable
    double x = u - 0.5;
    
    if (std::abs(x) < 0.42) {
        // ====================================================================
        // Central region: |u - 0.5| < 0.42  (i.e., 0.08 < u < 0.92)
        // ====================================================================
        double r = x * x;
        
        // Numerator: a₀ + a₁r + a₂r² + a₃r³
        static const double a0 =  2.50662823884;
        static const double a1 = -18.61500062529;
        static const double a2 =  41.39119773534;
        static const double a3 = -25.44106049637;
        
        // Denominator: 1 + b₀r + b₁r² + b₂r³ + b₃r⁴
        static const double b0 = -8.47351093090;
        static const double b1 =  23.08336743743;
        static const double b2 = -21.06224101826;
        static const double b3 =   3.13082909833;
        
        double num = a0 + r * (a1 + r * (a2 + r * a3));
        double den = 1.0 + r * (b0 + r * (b1 + r * (b2 + r * b3)));
        
        return x * num / den;
        
    } else {
        // ====================================================================
        // Tail regions: u ≤ 0.08 or u ≥ 0.92
        // ====================================================================
        double r = (x < 0.0) ? u : (1.0 - u);
        r = std::log(-std::log(r));
        
        // Polynomial: c₀ + c₁r + c₂r² + ... + c₈r⁸
        static const double c0 = 0.3374754822726147;
        static const double c1 = 0.9761690190917186;
        static const double c2 = 0.1607979714918209;
        static const double c3 = 0.0276438810333863;
        static const double c4 = 0.0038405729373609;
        static const double c5 = 0.0003951896511919;
        static const double c6 = 0.0000321767881768;
        static const double c7 = 0.0000002888167364;
        static const double c8 = 0.0000003960315187;
        
        double z = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + 
                   r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
        
        return (x < 0.0) ? -z : z;
    }
}


// ============================================================================
// * Tridiagonal System Solver Implementation
// ============================================================================

std::vector<double> ThomasAlgorithm::solve(const std::vector<double>& lower,
                                           const std::vector<double>& diag,
                                           const std::vector<double>& upper,
                                           const std::vector<double>& rhs)
{
    // Validate input sizes
    size_t n = rhs.size();
    
    if (lower.size() != n || diag.size() != n || upper.size() != n) {
        throw std::invalid_argument("TridiagonalSolver::solve: All input vectors must have the same size");
    }
    
    if (n == 0) {
        throw std::invalid_argument("TridiagonalSolver::solve: Cannot solve empty system");
    }
    
    if (n == 1) {
        // Trivial case: b[0]*x[0] = d[0]
        if (std::abs(diag[0]) < 1e-14) {
            throw std::runtime_error(
                "TridiagonalSolver::solve: Singular matrix (diagonal element is zero)");
        }
        return {rhs[0] / diag[0]};
    }

    // Diagonal dominance, M-matrices etc.
    
    // ========================================================================
    // THOMAS ALGORITHM (General Implementation)
    // ========================================================================
    
    // Temporary storage for modified coefficients
    std::vector<double> c_prime(n, 0.0);  // Modified upper diagonal
    std::vector<double> d_prime(n, 0.0);  // Modified RHS
    std::vector<double> x(n, 0.0);        // Solution vector
    
    // ------------------------------------------------------------------------
    /// PROCESS: Mx = r -> LUx = r -> (1) Ly = r -> (2) Ux = y
    // Step 1: Forward Elimination
    // ------------------------------------------------------------------------
    // Transform the system to upper triangular form
    // Original:  b[i]*x[i] + c[i]*x[i+1] = d[i] - a[i]*x[i-1]
    // Modified:  x[i] + c'[i]*x[i+1] = d'[i]
    
    // First row (i = 0): b[0]*x[0] + c[0]*x[1] = d[0]
    if (std::abs(diag[0]) < 1e-14) {
        throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row 0");
    }
    c_prime[0] = upper[0] / diag[0];         // γ1 = c1/b1
    d_prime[0] = rhs[0] / diag[0];           // ρ1 = r1/b1
    
    // Interior rows (i = 1, 2, ..., n-1)
    for (size_t i = 1; i < n; ++i) {
        // Denominator after eliminating x[i-1]
        double m = diag[i] - lower[i] * c_prime[i-1];
        //          b[i] -  a[i] * γ[i-1]
        if (std::abs(m) < 1e-14) {
            throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row " + std::to_string(i));
        }
        
        c_prime[i] = upper[i] / m;
        // c[i] = c[i+1] / (b[i] - a[i] * γ[i-1])
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i-1]) / m;
        // d'[i] = (d[i] - a[i] * ρ[i-1]) / (b[i] - a[i] * γ[i-1])
    }
    
    // ------------------------------------------------------------------------
    // Step 2: Back Substitution
    // ------------------------------------------------------------------------
    // Solve the upper triangular system from bottom to top
    // x[i] = d'[i] - c'[i]*x[i+1]
    
    // Last equation (i = n-1): x[n-1] = d'[n-1]
    x[n-1] = d_prime[n-1];
    
    // Work backwards (i = n-2, n-3, ..., 0)
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
        //     ρ[i+1]     - γ[i+1]     * x[i+2]
        //  (Code i=k corresponds to PDF row k+1)
    }
    
    return x;
}       // end of ThomasAlgorithm::solve


// ============================================================================
// * Numerical Derivatives Implementation
// ============================================================================

double NumericalDerivatives::firstDerivative(std::function<double(double)> f, double x, double h)
{
    // Scale-invariant step size: h_scaled = h * max(|x|, 1.0)
    // This ensures relative perturbation is consistent across different magnitudes
    double h_scaled = h * std::max(std::abs(x), 1.0);
    
    // Central difference: f'(x) = [f(x+h) - f(x-h)] / (2h)
    return (f(x + h_scaled) - f(x - h_scaled)) / (2.0 * h_scaled);
}

double NumericalDerivatives::secondDerivative(std::function<double(double)> f, double x, double h)
{
    // Scale-invariant step size
    double h_scaled = h * std::max(std::abs(x), 1.0);
    
    // Central second difference: f''(x) = [f(x+h) - 2f(x) + f(x-h)] / h^2
    return (f(x + h_scaled) - 2.0 * f(x) + f(x - h_scaled)) / (h_scaled * h_scaled);
}

// ============================================================================
// * COS Method Implementation
// ============================================================================

std::vector<double> Transforms::recoverCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x)
{
    /**
     * COS Method for CDF Recovery
     * 
     * F(x) ≈ 1/2 * F₀(x-a)/(b-a) + Σ[k=1 to N-1] Fₖ/ωₖ · sin(ωₖ(x-a))
     * 
     * where:
     *   ωₖ = kπ/(b-a)
     *   Fₖ = (2/(b-a)) · Re[φ(ωₖ) · exp(-iωₖa)]
     *   φ(ω) is the characteristic function
     */
    
    std::vector<double> cdf(x.size(), 0.0);
    
    // Compute F_k coefficients
    std::vector<double> F_k(N);
    std::vector<double> omega(N);
    
    for (size_t k = 0; k < N; ++k) {
        omega[k] = k * PI / (b - a);
        
        // F_k = (2/(b-a)) · Re[φ(ω_k) · exp(-iω_k·a)]
        std::complex<double> chf_val = chf(omega[k]);
        std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -omega[k] * a));
        F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
    }
    
    // Compute CDF at each point x
    for (size_t i = 0; i < x.size(); ++i) {
        double x_shifted = x[i] - a;
        
        // First term: F₀/2 · (x-a)/(b-a)
        cdf[i] = F_k[0] / 2.0 * x_shifted;
        
        // Sum over k = 1, 2, ..., N-1
        for (size_t k = 1; k < N; ++k) {
            cdf[i] += F_k[k] / omega[k] * std::sin(omega[k] * x_shifted);
        }
    }
    
    return cdf;
}

std::vector<double> Transforms::recoverPDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x)
{
    /**
     * COS Method for PDF Recovery
     * 
     * f(x) ≈ Σ[k=0 to N-1] Fₖ · cos(ωₖ(x-a))
     * 
     * where:
     *   ωₖ = kπ/(b-a)
     *   Fₖ = (2/(b-a)) · Re[φ(ωₖ) · exp(-iωₖa)]
     *   F₀ is halved (first term correction)
     */
    
    std::vector<double> pdf(x.size(), 0.0);
    
    // Compute F_k coefficients
    std::vector<double> F_k(N);
    std::vector<double> omega(N);
    
    for (size_t k = 0; k < N; ++k) {
        omega[k] = k * PI / (b - a);
        
        // F_k = (2/(b-a)) · Re[φ(ω_k) · exp(-iω_k·a)]
        std::complex<double> chf_val = chf(omega[k]);
        std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -omega[k] * a));
        F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
    }
    
    // First term is halved
    F_k[0] *= 0.5;
    
    // Compute PDF at each point x
    for (size_t i = 0; i < x.size(); ++i) {
        double x_shifted = x[i] - a;
        
        // Sum over k = 0, 1, ..., N-1
        for (size_t k = 0; k < N; ++k) {
            pdf[i] += F_k[k] * std::cos(omega[k] * x_shifted);
        }
    }
    
    return pdf;
}

std::pair<double, size_t> Transforms::invertCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, double p, size_t max_iter, double tol)
{
    /**
     * Newton-Raphson CDF Inversion for COS Method
     * 
     * Find x such that F(x) = p where F is CDF recovered by COS method.
     * 
     * Algorithm:
     * 1. Evaluate CDF at multiple points to get good initial guess
     * 2. Newton iteration: x := x - (F(x) - p) / f(x)
     *    where f(x) = PDF(x) = dF/dx
     * 3. Apply bounds and convergence checks
     * 
     * Reference: Standard numerical root-finding with COS method
     */
    
    // Clamp probability to [0, 1]
    p = std::max(0.0, std::min(1.0, p));
    
    // Handle edge cases
    if (p <= 0.0) return {a, 0};
    if (p >= 1.0) return {b, 0};
    
    // ========================================================================
    // STEP 1: Find good initial guess by evaluating CDF at multiple points
    // ========================================================================
    
    // Number of initial evaluation points (sweet spot: 25-30 based on testing)
    const size_t num_initial_points = 30;
    std::vector<double> x_initial(num_initial_points);
    
    // Create uniform grid over [a, b]
    for (size_t i = 0; i < num_initial_points; ++i) {
        x_initial[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(num_initial_points - 1);
    }
    
    // Evaluate CDF at initial points
    std::vector<double> cdf_initial = recoverCDF(a, b, N, chf, x_initial);
    
    // Find point closest to target probability p
    size_t best_idx = 0;
    double min_dist = std::abs(cdf_initial[0] - p);
    for (size_t i = 1; i < num_initial_points; ++i) {
        double dist = std::abs(cdf_initial[i] - p);
        if (dist < min_dist) {
            min_dist = dist;
            best_idx = i;
        }
    }
    
    // Initial guess for Newton-Raphson
    double x = x_initial[best_idx];
    
    // ========================================================================
    // STEP 2: Newton-Raphson iteration
    // ========================================================================
    
    for (size_t iter = 0; iter < max_iter; ++iter) {
        // Evaluate CDF at current point
        std::vector<double> x_vec = {x};
        std::vector<double> cdf_x_vec = recoverCDF(a, b, N, chf, x_vec);
        double cdf_x = cdf_x_vec[0];
        
        // Compute residual: f(x) = F(x) - p
        double fx = cdf_x - p;
        
        // Check convergence: |F(x) - p| < tol
        if (std::abs(fx) < tol) {
            return {x, iter};                   // Return x and number of iterations
        }
        
        // Evaluate PDF at current point (derivative of CDF)
        std::vector<double> pdf_x_vec = recoverPDF(a, b, N, chf, x_vec);
        double pdf_x = pdf_x_vec[0];
        
        // Safeguard against division by very small PDF
        // If PDF ≈ 0, we're likely at a tail where CDF is flat
        if (std::abs(pdf_x) < 1e-10) {
            pdf_x = (pdf_x >= 0.0) ? 1e-10 : -1e-10;
        }
        
        // Newton step: x_new = x - f(x) / f'(x) = x - (F(x) - p) / PDF(x)
        double x_new = x - fx / pdf_x;
        
        // Keep within bounds [a, b]
        x_new = std::max(a, std::min(b, x_new));
        
        // Check for convergence in x
        if (std::abs(x_new - x) < tol) {
            return {x_new, iter + 1};         // Return x and number of iterations
        }
        
        // Update for next iteration
        x = x_new;
    }
    
    // If max iterations reached, return current best estimate
    // This shouldn't happen often with good initial guess and reasonable parameters
    return {x, max_iter};                        // Return x and number of iterations
}

// ============================================================================
// Optimized Newton-Raphson with Coefficient Caching
// ============================================================================

Transforms::Coefficients Transforms::precomputeCoefficients(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf)
{
    /**
     * Precompute Fourier Coefficients for Fast CDF/PDF Evaluation
     * 
     * This is the expensive step: O(N) evaluations of the characteristic function.
     * Once computed, we can evaluate CDF/PDF at any point x in O(N) arithmetic operations.
     * 
     * F_k = (2/(b-a)) · Re[φ(ω_k) · exp(-iω_k·a)]
     * ω_k = kπ/(b-a)
     */
    
    Coefficients coeffs;
    coeffs.a = a;
    coeffs.b = b;
    coeffs.N = N;
    coeffs.F_k.resize(N);
    coeffs.omega.resize(N);
    
    // Precompute F_k and ω_k for k = 0, 1, ..., N-1
    for (size_t k = 0; k < N; ++k) {
        coeffs.omega[k] = k * PI / (b - a);
        
        // F_k = (2/(b-a)) · Re[φ(ω_k) · exp(-iω_k·a)]
        std::complex<double> chf_val = chf(coeffs.omega[k]);
        std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -coeffs.omega[k] * a));
        coeffs.F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
    }
    
    return coeffs;
}

double Transforms::evaluateCDF(const Coefficients& coeffs, double x)
{
    /**
     * Fast CDF Evaluation from Precomputed Coefficients
     * 
     * CDF(x) ≈ F₀/2 · (x-a)/(b-a) + Σ[k=1 to N-1] Fₖ/ωₖ · sin(ωₖ(x-a))
     * 
     * Complexity: O(N) arithmetic operations (no ChF evaluations!)
     */
    
    double x_shifted = x - coeffs.a;
    
    // First term: F₀/2 · (x-a)
    double cdf = coeffs.F_k[0] / 2.0 * x_shifted;
    
    // Sum over k = 1, 2, ..., N-1
    for (size_t k = 1; k < coeffs.N; ++k) {
        cdf += coeffs.F_k[k] / coeffs.omega[k] * std::sin(coeffs.omega[k] * x_shifted);
    }
    
    return cdf;
}

double Transforms::evaluatePDF(const Coefficients& coeffs, double x)
{
    /**
     * Fast PDF Evaluation from Precomputed Coefficients
     * 
     * PDF(x) ≈ Σ[k=0 to N-1] Fₖ · cos(ωₖ(x-a))
     * where F₀ is halved (first term correction)
     * 
     * Complexity: O(N) arithmetic operations (no ChF evaluations!)
     */
    
    double x_shifted = x - coeffs.a;
    
    // First term is halved
    double pdf = coeffs.F_k[0] / 2.0 * std::cos(coeffs.omega[0] * x_shifted);
    
    // Sum over k = 1, 2, ..., N-1
    for (size_t k = 1; k < coeffs.N; ++k) {
        pdf += coeffs.F_k[k] * std::cos(coeffs.omega[k] * x_shifted);
    }
    
    return pdf;
}

std::pair<double, size_t> Transforms::invertCDF_Optimized(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, double p, size_t max_iter, double tol)
{
    /**
     * Optimized Newton-Raphson CDF Inversion
     * 
     * Key optimization: Precompute F_k coefficients ONCE, then use fast evaluation
     * in Newton iterations.
     * 
     * Performance gain:
     * - Original: O(N × iterations × 2) ChF evaluations
     * - Optimized: O(N) ChF evaluations + O(N × iterations) arithmetic
     * 
     * Algorithm:
     * 1. Precompute F_k coefficients (expensive, but done once)
     * 2. Fast initial guess using cached coefficients
     * 3. Newton-Raphson with fast CDF/PDF evaluation (cheap)
     */
    
    // Clamp probability to [0, 1]
    p = std::max(0.0, std::min(1.0, p));
    
    // Handle edge cases
    if (p <= 0.0) return {a, 0};
    if (p >= 1.0) return {b, 0};
    
    // ========================================================================
    // STEP 1: Precompute coefficients (expensive, but only once!)
    // ========================================================================
    
    Coefficients coeffs = precomputeCoefficients(a, b, N, chf);
    
    // ========================================================================
    // STEP 2: Find good initial guess using fast evaluation
    // ========================================================================
    
    const size_t num_initial_points = 30;
    std::vector<double> x_initial(num_initial_points);
    
    // Create uniform grid over [a, b]
    for (size_t i = 0; i < num_initial_points; ++i) {
        x_initial[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(num_initial_points - 1);
    }
    
    // Fast CDF evaluation at initial points (no ChF evaluations!)
    size_t best_idx = 0;
    double min_dist = std::abs(evaluateCDF(coeffs, x_initial[0]) - p);
    for (size_t i = 1; i < num_initial_points; ++i) {
        double cdf_i = evaluateCDF(coeffs, x_initial[i]);
        double dist = std::abs(cdf_i - p);
        if (dist < min_dist) {
            min_dist = dist;
            best_idx = i;
        }
    }
    
    // Initial guess for Newton-Raphson
    double x = x_initial[best_idx];
    
    // ========================================================================
    // STEP 3: Optimized Newton-Raphson iteration (fast evaluation!)
    // ========================================================================
    
    for (size_t iter = 0; iter < max_iter; ++iter) {
        // Fast CDF evaluation (O(N) arithmetic operations)
        double cdf_x = evaluateCDF(coeffs, x);
        
        // Compute residual: f(x) = F(x) - p
        double fx = cdf_x - p;
        
        // Check convergence: |F(x) - p| < tol
        if (std::abs(fx) < tol) {
            return {x, iter};                   // Return x and number of iterations
        }
        
        // Fast PDF evaluation (O(N) arithmetic operations)
        double pdf_x = evaluatePDF(coeffs, x);
        
        // Safeguard against division by very small PDF
        if (std::abs(pdf_x) < 1e-10) {
            pdf_x = (pdf_x >= 0.0) ? 1e-10 : -1e-10;
        }
        
        // Newton step: x_new = x - (F(x) - p) / PDF(x)
        double x_new = x - fx / pdf_x;
        
        // Keep within bounds [a, b]
        x_new = std::max(a, std::min(b, x_new));
        
        // Check for convergence in x
        if (std::abs(x_new - x) < tol) {
            return {x_new, iter + 1};         // Return x and number of iterations
        }
        
        // Update for next iteration
        x = x_new;
    }
    
    // If max iterations reached, return current best estimate
    return {x, max_iter};                        // Return x and number of iterations
}

// ============================================================================
// CIR Sampler Implementation
// ============================================================================

double CIRSampler::sampleCIR(double kappa, double vbar, double sigma,
                             double v_s, double s, double t, double u)
{
    /**
    * Exact CIR Sampling via Noncentral Chi-Squared Distribution
    * 
    * V(t) | V(s) ~ (1/2c) · χ²(δ, λ)
    * 
    * Accepts pre-generated uniform u in (0,1) to enable:
    * - Antithetic sampling: pairs (u, 1-u)
    * 
    * @param u Uniform random variable for inverse CDF
    * 
    * Reference: Cox, Ingersoll & Ross (1985); Glasserman (2003)
    */
    
    double tau = t - s;

    // Degrees of freedom: δ = 4κθ/σ²
    double delta = 4.0 * kappa * vbar / (sigma * sigma);
     
    // Scaling factor: c = 2κ/(σ²(1-exp(-κτ)))
    double exp_term = std::exp(-kappa * tau);
    double c = 2.0 * kappa / (sigma * sigma * (1.0 - exp_term));
    
    // Noncentrality parameter: λ = 2c·exp(-κτ)·V(s)
    double lambda = 2.0 * c * exp_term * v_s;
    
    // Sample from noncentral chi-squared distribution using Boost
    boost::math::non_central_chi_squared_distribution<double> ncx2(delta, lambda);
    
    // Inverse CDF to get sample (using provided uniform)
    double chi_sq_sample = boost::math::quantile(ncx2, u);
    
    // Scale by 1/(2c) to get V(t)
    double V_t = chi_sq_sample / (2.0 * c);
    
    return V_t;
}

// ============================================================================
// * Characteristic Function of Integrated Variance Implementation
// ============================================================================

// TODO: As static method?
std::complex<double> ChFIntegratedVariance::compute(double omega, double kappa, double vbar, double sigma, double v_s, double v_t, double tau)
{
    /**
     * Characteristic Function of ∫[s,t] V(u)du for CIR Process
     * 
     * Direct translation of Python implementation
     * Equivalent to: scipy.special.iv() in Python
     * 
     * φ(ω) = E[exp(iω·∫V(u)du)] for V following CIR
     * 
     * Note: Uses Boost for real Bessel, custom implementation for complex
     * (Boost::Math doesn't support complex-valued Bessel functions)
     * 
     * Reference: Broadie & Kaya (2006), Equation (27)
     */
    
    using namespace std::complex_literals;
    // ! WARNING: Branch cut risk with std::sqrt if argument crosses negative real axis.
    //   For Integrated Variance, usually safe, but need to monitor for discontinuities.

    // R = √(κ² - 2γ²iω)  [Note: sigma ≡ gamma in Python notation]
    std::complex<double> R = std::sqrt(
        kappa * kappa - 2.0 * sigma * sigma * 1i * omega
    );
    
    // d = 4κθ/γ²  (degrees of freedom parameter)
    double d = 4.0 * kappa * vbar / (sigma * sigma);
    
    // Bessel function order: ν = d/2 - 1
    double nu = 0.5 * d - 1.0;
    
    // Common exponentials
    double exp_kappa_tau = std::exp(-kappa * tau);         // Real exponential
    std::complex<double> exp_R_tau = std::exp(-R * tau);   // Complex exponential
    
    // temp1 = R·exp(-τ(R-κ)/2)·(1-exp(-κτ))/(κ(1-exp(-Rτ)))
    std::complex<double> temp1 = R * std::exp(-tau / 2.0 * (R - kappa)) * 
                                 (1.0 - exp_kappa_tau) / 
                                 (kappa * (1.0 - exp_R_tau));
    
    // temp2 = exp((vu+vt)/γ² · [κ(1+exp(-κτ))/(1-exp(-κτ)) - R(1+exp(-Rτ))/(1-exp(-Rτ))])
    double kappa_part = kappa * (1.0 + exp_kappa_tau) / (1.0 - exp_kappa_tau);
    std::complex<double> R_part = R * (1.0 + exp_R_tau) / (1.0 - exp_R_tau);
    
    std::complex<double> temp2 = std::exp(
        (v_s + v_t) / (sigma * sigma) * (kappa_part - R_part)
    );
    
    // Bessel function arguments
    double sqrt_prod = std::sqrt(v_t * v_s);  // √(vt·vu)
    
    // z_R = 4R·√(vu·vt)·exp(-Rτ/2)/(γ²(1-exp(-Rτ)))  [complex]
    std::complex<double> z_R = sqrt_prod * 4.0 * R * std::exp(-R * tau / 2.0) /
                                (sigma * sigma * (1.0 - exp_R_tau));
    
    // z_kappa = 4κ·√(vu·vt)·exp(-κτ/2)/(γ²(1-exp(-κτ)))  [real]
    double z_kappa = sqrt_prod * 4.0 * kappa * std::exp(-kappa * tau / 2.0) /
                     (sigma * sigma * (1.0 - exp_kappa_tau));
    
    // Modified Bessel functions
    // For real argument: use Boost (fast, accurate)
    // For complex argument: use custom implementation (Boost doesn't support complex)
    
    // temp4 = I_ν(z_kappa) - real argument, use Boost
    double temp4 = boost::math::cyl_bessel_i(nu, z_kappa);
    
    // temp3 = I_ν(z_R) - complex argument, use custom implementation
    std::complex<double> temp3 = modifiedBesselI(nu, z_R);
    
    // Final characteristic function
    std::complex<double> chf = temp1 * temp2 * temp3 / temp4;
    
    return chf;
}

// Private helper: Modified Bessel function I_ν(z) for complex argument
std::complex<double> ChFIntegratedVariance::modifiedBesselI(double nu, std::complex<double> z)
{
    /**
     * Modified Bessel Function I_ν(z) for complex z
     * Cannot use full Boost::Math because it does not support complex-valued Bessel functions
     * So we need to implement that case ourselves
     * 
     * Uses power series for small |z|, asymptotic expansion for large |z|
     * I_ν(z) = Σ[k=0,∞] (z/2)^(2k+ν) / (k! · Γ(k+ν+1))
     */
    
    if (std::abs(z) < 1e-10) {
        return (nu == 0.0) ? 1.0 : 0.0;
    }
    
    // For large |z|, use asymptotic expansion: I_ν(z) ≈ exp(z)/√(2πz)
    // Asymptotic expansion for large |z|
    // NOTE: Threshold 30.0 is heuristic.
    //  could show the transition zone - against Python's
    if (std::abs(z) > 30.0) {
        return std::exp(z) / std::sqrt(2.0 * PI * z);
    }
    
    // Power series expansion
    std::complex<double> sum = 0.0;
    std::complex<double> term = std::pow(z / 2.0, nu) / std::tgamma(nu + 1.0);
    sum += term;
    
    for (int k = 1; k < 50; ++k) {
        term *= (z * z) / (4.0 * k * (k + nu));
        sum += term;
        if (std::abs(term / sum) < 1e-12) break;  // Convergence check
    }
    
    return sum;
}

// ============================================================================
// * Newton-Raphson Method
// ============================================================================

