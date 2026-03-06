#include <boost/math/distributions/non_central_chi_squared.hpp> // for sampling from non-central chi-squared distribution for CIRSampler
#include <cmath>
#include <cppfm/utils/Utils.h>
#include <numbers>
#include <stdexcept>

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

	if (std::abs(x) < 0.42)
	{
		// ====================================================================
		// Central region: |u - 0.5| < 0.42  (i.e., 0.08 < u < 0.92)
		// ====================================================================
		double r = x * x;

		// Numerator: a₀ + a₁r + a₂r² + a₃r³
		static const double a0 = 2.50662823884;
		static const double a1 = -18.61500062529;
		static const double a2 = 41.39119773534;
		static const double a3 = -25.44106049637;

		// Denominator: 1 + b₀r + b₁r² + b₂r³ + b₃r⁴
		static const double b0 = -8.47351093090;
		static const double b1 = 23.08336743743;
		static const double b2 = -21.06224101826;
		static const double b3 = 3.13082909833;

		double num = a0 + r * (a1 + r * (a2 + r * a3));
		double den = 1.0 + r * (b0 + r * (b1 + r * (b2 + r * b3)));

		return x * num / den;
	}
	else
	{
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

std::vector<double> ThomasAlgorithm::solve(const std::vector<double> &lower,
										   const std::vector<double> &diag,
										   const std::vector<double> &upper,
										   const std::vector<double> &rhs)
{
	// Validate input sizes
	size_t n = rhs.size();

	if (lower.size() != n || diag.size() != n || upper.size() != n)
	{
		throw std::invalid_argument("TridiagonalSolver::solve: All input vectors must have the same size");
	}

	if (n == 0)
	{
		throw std::invalid_argument("TridiagonalSolver::solve: Cannot solve empty system");
	}

	if (n == 1)
	{
		// Trivial case: b[0]*x[0] = d[0]
		if (std::abs(diag[0]) < 1e-14)
		{
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
	std::vector<double> c_prime(n, 0.0); // Modified upper diagonal
	std::vector<double> d_prime(n, 0.0); // Modified RHS
	std::vector<double> x(n, 0.0);		 // Solution vector

	// ------------------------------------------------------------------------
	/// PROCESS: Mx = r -> LUx = r -> (1) Ly = r -> (2) Ux = y
	// Step 1: Forward Elimination
	// ------------------------------------------------------------------------
	// Transform the system to upper triangular form
	// Original:  b[i]*x[i] + c[i]*x[i+1] = d[i] - a[i]*x[i-1]
	// Modified:  x[i] + c'[i]*x[i+1] = d'[i]

	// First row (i = 0): b[0]*x[0] + c[0]*x[1] = d[0]
	if (std::abs(diag[0]) < 1e-14)
	{
		throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row 0");
	}
	c_prime[0] = upper[0] / diag[0]; // γ1 = c1/b1
	d_prime[0] = rhs[0] / diag[0];	 // ρ1 = r1/b1

	// Interior rows (i = 1, 2, ..., n-1)
	for (size_t i = 1; i < n; ++i)
	{
		// Denominator after eliminating x[i-1]
		double m = diag[i] - lower[i] * c_prime[i - 1];
		//          b[i] -  a[i] * γ[i-1]
		if (std::abs(m) < 1e-14)
		{
			throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row " + std::to_string(i));
		}

		c_prime[i] = upper[i] / m;
		// c[i] = c[i+1] / (b[i] - a[i] * γ[i-1])
		d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / m;
		// d'[i] = (d[i] - a[i] * ρ[i-1]) / (b[i] - a[i] * γ[i-1])
	}

	// ------------------------------------------------------------------------
	// Step 2: Back Substitution
	// ------------------------------------------------------------------------
	// Solve the upper triangular system from bottom to top
	// x[i] = d'[i] - c'[i]*x[i+1]

	// Last equation (i = n-1): x[n-1] = d'[n-1]
	x[n - 1] = d_prime[n - 1];

	// Work backwards (i = n-2, n-3, ..., 0)
	for (int i = static_cast<int>(n) - 2; i >= 0; --i)
	{
		x[i] = d_prime[i] - c_prime[i] * x[i + 1];
		//     ρ[i+1]     - γ[i+1]     * x[i+2]
		//  (Code i=k corresponds to PDF row k+1)
	}

	return x;
} // end of ThomasAlgorithm::solve

// ============================================================================
// ! Numerical Derivatives Implementation
// ============================================================================

double NumericalDerivatives::firstDerivative(std::function<double(double)> f, double x, double h)
{
	// h is absolute — callers are responsible for scaling
	return (f(x + h) - f(x - h)) / (2.0 * h);
}

double NumericalDerivatives::secondDerivative(std::function<double(double)> f, double x, double h)
{
	// h is absolute
	return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
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
