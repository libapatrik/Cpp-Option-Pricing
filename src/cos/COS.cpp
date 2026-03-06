#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <cppfm/cos/COS.h>
#include <numbers>

constexpr double PI = std::numbers::pi;

// ============================================================================
// COS Method Implementation
// ============================================================================

std::vector<double> Transforms::recoverCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, const std::vector<double> &x)
{
	// F(x) ~ 1/2 * F_0*(x-a)/(b-a) + sum[k=1..N-1] F_k/omega_k * sin(omega_k*(x-a))
	// F_k = (2/(b-a)) * Re[phi(omega_k) * exp(-i*omega_k*a)]

	std::vector<double> cdf(x.size(), 0.0);

	std::vector<double> F_k(N);
	std::vector<double> omega(N);

	for (size_t k = 0; k < N; ++k)
	{
		omega[k] = k * PI / (b - a);

		std::complex<double> chf_val = chf(omega[k]);
		std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -omega[k] * a));
		F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
	}

	for (size_t i = 0; i < x.size(); ++i)
	{
		double x_shifted = x[i] - a;

		// first term: F_0/2 * (x-a)
		cdf[i] = F_k[0] / 2.0 * x_shifted;

		for (size_t k = 1; k < N; ++k)
		{
			cdf[i] += F_k[k] / omega[k] * std::sin(omega[k] * x_shifted);
		}
	}

	return cdf;
}

std::vector<double> Transforms::recoverPDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, const std::vector<double> &x)
{
	// f(x) ~ sum[k=0..N-1] F_k * cos(omega_k*(x-a)), F_0 halved

	std::vector<double> pdf(x.size(), 0.0);

	std::vector<double> F_k(N);
	std::vector<double> omega(N);

	for (size_t k = 0; k < N; ++k)
	{
		omega[k] = k * PI / (b - a);

		std::complex<double> chf_val = chf(omega[k]);
		std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -omega[k] * a));
		F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
	}

	// first term halved
	F_k[0] *= 0.5;

	for (size_t i = 0; i < x.size(); ++i)
	{
		double x_shifted = x[i] - a;

		for (size_t k = 0; k < N; ++k)
		{
			pdf[i] += F_k[k] * std::cos(omega[k] * x_shifted);
		}
	}

	return pdf;
}

std::pair<double, size_t> Transforms::invertCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, double p, size_t max_iter, double tol)
{
	// Newton-Raphson: find x s.t. F(x) = p
	// 1. Eval CDF at grid to get initial guess
	// 2. Newton: x := x - (F(x) - p) / f(x)

	p = std::max(0.0, std::min(1.0, p));

	if (p <= 0.0)
		return {a, 0};
	if (p >= 1.0)
		return {b, 0};

	// initial guess from grid search
	const size_t num_initial_points = 30;
	std::vector<double> x_initial(num_initial_points);

	for (size_t i = 0; i < num_initial_points; ++i)
	{
		x_initial[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(num_initial_points - 1);
	}

	std::vector<double> cdf_initial = recoverCDF(a, b, N, chf, x_initial);

	size_t best_idx = 0;
	double min_dist = std::abs(cdf_initial[0] - p);
	for (size_t i = 1; i < num_initial_points; ++i)
	{
		double dist = std::abs(cdf_initial[i] - p);
		if (dist < min_dist)
		{
			min_dist = dist;
			best_idx = i;
		}
	}

	double x = x_initial[best_idx];

	// Newton iterations
	for (size_t iter = 0; iter < max_iter; ++iter)
	{
		std::vector<double> x_vec = {x};
		std::vector<double> cdf_x_vec = recoverCDF(a, b, N, chf, x_vec);
		double cdf_x = cdf_x_vec[0];

		double fx = cdf_x - p;

		if (std::abs(fx) < tol)
		{
			return {x, iter};
		}

		std::vector<double> pdf_x_vec = recoverPDF(a, b, N, chf, x_vec);
		double pdf_x = pdf_x_vec[0];

		if (std::abs(pdf_x) < 1e-10)
		{
			pdf_x = (pdf_x >= 0.0) ? 1e-10 : -1e-10;
		}

		double x_new = x - fx / pdf_x;

		x_new = std::max(a, std::min(b, x_new));

		if (std::abs(x_new - x) < tol)
		{
			return {x_new, iter + 1};
		}

		x = x_new;
	}

	return {x, max_iter};
}

// ============================================================================
// Optimized Newton-Raphson with Coefficient Caching
// ============================================================================

Transforms::Coefficients Transforms::precomputeCoefficients(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf)
{
	// precompute F_k = (2/(b-a)) * Re[phi(omega_k) * exp(-i*omega_k*a)]
	// this is O(N) CF evaluations, then CDF/PDF eval is just arithmetic

	Coefficients coeffs;
	coeffs.a = a;
	coeffs.b = b;
	coeffs.N = N;
	coeffs.F_k.resize(N);
	coeffs.omega.resize(N);

	for (size_t k = 0; k < N; ++k)
	{
		coeffs.omega[k] = k * PI / (b - a);

		std::complex<double> chf_val = chf(coeffs.omega[k]);
		std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -coeffs.omega[k] * a));
		coeffs.F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);
	}

	return coeffs;
}

double Transforms::evaluateCDF(const Coefficients &coeffs, double x)
{
	// CDF(x) ~ F_0/2 * (x-a) + sum[k=1..N-1] F_k/omega_k * sin(omega_k*(x-a))

	double x_shifted = x - coeffs.a;

	double cdf = coeffs.F_k[0] / 2.0 * x_shifted;

	for (size_t k = 1; k < coeffs.N; ++k)
	{
		cdf += coeffs.F_k[k] / coeffs.omega[k] * std::sin(coeffs.omega[k] * x_shifted);
	}

	return cdf;
}

double Transforms::evaluatePDF(const Coefficients &coeffs, double x)
{
	// PDF(x) ~ sum[k=0..N-1] F_k * cos(omega_k*(x-a)), F_0 halved

	double x_shifted = x - coeffs.a;

	double pdf = coeffs.F_k[0] / 2.0 * std::cos(coeffs.omega[0] * x_shifted);

	for (size_t k = 1; k < coeffs.N; ++k)
	{
		pdf += coeffs.F_k[k] * std::cos(coeffs.omega[k] * x_shifted);
	}

	return pdf;
}

std::pair<double, size_t> Transforms::invertCDF_Optimized(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, double p, size_t max_iter, double tol)
{
	// same as invertCDF but precomputes coefficients once
	// O(N) CF evals + O(N * iters) arithmetic vs O(N * iters * 2) CF evals

	p = std::max(0.0, std::min(1.0, p));

	if (p <= 0.0)
		return {a, 0};
	if (p >= 1.0)
		return {b, 0};

	// precompute once
	Coefficients coeffs = precomputeCoefficients(a, b, N, chf);

	// initial guess from grid
	const size_t num_initial_points = 30;
	std::vector<double> x_initial(num_initial_points);

	for (size_t i = 0; i < num_initial_points; ++i)
	{
		x_initial[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(num_initial_points - 1);
	}

	// fast CDF eval at initial points (no CF evaluations!)
	size_t best_idx = 0;
	double min_dist = std::abs(evaluateCDF(coeffs, x_initial[0]) - p);
	for (size_t i = 1; i < num_initial_points; ++i)
	{
		double cdf_i = evaluateCDF(coeffs, x_initial[i]);
		double dist = std::abs(cdf_i - p);
		if (dist < min_dist)
		{
			min_dist = dist;
			best_idx = i;
		}
	}

	double x = x_initial[best_idx];

	// Newton iterations with fast evaluation
	for (size_t iter = 0; iter < max_iter; ++iter)
	{
		double cdf_x = evaluateCDF(coeffs, x);

		double fx = cdf_x - p;

		if (std::abs(fx) < tol)
		{
			return {x, iter};
		}

		double pdf_x = evaluatePDF(coeffs, x);

		if (std::abs(pdf_x) < 1e-10)
		{
			pdf_x = (pdf_x >= 0.0) ? 1e-10 : -1e-10;
		}

		double x_new = x - fx / pdf_x;

		x_new = std::max(a, std::min(b, x_new));

		if (std::abs(x_new - x) < tol)
		{
			return {x_new, iter + 1};
		}

		x = x_new;
	}

	return {x, max_iter};
}

// ============================================================================
// Characteristic Function of Integrated Variance
// ============================================================================

std::complex<double> ChFIntegratedVariance::compute(double omega, double kappa, double vbar, double sigma, double v_s, double v_t, double tau)
{
	// phi(omega) = E[exp(i*omega * integral V(u)du)] for CIR V
	// Broadie & Kaya (2006), Equation (27)

	using namespace std::complex_literals;
	// ! branch cut risk with std::sqrt if argument crosses negative real axis

	// R = sqrt(kappa^2 - 2*sigma^2*i*omega)
	std::complex<double> R = std::sqrt(
		kappa * kappa - 2.0 * sigma * sigma * 1i * omega);

	// d = 4*kappa*vbar/sigma^2 (degrees of freedom)
	double d = 4.0 * kappa * vbar / (sigma * sigma);

	// Bessel order nu = d/2 - 1
	double nu = 0.5 * d - 1.0;

	double exp_kappa_tau = std::exp(-kappa * tau);
	std::complex<double> exp_R_tau = std::exp(-R * tau);

	// temp1 = R*exp(-tau*(R-kappa)/2)*(1-exp(-kappa*tau)) / (kappa*(1-exp(-R*tau)))
	std::complex<double> temp1 = R * std::exp(-tau / 2.0 * (R - kappa)) *
								 (1.0 - exp_kappa_tau) /
								 (kappa * (1.0 - exp_R_tau));

	// temp2 = exp((vs+vt)/sigma^2 * [kappa*(1+exp(-kappa*tau))/(1-exp(-kappa*tau)) - R*(1+exp(-R*tau))/(1-exp(-R*tau))])
	double kappa_part = kappa * (1.0 + exp_kappa_tau) / (1.0 - exp_kappa_tau);
	std::complex<double> R_part = R * (1.0 + exp_R_tau) / (1.0 - exp_R_tau);

	std::complex<double> temp2 = std::exp(
		(v_s + v_t) / (sigma * sigma) * (kappa_part - R_part));

	// Bessel function arguments
	double sqrt_prod = std::sqrt(v_t * v_s);

	// z_R = 4*R*sqrt(vs*vt)*exp(-R*tau/2) / (sigma^2*(1-exp(-R*tau)))
	std::complex<double> z_R = sqrt_prod * 4.0 * R * std::exp(-R * tau / 2.0) /
							   (sigma * sigma * (1.0 - exp_R_tau));

	// z_kappa = 4*kappa*sqrt(vs*vt)*exp(-kappa*tau/2) / (sigma^2*(1-exp(-kappa*tau)))
	double z_kappa = sqrt_prod * 4.0 * kappa * std::exp(-kappa * tau / 2.0) /
					 (sigma * sigma * (1.0 - exp_kappa_tau));

	// I_nu(z_kappa) - real arg, use Boost
	double temp4 = boost::math::cyl_bessel_i(nu, z_kappa);

	// I_nu(z_R) - complex arg, custom implementation
	std::complex<double> temp3 = modifiedBesselI(nu, z_R);

	std::complex<double> chf = temp1 * temp2 * temp3 / temp4;

	return chf;
}

// I_nu(z) for complex z
// power series for small |z|, asymptotic for large |z|
std::complex<double> ChFIntegratedVariance::modifiedBesselI(double nu, std::complex<double> z)
{
	if (std::abs(z) < 1e-10)
	{
		return (nu == 0.0) ? 1.0 : 0.0;
	}

	// asymptotic: I_nu(z) ~ exp(z)/sqrt(2*pi*z) * [1 - (4*nu^2-1)/(8z)]
	if (std::abs(z) > 30.0)
	{
		double mu = 4.0 * nu * nu;
		return std::exp(z) / std::sqrt(2.0 * PI * z) * (1.0 - (mu - 1.0) / (8.0 * z));
	}

	// power series: I_nu(z) = sum (z/2)^(2k+nu) / (k! * Gamma(k+nu+1))
	std::complex<double> sum = 0.0;
	std::complex<double> term = std::pow(z / 2.0, nu) / std::tgamma(nu + 1.0);
	sum += term;

	for (int k = 1; k < 50; ++k)
	{
		term *= (z * z) / (4.0 * k * (k + nu));
		sum += term;
		if (std::abs(term / sum) < 1e-12)
			break;
	}

	return sum;
}

// ============================================================================
// Heston Characteristic Function
// ============================================================================

HestonCF::HestonCF(double kappa, double vbar, double sigma, double rho, double v0, double r, double T)
	: _kappa(kappa), _vbar(vbar), _sigma(sigma), _rho(rho), _v0(v0), _r(r), _T(T)
{
}

std::complex<double> HestonCF::operator()(double u) const
{
	// phi(u) = exp(C + D*v0 + i*u*r*T)
	// "Little Heston Trap" formulation, forces Re(d) > 0
	using namespace std::complex_literals;

	std::complex<double> iu = 1i * u;

	// d = sqrt((kappa - rho*sigma*i*u)^2 + sigma^2*(i*u + u^2))
	std::complex<double> term1 = _kappa - _rho * _sigma * iu;
	std::complex<double> term2 = _sigma * _sigma * (iu + u * u);
	std::complex<double> d = std::sqrt(term1 * term1 + term2);

	// force Re(d) > 0 for stable branch (Albrecher et al. 2007)
	if (std::real(d) < 0.0)
	{
		d = -d;
	}

	// g = (kappa - rho*sigma*i*u - d) / (kappa - rho*sigma*i*u + d)
	std::complex<double> g = (term1 - d) / (term1 + d);

	std::complex<double> exp_neg_dT = std::exp(-d * _T);

	// C = (kappa*vbar/sigma^2) * [(kappa - rho*sigma*i*u - d)*T - 2*ln((1 - g*e^{-dT})/(1-g))]
	std::complex<double> C = (_kappa * _vbar / (_sigma * _sigma)) * ((term1 - d) * _T - 2.0 * std::log((1.0 - g * exp_neg_dT) / (1.0 - g)));

	// D = ((kappa - rho*sigma*i*u - d)/sigma^2) * (1 - e^{-dT})/(1 - g*e^{-dT})
	std::complex<double> D = ((term1 - d) / (_sigma * _sigma)) *
							 (1.0 - exp_neg_dT) / (1.0 - g * exp_neg_dT);

	return std::exp(C + D * _v0 + iu * _r * _T);
}
