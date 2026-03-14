#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <cppfm/cos/COS.h>
#include <numbers>


// ============================================================================
// COS Method Implementation
// ============================================================================

std::vector<double> Transforms::recoverCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, const std::vector<double> &x)
{
	auto coeffs = precomputeCoefficients(a, b, N, chf);
	std::vector<double> cdf(x.size());
	for (size_t i = 0; i < x.size(); ++i)
		cdf[i] = evaluateCDF(coeffs, x[i]);
	return cdf;
}

std::vector<double> Transforms::recoverPDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, const std::vector<double> &x)
{
	auto coeffs = precomputeCoefficients(a, b, N, chf);
	std::vector<double> pdf(x.size());
	for (size_t i = 0; i < x.size(); ++i)
		pdf[i] = evaluatePDF(coeffs, x[i]);
	return pdf;
}

// ============================================================================
// Newton-Raphson with Coefficient Caching
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
		coeffs.omega[k] = k * std::numbers::pi / (b - a);

		std::complex<double> chf_val = chf(coeffs.omega[k]);
		std::complex<double> exp_term = std::exp(std::complex<double>(0.0, -coeffs.omega[k] * a));
		coeffs.F_k[k] = 2.0 / (b - a) * std::real(chf_val * exp_term);

		// adaptive truncation: stop when coefficients have decayed
		if (k > 16 && std::abs(coeffs.F_k[k]) < 1e-12 && std::abs(coeffs.F_k[k - 1]) < 1e-12)
		{
			coeffs.N = k + 1;
			coeffs.F_k.resize(coeffs.N);
			coeffs.omega.resize(coeffs.N);
			break;
		}
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

std::pair<double, size_t> Transforms::invertCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)> &chf, double p, size_t max_iter, double tol)
{
	// Newton-Raphson: find x s.t. F(x) = p
	// precompute coefficients once, then iterate with fast CDF/PDF evaluation

	p = std::max(0.0, std::min(1.0, p));

	if (p <= 0.0)
		return {a, 0};
	if (p >= 1.0)
		return {b, 0};

	// precompute once
	Coefficients coeffs = precomputeCoefficients(a, b, N, chf);

	// bisection initial guess: 5 CDF evals → bracket of (b-a)/32
	double lo = a, hi = b;
	for (int i = 0; i < 5; ++i)
	{
		double mid = 0.5 * (lo + hi);
		if (evaluateCDF(coeffs, mid) < p)
			lo = mid;
		else
			hi = mid;
	}
	double x = 0.5 * (lo + hi);

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
// Truncation Bounds
// ============================================================================

std::pair<double, double> computeTruncationBounds(double c1, double c2, double c4, double L)
{
	// Fang & Oosterlee (2008) Eq. 23
	double width = L * std::sqrt(c2 + std::sqrt(std::max(c4, 0.0)));
	return {c1 - width, c1 + width};
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

	// Hankel asymptotic: I_nu(z) ~ e^z/sqrt(2*pi*z) * sum_k a_k/z^k
	// a_k = prod_{j=0}^{k-1} (4*nu^2 - (2j+1)^2) / (8*(j+1))
	// 4 terms: error ~ O(1/z^4) ≈ 1e-7 at |z|=30
	if (std::abs(z) > 30.0 && std::real(z) > 0.0)
	{
		double mu = 4.0 * nu * nu;
		std::complex<double> iz = 1.0 / z;
		// recurrence: a_{k+1} = a_k * -(mu - (2k+1)^2) / (8*(k+1)) / z
		std::complex<double> a = 1.0;
		std::complex<double> s = a;
		for (int k = 0; k < 6; ++k)
		{
			double odd = 2.0 * k + 1.0;
			a *= -(mu - odd * odd) * iz / (8.0 * (k + 1));
			s += a;
		}
		return std::exp(z) / std::sqrt(2.0 * std::numbers::pi * z) * s;
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

// ============================================================================
// Heston Analytical Cumulants
// Le Floc'h (2018) — verified against algorithmic differentiation
// Corrects known error in Fang & Oosterlee (2008) Table 2 c2 formula
// ============================================================================

HestonCF::Cumulants HestonCF::cumulants() const
{
	// match Le Floc'h variable names for translation
	double lambda = _kappa;
	double ubar = _vbar;
	double u0 = _v0;
	double eta = _sigma;
	double rho = _rho;
	double term = _T;

	// degenerate: constant variance, no vol-of-vol
	if (lambda == 0.0 && eta == 0.0)
	{
		double c1 = _r * term - u0 * term / 2.0;
		double c2 = u0 * term;
		return {c1, c2, 0.0};
	}

	// no mean reversion
	if (lambda == 0.0)
	{
		double c1 = _r * term - u0 * term / 2.0;
		double c2 = u0 * term * (1.0 + term * term * eta * (eta * term / 12.0 - rho / 2.0));
		double eta2 = eta * eta;
		double eta3 = eta2 * eta;
		double eta4 = eta2 * eta2;
		double term2 = term * term;
		double term3 = term2 * term;
		double term4 = term3 * term;
		double term5 = term4 * term;
		double rho2 = rho * rho;
		double rho3 = rho2 * rho;
		double c4 = u0 * eta2 * term3 * (2.0 * rho2 - 2.0 * eta * term * rho + eta4 * term4 * 17.0 / 1680.0 + eta4 * term5 * rho2 * 11.0 / 20.0 - eta * term * rho3 / 2.0 + eta2 * term2 * 3.0 / 10.0 - eta3 * term3 * rho * 17.0 / 120.0 + 1.0);
		return {c1, c2, c4};
	}

	// general case
	double elt = std::exp(-lambda * term);

	double c1 = _r * term + (1.0 - elt) * (ubar - u0) / (2.0 * lambda) - ubar * term / 2.0;

	double lambda2 = lambda * lambda;
	double lambda3 = lambda * lambda2;
	double eta2 = eta * eta;
	double elt2 = elt * elt;

	double c2 = 1.0 / (8.0 * lambda3) * (2.0 * u0 * (lambda2 * 4.0 * (elt * term * rho * eta + 1.0 - elt) + lambda * (4.0 * eta * rho * (elt - 1.0) - 2.0 * elt * term * eta2) + eta2 * (1.0 - elt2)) + ubar * (8.0 * lambda3 * term + lambda2 * 8.0 * (-eta * rho * term * (1.0 + elt) + elt - 1.0) + lambda * (2.0 * eta2 * term * (1.0 + 2.0 * elt) + eta * rho * 16.0 * (1.0 - elt)) + eta2 * (-5.0 + 4.0 * elt + elt2)));

	// c4 — full analytical, separate theta and v0 contributions
	double lambda4 = lambda2 * lambda2;
	double lambda5 = lambda4 * lambda;
	double lambda6 = lambda4 * lambda2;
	double lambda7 = lambda3 * lambda4;
	double eta3 = eta * eta2;
	double eta4 = eta2 * eta2;
	double term2 = term * term;
	double term3 = term2 * term;
	double rho2 = rho * rho;
	double rho3 = rho2 * rho;
	double elt3 = elt2 * elt;
	double elt4 = elt2 * elt2;

	// theta (ubar) contribution
	double c4 = -1.0 / lambda7 * 2.0 * ubar * eta2 * (((term3 * rho3 * eta - 3.0 * term2 * rho2) * lambda6 - 3.0 * term * (term2 * rho2 * eta2 - 4.0 * term * rho * (rho2 + 1.0) * eta + 8.0 * rho2 + 2.0) * lambda5 / 2.0 + (3.0 * term3 * rho * eta3 / 4.0 - 21.0 * term2 * (rho2 + 3.0 / 14.0) * eta2 / 2.0 + (18.0 * term * rho3 + 24.0 * term * rho) * eta - 18.0 * rho2 - 3.0) * lambda4 - (term3 * eta3 - 42.0 * term2 * rho * eta2 + (240.0 * term * rho2 + 54.0 * term) * eta - 192.0 * rho3 - 192.0 * rho) * eta * lambda3 / 8.0 - 3.0 * eta2 * (term2 * eta2 - 35.0 * term / 2.0 * rho * eta + 40.0 * rho2 + 15.0 / 2.0) * lambda2 / 4.0 - 27.0 * eta3 * (term * eta - 20.0 * rho / 3.0) * lambda / 16.0 - 21.0 * eta4 / 16.0) * elt + ((-3.0 / 4.0 + 3.0 * term * rho * eta - 3.0 * term2 * rho2 * eta2 / 2.0) * lambda4 + 3.0 * eta * (term2 * rho * eta2 + (-4.0 * term * rho2 - 3.0 * term / 2.0) * eta + 4.0 * rho) * lambda3 / 2.0 - 3.0 * eta2 * (term2 * eta2 - 14.0 * term * rho * eta + 20.0 * rho2 + 6.0) * lambda2 / 8.0 + (-15.0 / 16.0 * term * eta4 + 9.0 * eta3 * rho / 2.0) * lambda - 21.0 * eta4 / 32.0) * elt2 + 3.0 * eta2 * ((term * rho * eta - 1.0) * lambda2 + (-term / 2.0 * eta2 + 2.0 * rho * eta) * lambda - eta2 / 2.0) * elt3 / 8.0 - 3.0 * eta4 * elt4 / 128.0 + (-6.0 * term * rho2 - 3.0 * term / 2.0) * lambda5 + ((6.0 * term * rho3 + 9.0 * term * rho) * eta + 18.0 * rho2 + 15.0 / 4.0) * lambda4 - 9.0 * eta * ((rho2 + 0.25) * term * eta + 8.0 * rho3 / 3.0 + 10.0 * rho / 3.0) * lambda3 + 15.0 * eta2 * (term * rho * eta + 10.0 * rho2 + 11.0 / 5.0) * lambda2 / 4.0 + (-33.0 / 2.0 * eta3 * rho - 15.0 / 32.0 * term * eta4) * lambda + 279.0 * eta4 / 128.0);

	// v0 contribution
	c4 += u0 / lambda7 * (2.0 * eta2 * (((term3 * rho3 * eta - 3.0 * term2 * rho2) * lambda6 - 3.0 * term * (term2 * rho2 * eta2 - 2.0 * term * rho * (rho2 + 2.0) * eta + 4.0 * rho2 + 2.0) * lambda5 / 2.0 + (3.0 * term3 * rho * eta3 / 4.0 - 6.0 * (rho2 + 3.0 / 8.0) * term2 * eta2 + 6.0 * term * rho * (rho2 + 2.0) * eta - 6.0 * rho2) * lambda4 - eta * (term3 * eta3 - 24.0 * term2 * rho * eta2 + (72.0 * term * rho2 + 18.0 * term) * eta - 48.0 * rho3) * lambda3 / 8.0 - 3.0 * eta2 * (term2 * eta2 - 7.0 * term * rho * eta - 3.0) * lambda2 / 8.0 - 3.0 * eta3 * (term * eta + 10.0 * rho) * lambda / 16.0 + 3.0 * eta4 / 8.0) * elt + ((6.0 * term * rho * eta - 3.0 * term2 * rho2 * eta2 - 3.0 / 2.0) * lambda4 + 3.0 * (term2 * rho * eta2 + (-3.0 * term * rho2 - 3.0 * term / 2.0) * eta + 3.0 * rho) * eta * lambda3 - 3.0 * eta2 * (term2 * eta2 - 10.0 * term * rho * eta + 12.0 * rho2 + 3.0) * lambda2 / 4.0 - 9.0 * eta3 * (term * eta - 10.0 * rho / 3.0) * lambda / 8.0 - 3.0 * eta4 / 8.0) * elt2 + 9.0 * eta2 * ((term * rho * eta - 1.0) * lambda2 + (-term / 2.0 * eta2 + 5.0 / 3.0 * rho * eta) * lambda - eta2 / 3.0) * elt3 / 8.0 - 3.0 * eta4 * elt4 / 32.0 - 6.0 * ((rho2 + 1.0 / 4.0) * lambda2 - 5.0 * lambda * rho * eta / 4.0 + 5.0 * eta2 / 16.0) * (lambda * rho * eta - eta2 / 4.0 - lambda2)));

	return {c1, c2, c4};
}

// ============================================================================
// ChFIntegratedVariance Numerical Cumulants
// c1, c2 via finite differences at h=7e-3; c4=0 (Bessel noise makes it unreliable)
// ============================================================================

ChFIntegratedVariance::Cumulants ChFIntegratedVariance::cumulants(
	double kappa, double vbar, double sigma, double v_s, double v_t, double tau)
{
	double h = 7e-3;

	auto lnphi = [&](double omega) -> std::complex<double>
	{
		return std::log(compute(omega, kappa, vbar, sigma, v_s, v_t, tau));
	};

	std::complex<double> lp1 = lnphi(h);
	std::complex<double> lm1 = lnphi(-h);

	// c1 = Im[(lp1 - lm1) / (2h)]  (κ₁ = (-i)·d/du ln φ, Re(-iz) = Im(z))
	double c1 = std::imag(lp1 - lm1) / (2.0 * h);

	// c2 = -Re[(lp1 + lm1) / h^2]  (using ln(phi(0)) = 0)
	double c2 = -std::real(lp1 + lm1) / (h * h);

	// Bessel noise can make |phi| > 1 → Re(ln phi) > 0 → c2 < 0
	// fall back to Var[∫V ds] ≈ σ²·V_avg·τ³/3 (CIR moment approximation)
	if (c2 <= 0.0)
	{
		double V_avg = 0.5 * (v_s + v_t);
		c2 = sigma * sigma * V_avg * tau * tau * tau / 3.0;
	}

	return {c1, c2, 0.0};
}
