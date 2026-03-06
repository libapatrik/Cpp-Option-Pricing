#include <cmath>
#include <cppfm/pricers/COSPricer.h>
#include <numbers>

constexpr double PI = std::numbers::pi;

// ============================================================================
// COS Option Pricer
// ============================================================================

double COSPricer::chi(size_t k, double a, double b, double c, double d)
{
	// chi_k(c,d) from Fang & Oosterlee Eq. (22)
	// k=0: chi_0 = e^d - e^c
	// k!=0: [e^d(cos(w(d-a)) + w*sin(w(d-a))) - e^c(cos(w(c-a)) + w*sin(w(c-a)))] / (1 + w^2)
	double bma = b - a;

	if (k == 0)
	{
		return std::exp(d) - std::exp(c);
	}

	double w = k * PI / bma;
	double w2 = w * w;

	double d_shifted = d - a;
	double c_shifted = c - a;

	double result = (std::exp(d) * (std::cos(w * d_shifted) + w * std::sin(w * d_shifted)) -
					 std::exp(c) * (std::cos(w * c_shifted) + w * std::sin(w * c_shifted))) /
					(1.0 + w2);

	return result;
}

double COSPricer::psi(size_t k, double a, double b, double c, double d)
{
	// psi_k(c,d) from Fang & Oosterlee Eq. (23)
	// k=0: psi_0 = d - c
	// k!=0: [sin(w(d-a)) - sin(w(c-a))] / w
	double bma = b - a;

	if (k == 0)
	{
		return d - c;
	}

	double w = k * PI / bma;
	double d_shifted = d - a;
	double c_shifted = c - a;

	return (std::sin(w * d_shifted) - std::sin(w * c_shifted)) / w;
}

double COSPricer::price(double S0, double K, double r, double T,
						const std::function<std::complex<double>(double)> &chf,
						OptionType type, size_t N, double L, double sigma)
{
	// V = K*e^{-rT} * sum' Re[H_k] * V_k
	// put-call parity for OTM options (numerical stability)
	bool isCall = (type == OptionType::Call);

	if (T <= 0.0)
	{
		return isCall ? std::max(S0 - K, 0.0) : std::max(K - S0, 0.0);
	}

	bool isOTM = (isCall && K > S0) || (!isCall && K < S0);
	bool priceAsCall = isOTM ? !isCall : isCall;

	double x0 = std::log(S0 / K);

	double vol = (sigma > 0.0) ? sigma : 0.25;
	double halfWidth = L * vol * std::sqrt(T);
	double a = x0 - halfWidth;
	double b = x0 + halfWidth;
	double bma = b - a;

	double sum = 0.0;

	for (size_t k = 0; k < N; ++k)
	{
		double w = k * PI / bma;

		std::complex<double> phi_k = chf(w);
		std::complex<double> H_k = phi_k * std::exp(std::complex<double>(0.0, w * (x0 - a)));

		double V_k;
		if (priceAsCall)
		{
			V_k = (2.0 / bma) * (chi(k, a, b, 0.0, b) - psi(k, a, b, 0.0, b));
		}
		else
		{
			V_k = (2.0 / bma) * (-chi(k, a, b, a, 0.0) + psi(k, a, b, a, 0.0));
		}

		// first term halved (Fourier convention)
		if (k == 0)
		{
			V_k *= 0.5;
		}

		sum += std::real(H_k) * V_k;
	}

	double rawPrice = std::max(0.0, K * std::exp(-r * T) * sum);

	// put-call parity if we priced the opposite option
	if (isOTM)
	{
		double pv_strike = K * std::exp(-r * T);
		if (isCall)
		{
			rawPrice = rawPrice + S0 - pv_strike;
		}
		else
		{
			rawPrice = rawPrice - S0 + pv_strike;
		}
	}

	return std::max(0.0, rawPrice);
}

double COSPricer::callPrice(double S0, double K, double r, double T,
							const std::function<std::complex<double>(double)> &chf,
							size_t N, double L, double sigma)
{
	return price(S0, K, r, T, chf, OptionType::Call, N, L, sigma);
}

double COSPricer::putPrice(double S0, double K, double r, double T,
						   const std::function<std::complex<double>(double)> &chf,
						   size_t N, double L, double sigma)
{
	return price(S0, K, r, T, chf, OptionType::Put, N, L, sigma);
}

// =============================================================================
// Vectorized COS Pricing - Multiple strikes at same maturity
// =============================================================================

std::vector<double> COSPricer::prices(double S0, const std::vector<double> &strikes,
									  double r, double T,
									  const std::function<std::complex<double>(double)> &chf,
									  OptionType type, size_t N, double L, double sigma)
{
	// precompute phi(omega_k) once, reuse for all strikes
	// O(M * N) CF evals -> O(N) CF evals

	size_t M = strikes.size();
	std::vector<double> result(M);

	if (M == 0)
		return result;

	bool isCall = (type == OptionType::Call);
	if (T <= 0.0)
	{
		for (size_t j = 0; j < M; ++j)
		{
			result[j] = isCall ? std::max(S0 - strikes[j], 0.0)
							   : std::max(strikes[j] - S0, 0.0);
		}
		return result;
	}

	// common domain [a, b] covering all strikes
	double x0_min = std::numeric_limits<double>::max();
	double x0_max = std::numeric_limits<double>::lowest();

	for (double K : strikes)
	{
		double x0 = std::log(S0 / K);
		x0_min = std::min(x0_min, x0);
		x0_max = std::max(x0_max, x0);
	}

	double vol = (sigma > 0.0) ? sigma : 0.25;
	double sqrtT = std::sqrt(T);
	double halfWidth = L * vol * sqrtT;

	double a = x0_min - halfWidth;
	double b = x0_max + halfWidth;
	double bma = b - a;

	// precompute phi(omega_k)
	std::vector<std::complex<double>> phi_vals(N);
	for (size_t k = 0; k < N; ++k)
	{
		double w = k * PI / bma;
		phi_vals[k] = chf(w);
	}

	// precompute chi/psi for common bounds
	std::vector<double> chi_call(N), psi_call(N);
	std::vector<double> chi_put(N), psi_put(N);

	for (size_t k = 0; k < N; ++k)
	{
		chi_call[k] = chi(k, a, b, 0.0, b);
		psi_call[k] = psi(k, a, b, 0.0, b);
		chi_put[k] = chi(k, a, b, a, 0.0);
		psi_put[k] = psi(k, a, b, a, 0.0);
	}

	double discount = std::exp(-r * T);

	// price each strike
	for (size_t j = 0; j < M; ++j)
	{
		double K = strikes[j];
		double x0 = std::log(S0 / K);

		bool isOTM = (isCall && K > S0) || (!isCall && K < S0);
		bool priceAsCall = isOTM ? !isCall : isCall;

		double sum = 0.0;
		for (size_t k = 0; k < N; ++k)
		{
			double w = k * PI / bma;

			std::complex<double> H_k = phi_vals[k] *
									   std::exp(std::complex<double>(0.0, w * (x0 - a)));

			double V_k;
			if (priceAsCall)
			{
				V_k = (2.0 / bma) * (chi_call[k] - psi_call[k]);
			}
			else
			{
				V_k = (2.0 / bma) * (-chi_put[k] + psi_put[k]);
			}

			if (k == 0)
				V_k *= 0.5;

			sum += std::real(H_k) * V_k;
		}

		double rawPrice = std::max(0.0, K * discount * sum);

		if (isOTM)
		{
			double pv_strike = K * discount;
			if (isCall)
			{
				rawPrice = rawPrice + S0 - pv_strike;
			}
			else
			{
				rawPrice = rawPrice - S0 + pv_strike;
			}
		}

		result[j] = std::max(0.0, rawPrice);
	}

	return result;
}

std::vector<double> COSPricer::callPrices(double S0, const std::vector<double> &strikes,
										  double r, double T,
										  const std::function<std::complex<double>(double)> &chf,
										  size_t N, double L, double sigma)
{
	return prices(S0, strikes, r, T, chf, OptionType::Call, N, L, sigma);
}

std::vector<double> COSPricer::putPrices(double S0, const std::vector<double> &strikes,
										 double r, double T,
										 const std::function<std::complex<double>(double)> &chf,
										 size_t N, double L, double sigma)
{
	return prices(S0, strikes, r, T, chf, OptionType::Put, N, L, sigma);
}
