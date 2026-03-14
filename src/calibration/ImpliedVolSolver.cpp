//
// Brent + Newton implied vol solver
//

#include <algorithm>
#include <cmath>
#include <cppfm/calibration/ImpliedVolSolver.h>
#include <cppfm/pricers/BlackScholesFormulas.h>

double ImpliedVolSolver::solve(double price, double S0, double K, double r, double T,
							   bool isCall, double tol, size_t maxIter)
{
	if (T <= 0.0)
		return 0.0;

	// intrinsic bound
	double df = std::exp(-r * T);
	double intrinsic = isCall ? std::max(S0 - K * df, 0.0) : std::max(K * df - S0, 0.0);
	if (price <= intrinsic + tol)
		return 0.0;

	auto type = isCall ? Option::Type::Call : Option::Type::Put;
	const double lo = 1e-6, hi = 5.0;

	// Newton with bisection fallback
	double sigma = 0.3;
	double a = lo, b = hi;
	double fa = BlackScholesFormulas::price(S0, K, r, a, T, type) - price;

	for (size_t i = 0; i < maxIter; ++i)
	{
		double p = BlackScholesFormulas::price(S0, K, r, sigma, T, type);
		double diff = p - price;

		if (std::abs(diff) < tol)
			return sigma;

		// try Newton step
		double v = BlackScholesFormulas::vega(S0, K, r, sigma, T);
		double newSigma = sigma;
		if (v > 1e-12)
		{
			newSigma = sigma - diff / v;
		}

		// accept Newton if in bounds, otherwise bisect
		if (newSigma > lo && newSigma < hi)
		{
			sigma = newSigma;
		}
		else
		{
			if (diff * fa < 0)
			{
				b = sigma;
			}
			else
			{
				a = sigma;
				fa = diff;
			}
			sigma = 0.5 * (a + b);
		}
	}

	return sigma;
}
