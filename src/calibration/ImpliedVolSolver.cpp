#include <algorithm>
#include <cmath>
#include <cppfm/calibration/ImpliedVolSolver.h>
#include <cppfm/pricers/BlackScholesFormulas.h>

// TODO: Add LetsBeRational by Jaeckel

double ImpliedVolSolver::solve(double price, double S0, double K, double r, double T,
							   bool isCall, double tol, size_t maxIter)
{
	// Brent's method with Newton acceleration for IV inversion
	if (T <= 0.0)
		return 0.0;

	double intrinsic = isCall ? std::max(S0 - K * std::exp(-r * T), 0.0)
							  : std::max(K * std::exp(-r * T) - S0, 0.0);
	if (price <= intrinsic + 1e-10)
	{
		return 0.001;
	}

	// BS price minus target
	auto bsObj = [&](double vol) -> double
	{
		double c = BlackScholesFormulas::callPrice(S0, K, r, vol, T);
		double p = isCall ? c : (c - S0 + K * std::exp(-r * T));
		return p - price;
	};

	// bracket [lo, hi]
	double lo = 0.001, hi = 5.0;
	double flo = bsObj(lo), fhi = bsObj(hi);

	// degenerate: price above BS(5.0) or below BS(0.001)
	if (flo * fhi > 0.0)
		return (flo > 0.0) ? lo : hi;

	// Brent's method
	double a = lo, b = hi, fa = flo, fb = fhi;
	double c = a, fc = fa;
	double d = b - a, e = d;

	for (size_t iter = 0; iter < maxIter; ++iter)
	{
		if (std::abs(fb) < tol)
			return b;

		if (fb * fc > 0.0)
		{
			c = a;
			fc = fa;
			d = b - a;
			e = d;
		}
		if (std::abs(fc) < std::abs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		double tolBrent = 2.0 * 1e-15 * std::abs(b) + 0.5 * tol;
		double m = 0.5 * (c - b);

		if (std::abs(m) <= tolBrent)
			return b;

		// try Newton if vega is decent
		double vega = BlackScholesFormulas::vega(S0, K, r, b, T);
		bool useNewton = false;
		double newtonStep = 0.0;
		if (vega > 1e-10)
		{
			newtonStep = -fb / vega;
			double candidate = b + newtonStep;
			// accept Newton only if it stays in bracket
			if (candidate > std::min(b, c) && candidate < std::max(b, c))
				useNewton = true;
		}

		if (useNewton)
		{
			d = newtonStep;
			e = d;
		}
		else if (std::abs(e) >= tolBrent && std::abs(fa) > std::abs(fb))
		{
			// inverse quadratic interpolation or secant
			double s = fb / fa;
			double p, q;
			if (std::abs(a - c) < 1e-15)
			{
				p = 2.0 * m * s;
				q = 1.0 - s;
			}
			else
			{
				double qq = fa / fc;
				double rr = fb / fc;
				p = s * (2.0 * m * qq * (qq - rr) - (b - a) * (rr - 1.0));
				q = (qq - 1.0) * (rr - 1.0) * (s - 1.0);
			}
			if (p > 0.0)
				q = -q;
			else
				p = -p;
			if (2.0 * p < std::min(3.0 * m * q - std::abs(tolBrent * q), std::abs(e * q)))
			{
				e = d;
				d = p / q;
			}
			else
			{
				d = m;
				e = m;
			}
		}
		else
		{
			d = m;
			e = m;
		}

		a = b;
		fa = fb;
		if (std::abs(d) > tolBrent)
			b += d;
		else
			b += (m > 0.0 ? tolBrent : -tolBrent);
		fb = bsObj(b);
	}

	return b;
}
