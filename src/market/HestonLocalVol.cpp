#include <cmath>
#include <cppfm/cos/COS.h>
#include <cppfm/market/HestonLocalVol.h>
#include <cppfm/pricers/COSPricer.h>

HestonLocalVol::HestonLocalVol(double kappa, double vbar, double sigma, double rho,
							   double v0, double r, double S0)
	: _kappa(kappa), _vbar(vbar), _sigma(sigma), _rho(rho), _v0(v0), _r(r), _S0(S0)
{
}

void HestonLocalVol::setCOSParams(size_t N, double L)
{
	_N = N;
	_L = L;
}

double HestonLocalVol::callPrice(double K, double T) const
{
	// delegate to batch method
	return callPrices({K}, T)[0];
}

double HestonLocalVol::localVariance(double K, double T) const
{
	double lv = localVol(K, T);
	return lv * lv;
}

double HestonLocalVol::localVol(double K, double T) const
{
	return localVolsAtTime({K}, T)[0];
}

// =============================================================================
// Vectorized Batch Computation
// =============================================================================

std::vector<double> HestonLocalVol::callPrices(const std::vector<double> &strikes, double T) const
{
	HestonCF cf(_kappa, _vbar, _sigma, _rho, _v0, _r, T);
	double vol_estimate = std::sqrt(_v0);
	return COSPricer::callPrices(_S0, strikes, _r, T, cf, _N, _L, vol_estimate);
}

std::vector<double> HestonLocalVol::localVolsAtTime(const std::vector<double> &spots, double T) const
{
	// batch Dupire local vol using vectorized COS pricing
	// for each spot, need 5 strikes at T (for dCdK, d2CdK2)
	// and central strike at T+-h (for dCdT)
	// Richardson extrapolation on all derivatives

	size_t M = spots.size();
	std::vector<double> result(M);

	if (M == 0)
		return result;

	if (T <= 1e-6)
	{
		double sqrtV0 = std::sqrt(_v0);
		std::fill(result.begin(), result.end(), sqrtV0);
		return result;
	}

	// build strike grids
	std::vector<double> h1_vec(M), h2_vec(M);
	std::vector<double> allStrikes_T;
	allStrikes_T.reserve(5 * M);

	for (size_t j = 0; j < M; ++j)
	{
		double K = spots[j];
		double h1 = K * 0.005; // 0.5% of strike
		double h2 = h1 / 2.0;

		h1_vec[j] = h1;
		h2_vec[j] = h2;

		// 5 strikes per spot: K-h1, K-h2, K, K+h2, K+h1
		allStrikes_T.push_back(K - h1); // index 5*j + 0
		allStrikes_T.push_back(K - h2); // index 5*j + 1
		allStrikes_T.push_back(K);		// index 5*j + 2
		allStrikes_T.push_back(K + h2); // index 5*j + 3
		allStrikes_T.push_back(K + h1); // index 5*j + 4
	}

	// central strikes for dCdT
	std::vector<double> centralStrikes(M);
	for (size_t j = 0; j < M; ++j)
	{
		centralStrikes[j] = spots[j];
	}

	// time derivative step sizes
	double hT1 = std::max(T * 0.02, 0.002); // ~2% or minimum 2 days
	double hT2 = hT1 / 2.0;

	bool useForwardDiff = (T - hT1 <= 0.0);
	if (useForwardDiff)
	{
		hT1 = T * 0.5;
	}

	// batch price using vectorized COS
	double vol_estimate = std::sqrt(_v0);

	// prices at T (for dCdK, d2CdK2)
	HestonCF cf_T(_kappa, _vbar, _sigma, _rho, _v0, _r, T);
	std::vector<double> prices_T = COSPricer::callPrices(
		_S0, allStrikes_T, _r, T, cf_T, _N, _L, vol_estimate);

	// prices for dCdT
	std::vector<double> prices_Tph1, prices_Tmh1, prices_Tph2, prices_Tmh2;

	if (useForwardDiff)
	{
		HestonCF cf_Tph1(_kappa, _vbar, _sigma, _rho, _v0, _r, T + hT1);
		prices_Tph1 = COSPricer::callPrices(
			_S0, centralStrikes, _r, T + hT1, cf_Tph1, _N, _L, vol_estimate);
	}
	else
	{
		// Richardson extrapolation for dCdT
		HestonCF cf_Tph1(_kappa, _vbar, _sigma, _rho, _v0, _r, T + hT1);
		HestonCF cf_Tmh1(_kappa, _vbar, _sigma, _rho, _v0, _r, T - hT1);
		HestonCF cf_Tph2(_kappa, _vbar, _sigma, _rho, _v0, _r, T + hT2);
		HestonCF cf_Tmh2(_kappa, _vbar, _sigma, _rho, _v0, _r, T - hT2);

		prices_Tph1 = COSPricer::callPrices(
			_S0, centralStrikes, _r, T + hT1, cf_Tph1, _N, _L, vol_estimate);
		prices_Tmh1 = COSPricer::callPrices(
			_S0, centralStrikes, _r, T - hT1, cf_Tmh1, _N, _L, vol_estimate);
		prices_Tph2 = COSPricer::callPrices(
			_S0, centralStrikes, _r, T + hT2, cf_Tph2, _N, _L, vol_estimate);
		prices_Tmh2 = COSPricer::callPrices(
			_S0, centralStrikes, _r, T - hT2, cf_Tmh2, _N, _L, vol_estimate);
	}

	// Dupire local vol for each spot
	for (size_t j = 0; j < M; ++j)
	{
		double K = spots[j];
		double h1 = h1_vec[j];
		double h2 = h2_vec[j];

		double C_Kmh1 = prices_T[5 * j + 0];
		double C_Kmh2 = prices_T[5 * j + 1];
		double C_K = prices_T[5 * j + 2];
		double C_Kph2 = prices_T[5 * j + 3];
		double C_Kph1 = prices_T[5 * j + 4];

		// dCdK with Richardson extrapolation
		double dCdK_h1 = (C_Kph1 - C_Kmh1) / (2.0 * h1);
		double dCdK_h2 = (C_Kph2 - C_Kmh2) / (2.0 * h2);
		double delta_K = (4.0 * dCdK_h2 - dCdK_h1) / 3.0;

		// d2CdK2 with Richardson extrapolation
		double d2CdK2_h1 = (C_Kph1 - 2.0 * C_K + C_Kmh1) / (h1 * h1);
		double d2CdK2_h2 = (C_Kph2 - 2.0 * C_K + C_Kmh2) / (h2 * h2);
		double gamma_K = (4.0 * d2CdK2_h2 - d2CdK2_h1) / 3.0;

		// dCdT
		double theta;
		if (useForwardDiff)
		{
			theta = (prices_Tph1[j] - C_K) / hT1;
		}
		else
		{
			double dCdT_h1 = (prices_Tph1[j] - prices_Tmh1[j]) / (2.0 * hT1);
			double dCdT_h2 = (prices_Tph2[j] - prices_Tmh2[j]) / (2.0 * hT2);
			theta = (4.0 * dCdT_h2 - dCdT_h1) / 3.0;
		}

		// Dupire formula
		double numerator = 2.0 * (theta + _r * K * delta_K);
		double denominator = K * K * gamma_K;

		double localVar;
		if (denominator < 1e-12)
		{
			localVar = _vbar; // fallback for deep OTM
		}
		else
		{
			localVar = numerator / denominator;
		}

		// sanity bounds
		localVar = std::max(localVar, 1e-6);
		localVar = std::min(localVar, 4.0); // 200% vol cap

		result[j] = std::sqrt(localVar);
	}

	return result;
}
