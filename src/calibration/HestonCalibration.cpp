#include <cmath>
#include <cppfm/calibration/HestonCalibration.h>
#include <cppfm/calibration/ImpliedVolSolver.h>
#include <cppfm/cos/COS.h>
#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cppfm/pricers/COSPricer.h>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

// ============================================================================
// File-local helpers
// ============================================================================

static std::vector<double> linspace(double lo, double hi, int n)
{
	std::vector<double> v(n);
	if (n == 1)
	{
		v[0] = 0.5 * (lo + hi);
		return v;
	}
	for (int i = 0; i < n; ++i)
		v[i] = lo + i * (hi - lo) / (n - 1);
	return v;
}

static void validateSlices(const std::vector<HestonSliceData> &slices, const char *caller)
{
	if (slices.empty())
		throw std::invalid_argument(std::string(caller) + ": no slices");
	for (auto &s : slices)
	{
		if (s.strikes.size() != s.marketIVs.size())
			throw std::invalid_argument(std::string(caller) + ": strikes/IVs size mismatch");
		if (s.strikes.empty())
			throw std::invalid_argument(std::string(caller) + ": empty slice");
	}
}

struct PrecomputedMarket
{
	std::vector<std::vector<double>> prices;
	std::vector<std::vector<double>> vegas;
};

// vega-weighted price residual: (C_model - C_mkt) / vega ≈ sigma_model - sigma_mkt
// avoids IV inversion inside LM → smooth Jacobian for all 5 params
static PrecomputedMarket precomputeMarketData(const std::vector<HestonSliceData> &slices,
											  double S0, double r)
{
	PrecomputedMarket mkt;
	mkt.prices.resize(slices.size());
	mkt.vegas.resize(slices.size());
	for (size_t s = 0; s < slices.size(); ++s)
	{
		auto &slice = slices[s];
		mkt.prices[s].resize(slice.strikes.size());
		mkt.vegas[s].resize(slice.strikes.size());
		for (size_t i = 0; i < slice.strikes.size(); ++i)
		{
			double iv = slice.marketIVs[i];
			mkt.prices[s][i] = BlackScholesFormulas::callPrice(S0, slice.strikes[i], r, iv, slice.T);
			mkt.vegas[s][i] = BlackScholesFormulas::vega(S0, slice.strikes[i], r, iv, slice.T);
			// floor vega so deep OTM doesn't blow up
			if (mkt.vegas[s][i] < 1e-6)
				mkt.vegas[s][i] = 1e-6;
		}
	}
	return mkt;
}

struct HestonGrid
{
	std::vector<double> v0, kappa, vbar, sig, rho;
};

// grid ranges biased for equity
static HestonGrid makeDefaultGrid(int gridPoints)
{
	HestonGrid g;
	g.v0 = linspace(0.01, 0.15, gridPoints);
	g.kappa = linspace(0.5, 5.0, gridPoints);
	g.vbar = linspace(0.01, 0.15, gridPoints);
	g.sig = linspace(0.1, 1.0, gridPoints);
	g.rho = linspace(-0.9, -0.1, gridPoints);
	return g;
}

// ============================================================================
// HestonParams
// ============================================================================

HestonParams::HestonParams(double v0, double kappa, double vbar, double sigma_v, double rho)
	: v0(v0), kappa(kappa), vbar(vbar), sigma_v(sigma_v), rho(rho) {}

HestonParams::HestonParams(const std::vector<double> &v)
{
	fromVector(v);
}

std::vector<double> HestonParams::toVector() const
{
	return {v0, kappa, vbar, sigma_v, rho};
}

void HestonParams::fromVector(const std::vector<double> &v)
{
	if (v.size() != 5)
		throw std::invalid_argument("HestonParams: need 5 values");
	v0 = v[0];
	kappa = v[1];
	vbar = v[2];
	sigma_v = v[3];
	rho = v[4];
}

std::vector<double> HestonParams::lowerBounds()
{
	// v0, kappa, vbar, sigma_v, rho
	return {1e-4, 1e-3, 1e-4, 1e-2, -0.999};
}

std::vector<double> HestonParams::upperBounds()
{
	return {1.0, 20.0, 1.0, 5.0, 0.999};
}

bool HestonParams::satisfiesFellerCondition() const
{
	return 2.0 * kappa * vbar >= sigma_v * sigma_v;
}

// ============================================================================
// calibrateHeston
// ============================================================================

HestonCalibrationResult calibrateHeston(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	const HestonParams &guess,
	const LMOptions &opts)
{
	validateSlices(slices, "calibrateHeston");

	size_t totalResiduals = 0;
	for (auto &s : slices)
		totalResiduals += s.strikes.size();

	auto mkt = precomputeMarketData(slices, S0, r);

	auto residualFn = [&](const std::vector<double> &x) -> std::vector<double>
	{
		HestonParams p(x);
		std::vector<double> res;
		res.reserve(totalResiduals);

		for (size_t s = 0; s < slices.size(); ++s)
		{
			auto &slice = slices[s];
			HestonCF cf(p.kappa, p.vbar, p.sigma_v, p.rho, p.v0, r, slice.T);
			auto chfFunc = [&cf](double u)
			{ return cf(u); };
			auto cum = cf.cumulants();

			auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
												chfFunc, 256, 10.0,
												cum.c1, cum.c2, cum.c4);

			for (size_t i = 0; i < slice.strikes.size(); ++i)
				res.push_back((prices[i] - mkt.prices[s][i]) / mkt.vegas[s][i]);
		}
		return res;
	};

	LevenbergMarquardt lm;
	auto lmResult = lm.solve(residualFn, guess.toVector(),
							 HestonParams::lowerBounds(),
							 HestonParams::upperBounds(),
							 nullptr, opts);

	HestonCalibrationResult result;
	result.params.fromVector(lmResult.params);
	result.iterations = lmResult.iterations;
	result.converged = lmResult.converged;
	result.message = lmResult.message;

	// RMSE in IV units for reporting (safe — final params produce good prices)
	HestonParams finalP(lmResult.params);
	double sumSq = 0;
	for (auto &slice : slices)
	{
		HestonCF cf(finalP.kappa, finalP.vbar, finalP.sigma_v, finalP.rho, finalP.v0, r, slice.T);
		auto chfFunc = [&cf](double u)
		{ return cf(u); };
		auto cum = cf.cumulants();
		auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
											chfFunc, 256, 10.0,
											cum.c1, cum.c2, cum.c4);

		double sliceSumSq = 0;
		for (size_t i = 0; i < slice.strikes.size(); ++i)
		{
			double modelIV = ImpliedVolSolver::solve(prices[i], S0, slice.strikes[i],
													 r, slice.T, true);
			double diff = modelIV - slice.marketIVs[i];
			sliceSumSq += diff * diff;
			sumSq += diff * diff;
		}
		result.sliceRmse.push_back(std::sqrt(sliceSumSq / slice.strikes.size()));
	}
	result.rmse = std::sqrt(sumSq / totalResiduals);

	return result;
}

// ============================================================================
// hestonGridSearch
// ============================================================================

HestonParams hestonGridSearch(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	int gridPoints)
{
	validateSlices(slices, "hestonGridSearch");

	auto grid = makeDefaultGrid(gridPoints);
	auto mkt = precomputeMarketData(slices, S0, r);

	double bestSSE = std::numeric_limits<double>::max();
	HestonParams best;
	// brute force grid search
	for (double v0 : grid.v0)
		for (double kappa : grid.kappa)
			for (double vbar : grid.vbar)
				for (double sig : grid.sig)
					for (double rho : grid.rho)
					{
						double sse = 0;

						for (size_t s = 0; s < slices.size(); ++s)
						{
							auto &slice = slices[s];
							HestonCF cf(kappa, vbar, sig, rho, v0, r, slice.T);
							auto chfFunc = [&cf](double u)
							{ return cf(u); };
							auto cum = cf.cumulants();

							auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
																chfFunc, 96, 10.0,
																cum.c1, cum.c2, cum.c4);

							for (size_t i = 0; i < slice.strikes.size(); ++i)
							{
								double diff = (prices[i] - mkt.prices[s][i]) / mkt.vegas[s][i];
								sse += diff * diff;
							}
						}

						if (sse < bestSSE)
						{
							bestSSE = sse;
							best = HestonParams(v0, kappa, vbar, sig, rho);
						}
					}

	return best;
}

// parallel grid search
HestonParams hestonGridSearchParallel(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	int gridPoints)
{
	validateSlices(slices, "hestonGridSearchParallel");

	auto grid = makeDefaultGrid(gridPoints);
	auto mkt = precomputeMarketData(slices, S0, r);

	size_t n0 = grid.v0.size(), n1 = grid.kappa.size(), n2 = grid.vbar.size(),
		   n3 = grid.sig.size(), n4 = grid.rho.size();
	size_t total = n0 * n1 * n2 * n3 * n4;

	unsigned nThreads = std::thread::hardware_concurrency();
	if (nThreads == 0)
		nThreads = 4;

	struct LocalBest
	{
		double sse = std::numeric_limits<double>::max();
		HestonParams params;
	};
	std::vector<LocalBest> threadBests(nThreads);
	std::vector<std::thread> threads;

	auto worker = [&](unsigned tid)
	{
		auto &local = threadBests[tid];

		for (size_t idx = tid; idx < total; idx += nThreads)
		{
			// decode flat index
			size_t rem = idx;
			size_t i4 = rem % n4;
			rem /= n4;
			size_t i3 = rem % n3;
			rem /= n3;
			size_t i2 = rem % n2;
			rem /= n2;
			size_t i1 = rem % n1;
			rem /= n1;
			size_t i0 = rem;

			double v0 = grid.v0[i0];
			double kappa = grid.kappa[i1];
			double vbar = grid.vbar[i2];
			double sig = grid.sig[i3];
			double rho = grid.rho[i4];

			double sse = 0;

			for (size_t s = 0; s < slices.size(); ++s)
			{
				auto &slice = slices[s];
				HestonCF cf(kappa, vbar, sig, rho, v0, r, slice.T);
				auto chfFunc = [&cf](double u)
				{ return cf(u); };
				auto cum = cf.cumulants();
				auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T, chfFunc, 96, 10.0,
													cum.c1, cum.c2, cum.c4);

				for (size_t i = 0; i < slice.strikes.size(); ++i)
				{
					double diff = (prices[i] - mkt.prices[s][i]) / mkt.vegas[s][i];
					sse += diff * diff;
				}
			}

			if (sse < local.sse)
			{
				local.sse = sse;
				local.params = HestonParams(v0, kappa, vbar, sig, rho);
			}
		}
	};
	for (unsigned t = 0; t < nThreads; ++t)
	{
		threads.emplace_back(worker, t);
	}
	for (auto &t : threads)
	{
		t.join();
	}
	// reduce thread-local bests
	double bestSSE = std::numeric_limits<double>::max();
	HestonParams best;
	for (auto &lb : threadBests)
	{
		if (lb.sse < bestSSE)
		{
			bestSSE = lb.sse;
			best = lb.params;
		}
	}
	return best;
}
