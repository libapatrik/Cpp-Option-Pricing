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
	if (slices.empty())
		throw std::invalid_argument("calibrateHeston: no slices");

	size_t totalResiduals = 0;
	for (auto &s : slices)
	{
		if (s.strikes.size() != s.marketIVs.size())
			throw std::invalid_argument("calibrateHeston: strikes/IVs size mismatch");
		if (s.strikes.empty())
			throw std::invalid_argument("calibrateHeston: empty slice");
		totalResiduals += s.strikes.size();
	}

	// precompute market prices and vegas from market IVs (constant during LM)
	// vega-weighted price residual: (C_model - C_mkt) / vega ≈ sigma_model - sigma_mkt
	// avoids IV inversion inside LM → smooth Jacobian for all 5 params
	std::vector<std::vector<double>> mktPrices(slices.size());
	std::vector<std::vector<double>> vegas(slices.size());
	for (size_t s = 0; s < slices.size(); ++s)
	{
		auto &slice = slices[s];
		mktPrices[s].resize(slice.strikes.size());
		vegas[s].resize(slice.strikes.size());
		for (size_t i = 0; i < slice.strikes.size(); ++i)
		{
			double iv = slice.marketIVs[i];
			mktPrices[s][i] = BlackScholesFormulas::callPrice(S0, slice.strikes[i], r, iv, slice.T);
			vegas[s][i] = BlackScholesFormulas::vega(S0, slice.strikes[i], r, iv, slice.T);
			// floor vega so deep OTM doesn't blow up
			if (vegas[s][i] < 1e-6)
				vegas[s][i] = 1e-6;
		}
	}

	auto residualFn = [&](const std::vector<double> &x) -> std::vector<double>
	{
		HestonParams p(x);
		std::vector<double> res;
		res.reserve(totalResiduals);

		// COS domain needs to cover the wider of initial and long-term vol
		double sigmaHint = std::sqrt(std::max(p.v0, p.vbar));

		for (size_t s = 0; s < slices.size(); ++s)
		{
			auto &slice = slices[s];
			HestonCF cf(p.kappa, p.vbar, p.sigma_v, p.rho, p.v0, r, slice.T);
			auto chfFunc = [&cf](double u)
			{ return cf(u); };

			auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
												chfFunc, 256, 10.0, sigmaHint);

			for (size_t i = 0; i < slice.strikes.size(); ++i)
				res.push_back((prices[i] - mktPrices[s][i]) / vegas[s][i]);
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
		double sigmaHint = std::sqrt(std::max(finalP.v0, finalP.vbar));
		HestonCF cf(finalP.kappa, finalP.vbar, finalP.sigma_v, finalP.rho, finalP.v0, r, slice.T);
		auto chfFunc = [&cf](double u)
		{ return cf(u); };
		auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
											chfFunc, 256, 10.0, sigmaHint);

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
	if (slices.empty())
		throw std::invalid_argument("hestonGridSearch: no slices");
	for (auto &s : slices)
	{
		if (s.strikes.size() != s.marketIVs.size())
			throw std::invalid_argument("hestonGridSearch: strikes/IVs size mismatch");
		if (s.strikes.empty())
			throw std::invalid_argument("hestonGridSearch: empty slice");
	}

	auto linspace = [](double lo, double hi, int n) -> std::vector<double>
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
	};

	// grid ranges biased for equity
	auto v0Grid = linspace(0.01, 0.15, gridPoints);
	auto kappaGrid = linspace(0.5, 5.0, gridPoints);
	auto vbarGrid = linspace(0.01, 0.15, gridPoints);
	auto sigGrid = linspace(0.1, 1.0, gridPoints);
	auto rhoGrid = linspace(-0.9, -0.1, gridPoints);

	// precompute market prices and vegas (same as calibrateHeston)
	std::vector<std::vector<double>> mktPrices(slices.size());
	std::vector<std::vector<double>> vegas(slices.size());
	for (size_t s = 0; s < slices.size(); ++s)
	{
		auto &slice = slices[s];
		mktPrices[s].resize(slice.strikes.size());
		vegas[s].resize(slice.strikes.size());
		for (size_t i = 0; i < slice.strikes.size(); ++i)
		{
			double iv = slice.marketIVs[i];
			mktPrices[s][i] = BlackScholesFormulas::callPrice(S0, slice.strikes[i], r, iv, slice.T);
			vegas[s][i] = BlackScholesFormulas::vega(S0, slice.strikes[i], r, iv, slice.T);
			if (vegas[s][i] < 1e-6)
				vegas[s][i] = 1e-6;
		}
	}

	double bestSSE = std::numeric_limits<double>::max();
	HestonParams best;
	// brute force grid search
	// TODO: Parallelize it?
	for (double v0 : v0Grid)
		for (double kappa : kappaGrid)
			for (double vbar : vbarGrid)
				for (double sig : sigGrid)
					for (double rho : rhoGrid)
					{
						double sse = 0;
						double sigmaHint = std::sqrt(std::max(v0, vbar));

						for (size_t s = 0; s < slices.size(); ++s)
						{
							auto &slice = slices[s];
							HestonCF cf(kappa, vbar, sig, rho, v0, r, slice.T);
							auto chfFunc = [&cf](double u)
							{ return cf(u); };

							auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T,
																chfFunc, 128, 10.0, sigmaHint);

							for (size_t i = 0; i < slice.strikes.size(); ++i)
							{
								double diff = (prices[i] - mktPrices[s][i]) / vegas[s][i];
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

// Parallel version:
// NOTE: to be called by the use either prove own guess or grid-search first for params
HestonParams hestonGridSearchParallel(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	int gridPoints)
{
	if (slices.empty())
		throw std::invalid_argument("hestonGridSearch: no slices");
	for (auto &s : slices)
	{
		if (s.strikes.size() != s.marketIVs.size())
			throw std::invalid_argument("hestonGridSearch: strikes/IVs size mismatch");
		if (s.strikes.empty())
			throw std::invalid_argument("hestonGridSearch: empty slice");
	}

	auto linspace = [](double lo, double hi, int n) -> std::vector<double>
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
	};

	// grid ranges biased for equity
	auto v0Grid = linspace(0.01, 0.15, gridPoints);
	auto kappaGrid = linspace(0.5, 5.0, gridPoints);
	auto vbarGrid = linspace(0.01, 0.15, gridPoints);
	auto sigGrid = linspace(0.1, 1.0, gridPoints);
	auto rhoGrid = linspace(-0.9, -0.1, gridPoints);

	// precompute market prices and vegas (same as calibrateHeston)
	std::vector<std::vector<double>> mktPrices(slices.size());
	std::vector<std::vector<double>> vegas(slices.size());
	for (size_t s = 0; s < slices.size(); ++s)
	{
		auto &slice = slices[s];
		mktPrices[s].resize(slice.strikes.size());
		vegas[s].resize(slice.strikes.size());
		for (size_t i = 0; i < slice.strikes.size(); ++i)
		{
			double iv = slice.marketIVs[i];
			mktPrices[s][i] = BlackScholesFormulas::callPrice(S0, slice.strikes[i], r, iv, slice.T);
			vegas[s][i] = BlackScholesFormulas::vega(S0, slice.strikes[i], r, iv, slice.T);
			if (vegas[s][i] < 1e-6)
				vegas[s][i] = 1e-6;
		}
	}

	size_t n0 = v0Grid.size(), n1 = kappaGrid.size(), n2 = vbarGrid.size(), n3 = sigGrid.size(), n4 = rhoGrid.size();
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
			// decod flat index
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

			double v0 = v0Grid[i0];
			double kappa = kappaGrid[i1];
			double vbar = vbarGrid[i2];
			double sig = sigGrid[i3];
			double rho = rhoGrid[i4];

			double sse = 0;
			double sigmaHint = std::sqrt(std::max(v0, vbar));

			for (size_t s = 0; s < slices.size(); ++s)
			{
				auto &slice = slices[s];
				HestonCF cf(kappa, vbar, sig, rho, v0, r, slice.T);
				auto chfFunc = [&cf](double u)
				{ return cf(u); };

				auto prices = COSPricer::callPrices(S0, slice.strikes, r, slice.T, chfFunc, 128, 10.0, sigmaHint);

				for (size_t i = 0; i < slice.strikes.size(); ++i)
				{
					double diff = (prices[i] - mktPrices[s][i]) / vegas[s][i];
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
};