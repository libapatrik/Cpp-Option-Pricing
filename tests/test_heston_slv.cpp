#include "../PathSimulator2D.h"
#include "../VolatilitySurface.h"
#include "../BlackScholesFormulas.h"
#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <vector>

class HestonSLVTest : public ::testing::Test {
protected:
  FlatDiscountCurve discountCurve{0.05};
  HestonModel model{100.0, discountCurve, 0.04, 2.0, 0.04, 0.3, -0.7};

  std::vector<double> times = {0.0, 0.25, 0.5, 0.75, 1.0};

  std::vector<double> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
  std::vector<double> maturities = {0.25, 0.5, 0.75, 1.0, 1.5};
  std::vector<std::vector<double>> volatilities = {
      {0.25, 0.23, 0.20, 0.22, 0.24},
      {0.24, 0.22, 0.19, 0.21, 0.23},
      {0.23, 0.21, 0.18, 0.20, 0.22},
      {0.22, 0.20, 0.17, 0.19, 0.21},
      {0.21, 0.19, 0.16, 0.18, 0.20}};

  size_t numPaths = 1000;
  size_t numBins = 20;
  size_t seed = 42;
};


TEST_F(HestonSLVTest, BasicConstruction) {
  VolatilitySurface volSurface(strikes, maturities, volatilities, discountCurve);
  HestonSLVPathSimulator2D simulator(model, volSurface, times, numPaths, numBins, seed);
  SUCCEED();
}

TEST_F(HestonSLVTest, SimulateAllPaths) {
  VolatilitySurface volSurface(strikes, maturities, volatilities, discountCurve);

  size_t localNumPaths = 500;
  HestonSLVPathSimulator2D simulator(model, volSurface, times, localNumPaths);

  auto terminalValues = simulator.simulateAllPaths();

  EXPECT_EQ(terminalValues.size(), localNumPaths);

  for (const auto &[S_T, V_T] : terminalValues) {
    EXPECT_GT(S_T, 0.0) << "Terminal spot must be positive";
    EXPECT_GE(V_T, 0.0) << "Terminal variance must be non-negative";
  }

  double sumS = 0.0;
  for (const auto &[S_T, V_T] : terminalValues) {
    sumS += S_T;
  }
  double meanS = sumS / localNumPaths;
  double T = times.back();
  double expectedMean = 100.0 * std::exp(0.05 * T);

  EXPECT_NEAR(meanS, expectedMean, expectedMean * 0.15)
      << "Mean terminal spot should be close to forward price";
}


// ============================================================================
// Stoep Paper Comparison - Figures 4.1-4.2
// ============================================================================

class StoepComparisonTest : public ::testing::Test {
protected:
    // Stoep paper Table 2 parameters
    double S0 = 1.0;
    double r = 0.0;
    double v0 = 0.04;
    double kappa = 0.5;
    double vbar = 0.04;
    double sigma_v = 1.0;
    double rho = -0.9;

    FlatDiscountCurve discountCurve{r};

    size_t numPaths = 20000;
    size_t numBins = 25;
    size_t seed = 42;
    size_t stepsPerYear = 50;

    std::vector<double> createTimeGrid(double T) {
        size_t numSteps = static_cast<size_t>(std::ceil(T * stepsPerYear));
        std::vector<double> times(numSteps + 1);
        for (size_t i = 0; i <= numSteps; ++i) {
            times[i] = T * static_cast<double>(i) / numSteps;
        }
        return times;
    }

    std::unique_ptr<VolatilitySurface> createVolSurface() {
        std::vector<double> K = {0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2};
        std::vector<double> T = {0.25, 0.5, 1.0, 2.0, 5.0};

        // Simple skewed surface
        std::vector<std::vector<double>> vols(T.size());
        for (size_t t = 0; t < T.size(); ++t) {
            vols[t].resize(K.size());
            for (size_t k = 0; k < K.size(); ++k) {
                double atmVol = std::sqrt(v0);
                double logK = std::log(K[k]);
                vols[t][k] = atmVol - 0.15 * logK + 0.02 / std::sqrt(T[t]) * logK * logK;
                vols[t][k] = std::max(0.05, std::min(0.80, vols[t][k]));
            }
        }
        return std::make_unique<VolatilitySurface>(K, T, vols, discountCurve);
    }

    double priceCallMC(const std::vector<std::pair<double, double>>& terminals, double K, double T) {
        double sum = 0.0;
        for (const auto& [S_T, V_T] : terminals) {
            sum += std::max(S_T - K, 0.0);
        }
        return discountCurve.discount(T) * sum / terminals.size();
    }

    double priceCallHeston(double K, double T, const HestonModel& model, size_t nPaths) {
        auto times = createTimeGrid(T);
        QEPathSimulator2D sim(times, model, seed);

        double sum = 0.0;
        for (size_t p = 0; p < nPaths; ++p) {
            auto [path, var] = sim.paths();
            sum += std::max(path.back() - K, 0.0);
        }
        return discountCurve.discount(T) * sum / nPaths;
    }
};


TEST_F(StoepComparisonTest, ImpliedVolTable) {
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "  STOEP PAPER - IMPLIED VOL COMPARISON\n";
    std::cout << "  Heston: kappa=" << kappa << " vbar=" << vbar
              << " sigma_v=" << sigma_v << " rho=" << rho << "\n";
    std::cout << "================================================================\n\n";

    auto volSurface = createVolSurface();
    HestonModel heston(S0, discountCurve, v0, kappa, vbar, sigma_v, rho);

    std::vector<double> testT = {0.5, 1.0, 2.0};
    std::vector<double> testK = {0.9, 0.95, 1.0, 1.05, 1.1};

    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(6) << "T"
              << std::setw(8) << "K"
              << std::setw(10) << "Market"
              << std::setw(10) << "SLV"
              << std::setw(10) << "Heston"
              << std::setw(12) << "SLV err"
              << std::setw(12) << "Hes err"
              << "\n";
    std::cout << std::string(68, '-') << "\n";

    for (double T : testT) {
        auto times = createTimeGrid(T);
        HestonSLVPathSimulator2D slvSim(heston, *volSurface, times, numPaths, numBins, seed);
        auto slvTerminals = slvSim.simulateAllPaths();

        for (double K : testK) {
            double mktVol = volSurface->impliedVolatility(K, T);

            double slvPrice = priceCallMC(slvTerminals, K, T);
            double slvVol = BlackScholesFormulas::impliedVolatility(S0, K, r, T, slvPrice, Option::Type::Call, mktVol);

            double hesPrice = priceCallHeston(K, T, heston, numPaths / 5);
            double hesVol = BlackScholesFormulas::impliedVolatility(S0, K, r, T, hesPrice, Option::Type::Call, mktVol);

            double slvErr = (slvVol - mktVol) * 10000;
            double hesErr = (hesVol - mktVol) * 10000;

            std::cout << std::setw(6) << T
                      << std::setw(8) << K
                      << std::setw(9) << mktVol * 100 << "%"
                      << std::setw(9) << slvVol * 100 << "%"
                      << std::setw(9) << hesVol * 100 << "%"
                      << std::setw(10) << std::showpos << slvErr << "bp"
                      << std::setw(10) << hesErr << "bp"
                      << std::noshowpos << "\n";
        }
        std::cout << "\n";
    }
}

// Verify the martingale property
TEST_F(StoepComparisonTest, TerminalStatistics) {
    auto volSurface = createVolSurface();
    HestonModel heston(S0, discountCurve, v0, kappa, vbar, sigma_v, rho);

    double T = 1.0;
    auto times = createTimeGrid(T);

    HestonSLVPathSimulator2D sim(heston, *volSurface, times, numPaths, numBins, seed);
    auto terminals = sim.simulateAllPaths();

    double sumS = 0.0, sumV = 0.0;
    for (const auto& [S, V] : terminals) {
        sumS += S;
        sumV += V;
    }

    double meanS = sumS / terminals.size();
    double meanV = sumV / terminals.size();

    std::cout << "\nTerminal stats (T=1):\n";
    std::cout << "  E[S_T] = " << meanS << " (fwd = " << S0 * std::exp(r * T) << ")\n";
    std::cout << "  E[V_T] = " << meanV << " (vbar = " << vbar << ")\n\n";

    EXPECT_NEAR(meanS, S0 * std::exp(r * T), 0.1);
}

