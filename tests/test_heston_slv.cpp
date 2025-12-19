//
// Created by Patrik  Liba on 18/12/2025.
//
//  tests/test_heston_slv.cpp
#include "../PathSimulator2D.h"
#include "../VolatilitySurface.h"
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

class HestonSLVTest : public ::testing::Test {
protected:
  // Create Heston model with realistic parameters
  // HestonModel(spot, discountCurve, v0, kappa, vbar, sigma_v, rho)
  FlatDiscountCurve discountCurve{0.05};
  HestonModel model{100.0, discountCurve, 0.04, 2.0, 0.04, 0.3, -0.7};

  // Time grid for simulation
  std::vector<double> times = {0.0, 0.25, 0.5, 0.75, 1.0};

  // Create volatility surface
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
  VolatilitySurface volSurface(strikes, maturities, volatilities,
                               discountCurve);

  HestonSLVPathSimulator2D simulator(model, volSurface, times, numPaths,
                                     numBins, seed);

  // Basic construction test - just ensure no crash
  SUCCEED();
}

TEST_F(HestonSLVTest, SimulateAllPaths) {
  VolatilitySurface volSurface(strikes, maturities, volatilities,
                               discountCurve);

  size_t localNumPaths = 500;
  HestonSLVPathSimulator2D simulator(model, volSurface, times, localNumPaths);

  // Simulate all paths
  auto terminalValues = simulator.simulateAllPaths();

  // Check we got the expected number of paths
  EXPECT_EQ(terminalValues.size(), localNumPaths);

  // Check all terminal values are positive (S > 0, V > 0)
  for (const auto &[S_T, V_T] : terminalValues) {
    EXPECT_GT(S_T, 0.0) << "Terminal spot must be positive";
    EXPECT_GE(V_T, 0.0) << "Terminal variance must be non-negative";
  }

  // Compute mean terminal spot - should be close to S0 * exp(rT) for
  // risk-neutral pricing
  double sumS = 0.0;
  for (const auto &[S_T, V_T] : terminalValues) {
    sumS += S_T;
  }
  double meanS = sumS / localNumPaths;
  double T = times.back();
  double expectedMean = 100.0 * std::exp(0.05 * T); // S0 * exp(rT)

  // Allow some Monte Carlo error (within 10% for small sample)
  EXPECT_NEAR(meanS, expectedMean, expectedMean * 0.15)
      << "Mean terminal spot should be close to forward price";
}
