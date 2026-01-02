//
// Created by Patrik  Liba on 27/11/2025.
//

#ifndef CPPFM_PATHSIMULATOR2D_H
#define CPPFM_PATHSIMULATOR2D_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Model.h"
#include "VolatilitySurface.h"
#include <random>
#include <vector>

/// GOAL: Andersen -> Stoep -> Hammouda QMC
/**
 *  NOTES:
 *
 *  Diagram of discretization schemes: (Eq = equation from the Andersen paper)
 *  X:   Euler    Milstein     BK        Eq:33            Eq:33
 *        |          |         |          |                 |
 *  V:   Euler    Milstein     BK     TG (Eq:13)     QE (Eq:23-26)
 *
 *  TODO:
 *  1. COS-method: given ChF recover CDF + PDF
 *  2. Newton method to invert CDF to get integrated variance sample
 *  NOTE: MC-Optimization: How much better can we make the MC simulation?
 *        Control variates (CV) vs. Antithetic sampling (AS) vs. Quasi-MC
 *        CV +AS, CV + QMC works, but AS + QMC doesn't work well
 *        Multilevel MC + CV + AS (per level)
 *        3. DONE: Antithetic sampling
 *        4. TODO: Control variates choice on top of antithetic sampling
 *
 *
 *  MAKE:
 *  (1) Could pre-compute and cache the repeated parameters
 *  (2) Adaptive γ [Dufresne 1998]
 *  (3)
 */

/**
 * Monte Carlo Simulation Methods for Heston Model
 * Current Implementation Status:
 *
 *   ├─ Scheme 1: Euler-Euler [DONE]
 *   │  ├─ V: Euler with full truncation (eq. 6-7)
 *   │  └─ X: Euler (eq. 6)
 *   │
 *   ├─ Scheme 2: TG (Andersen) [DONE]
 *   │  ├─ V: Truncated Gaussian (eq. 13)
 *   │  └─ X: Eqaution 33
 *   │
 *   ├─ Scheme 3: QE (Andersen) [DONE]
 *   │  ├─ V: Quadratic-Exponential (eq. 23/26)
 *   │  └─ X: Equation 33
 *   │
 *   └─ BK Family (Broadie-Kaya - all use Equation 11 for X):
 *      │
 *      ├─ BKExact [DONE] - Original Broadie-Kaya
 *      │  ├─ V: Exact CIR (non-central χ²)
 *      │  ├─ ∫V: Exact via ChF + COS/Newton
 *      │  └─ X: Equation 11 (exact)
 *      │
 *      ├─ BKTG [DONE] - Hybrid scheme
 *      │  ├─ V: Truncated Gaussian (eq. 13)
 *      │  ├─ ∫V: Approximation (32) - γ₁V(t) + γ₂V(t+Δ)
 *      │  └─ X: Equation 11
 *      │
 *      └─ BKQE [DONE] - Hybrid scheme
 *         ├─ V: Quadratic-Exponential (eq. 23/26)
 *         ├─ ∫V: Approximation (32) - γ₁V(t) + γ₂V(t+Δ)
 *         └─ X: Equation 11

 *
 * Future Work:
 *   - IM-IJK (Kahl-Jäckel): Implicit Milstein for V, IJK for X
 *   - Multilevel Monte Carlo (MLMC) with any BK variant
 *   - Quasi-Monte Carlo (QMC) integration
 */

// Abstract class for different discretization schemes for 2D models (e.g.,
// Heston)
class PathSimulator2D
{
public:
  PathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model,
                  size_t randomSeed);
  virtual ~PathSimulator2D();
  virtual std::pair<double, double>
  nextStep(size_t timeIndex, double assetPrice, double variance) const = 0;

  // Generate full paths for both S and V
  std::pair<std::vector<double>, std::vector<double>> paths() const;

protected:
  bool timeStepsSanityCheck() const;

  // Generate correlated normal random variables with correlation from model
  // NOTE: Uses antithetic sampling internally for variance reduction
  std::pair<double, double> generateCorrelatedNormals() const;

  double generateUniform() const;
  double generateStandardNormalSLV() const;
  double generateStandardNormal() const;

  // Thread-local RNGs for parallel simulation
  // ? faster with PCG or xoroshiro123++
  mutable std::vector<std::default_random_engine> _threadRNGs;
  // initialize thread-local RNGs (call once before parallel region)
  void initThreadLocalRNGs() const;
  // Generate random using thread's own RNG
  double generateStandardNormalParallel(int threadId) const;

  /**
   * Correlation-preserving discretization for X (Andersen Eq. 33)
   * Common to both TG and QE schemes
   *
   * Notation: X(t) = ln(S(t)) is the log-price
   *
   * @param X_t Current log-price X(t) = ln(S(t))
   * @param V_t Current variance V(t)
   * @param V_next Next variance V(t+Δ)
   * @param dt Time step Δt
   * @param Z Independent standard normal for X
   * @param gamma1 γ₁ weight for V(t)
   * @param gamma2 γ₂ weight for V(t+Δ)
   * @return X(t+Δ) = ln(S(t+Δ))
   */
  double stepLogPriceEq33(double X_t, double V_t, double V_next, double dt,
                          double Z, double gamma1, double gamma2) const;

  /**
   * Shared variance discretization methods (DRY principle)
   * These are used by multiple schemes to avoid code duplication
   */

  /**
   * Truncated Gaussian (TG) Scheme for Variance (Andersen Eq. 13)
   * Used by: TGPathSimulator2D, BKTGPathSimulator2D
   *
   * @param V_t Current variance V(t)
   * @param dt Time step Δt
   * @param Z_V Standard normal for variance
   * @return V(t+Δ)
   */
  double stepVarianceTG(double V_t, double dt, double Z_V) const;

  /**
   * Quadratic-Exponential (QE) Scheme for Variance (Andersen Eq. 23/26)
   * Used by: QEPathSimulator2D, BKQEPathSimulator2D
   *
   * @param V_t Current variance V(t)
   * @param dt Time step Δt
   * @param Z_V Standard normal for variance
   * @param psi_c Threshold for switching between quadratic/exponential
   * @return V(t+Δ)
   */
  double stepVarianceQE(double V_t, double dt, double Z_V, double psi_c) const;

  const std::vector<double> &_timeSteps;
  const Model2D
      *_modelPtr; // Pointer to base class Model2D - for polymorphic 2D models
  size_t _randomSeed;
  mutable std::default_random_engine _randomEngine;

  // Antithetic sampling support (PATH-LEVEL for proper variance reduction)
  // For N paths, generate N/2 with fresh randoms, then N/2 with negated randoms
  mutable bool _antitheticMode; // true = using cached negated values
  mutable std::vector<double>
      _randomCache;           // Cache ALL randoms from previous path
  mutable size_t _cacheIndex; // Current position in cache during replay
};

/** EULER as a pair for (X, V)
 * Full Truncation Scheme for Heston Model (Equations 6-7)
 *
 * Notation: X(t) = ln(S(t)) is the log-price
 *
 * X̂(t + Δ) = X̂(t) - (1/2)V̂(t)⁺Δ + √(V̂(t)⁺)·Z_X·√Δ
 * V̂(t + Δ) = V̂(t) + κ(θ - V̂(t)⁺)Δ + ε·√(V̂(t)⁺)·Z_V·√Δ
 *
 * where x⁺ = max(x, 0), and Z_X, Z_V are correlated standard normals with
 * correlation ρ.
 *
 * This scheme prevents negative variance by using V⁺ in all terms.
 * When V goes negative, it reverts back deterministically via drift term κθΔ.
 */
class EulerPathSimulator2D : public PathSimulator2D
{
public:
  EulerPathSimulator2D(const std::vector<double> &timeSteps,
                       const Model2D &model, size_t randomSeed);
  /**
   * Perform one time step using Full Truncation discretization
   * @param timeIndex Current time index in _timeSteps
   * @param assetPrice Current asset price S(t)
   * @param variance Current variance V(t)
   * @return Pair of (S(t+Δ), V(t+Δ))
   */
  std::pair<double, double> nextStep(size_t timeIndex, double assetPrice,
                                     double variance) const override;
};

class MilsteinPathSimulator2D : public PathSimulator2D
{
};

enum class NewtonMethod
{
  Original, // my implementation;
  Optimized // Optimized version - cached coefficients (as default)
};

// ============================================================================
// BK Scheme Base: Equation 11 (Exact SDE Representation)
// ============================================================================

/**
 * Abstract base class for all Broadie-Kaya (BK) schemes
 *
 * Implements Equation 11 from Andersen (2008):
 * ln X(t+Δ) = ln X(t) + r·Δt + (ρ/ε)(V(t+Δ) - V(t) - κθΔ)
 *                             + (κρ/ε - 1/2)∫V(s)ds
 *                             + √(1-ρ²)·√(∫V)·Z
 *
 * This is the EXACT representation. Subclasses differ in how they obtain:
 * 1. V(t+Δ): Exact sampling vs. discretization schemes
 * 2. ∫V: Exact sampling vs. approximation (Equation 32)
 */
class BKSchemeBase : public PathSimulator2D
{
public:
  BKSchemeBase(const std::vector<double> &timeSteps, const Model2D &model,
               size_t randomSeed);
  virtual ~BKSchemeBase() = default;

  /**
   * Common nextStep for all BK variants
   * Delegates to generateVarianceAndIntegral(), then applies Equation 11
   */
  std::pair<double, double> nextStep(size_t timeIndex, double assetPrice,
                                     double variance) const override final;

protected:
  /**
   * Broadie-Kaya Equation 11 (Exact SDE Representation)
   *
   * From Andersen (2008) Equation 11:
   * ln X(t+Δ) = ln X(t) + r·Δt + (ρ/ε)(V(t+Δ) - V(t) - κθΔ)
   *                             + (κρ/ε - 1/2)∫V(s)ds
   *                             + √(1-ρ²)·√(∫V)·Z
   *
   * Shared by ALL BK variants (exact and approximate)
   */
  double stepLogPriceEq11(double X_t, double V_t, double V_next,
                          double integratedV, double dt, double Z) const;

  /**
   * Pure virtual: Generate V(t+Δ) and ∫V
   *
   * Different strategies:
   * - BKExact: Exact sampling of both
   * - BKApproximateScheme: V via TG/QE, ∫V via approximation (32)
   */
  virtual std::pair<double, double>
  generateIntegratedVariance(double V_t, double dt) const = 0;
};

// ============================================================================
// BK Approximate Scheme: Equations 32 & 33
// ============================================================================

/**
 * Abstract base for BK schemes using approximation (32) for integrated
 * variance: ∫[t,t+Δ] V(u)du ≈ γ₁V(t) + γ₂V(t+Δ)
 *
 * This leads to efficient computation via Equation 33 with K coefficients.
 *
 * Subclasses (BKTG, BKQE) differ ONLY in how they generate V(t+Δ).
 */
class BKApproximateScheme : public BKSchemeBase
{
public:
  BKApproximateScheme(const std::vector<double> &timeSteps,
                      const Model2D &model, size_t randomSeed,
                      double gamma1 = 0.5, double gamma2 = 0.5);

protected:
  /**
   * Implements approximation (32): ∫V ≈ γ₁V(t) + γ₂V(t+Δ)
   * Uses trapezoidal rule by default (γ₁ = γ₂ = 0.5)
   */
  std::pair<double, double>
  generateIntegratedVariance(double V_t, double dt) const override final;
  /**
   * Pure virtual: Generate V(t+Δ)
   * Each subclass implements its own variance discretization
   */
  virtual double generateNextVariance(double V_t, double dt) const = 0;

  double _gamma1; // γ₁: weight for V(t) in approximation (32)
  double _gamma2; // γ₂: weight for V(t+Δ) in approximation (32)
};

// ============================================================================
// BKTG: TG Variance + Approximation (32)
// ============================================================================

/**
 * BKTG: Truncated Gaussian variance with approximation (32)
 *
 * Differs from BKQE ONLY in variance generation (TG vs QE)
 */
class BKTGPathSimulator2D : public BKApproximateScheme
{
public:
  BKTGPathSimulator2D(const std::vector<double> &timeSteps,
                      const Model2D &model, size_t randomSeed,
                      double gamma1 = 0.5, double gamma2 = 0.5);

protected:
  double generateNextVariance(double V_t, double dt) const override;
};

// ============================================================================
// BKQE: QE Variance + Approximation (32)
// ============================================================================

/**
 * BKQE: Quadratic-Exponential variance with approximation (32)
 *
 * Differs from BKTG ONLY in variance generation (QE vs TG)
 */
class BKQEPathSimulator2D : public BKApproximateScheme
{
public:
  BKQEPathSimulator2D(const std::vector<double> &timeSteps,
                      const Model2D &model, size_t randomSeed,
                      double psi_c = 1.5, double gamma1 = 0.5,
                      double gamma2 = 0.5);

protected:
  double generateNextVariance(double V_t, double dt) const override;

private:
  double _psi_c; // Threshold for QE switching rule
};

// ============================================================================
// BKExact: Exact CIR + Exact ∫V (Original Broadie-Kaya)
// ============================================================================

/**
 * BKExact: Original Broadie-Kaya with exact sampling
 *
 * Bypasses approximation (32) entirely:
 * - V(t+Δ): Exact CIR sampling via non-central χ^2
 * - ∫V: Exact sampling via characteristic function + COS/Newton inversion
 *
 * Most accurate but slowest method
 */
class BKExactPathSimulator2D : public BKSchemeBase
{
public:
  BKExactPathSimulator2D(const std::vector<double> &timeSteps,
                         const Model2D &model, size_t randomSeed,
                         NewtonMethod newtonMethod = NewtonMethod::Optimized);

protected:
  std::pair<double, double>
  generateIntegratedVariance(double V_t, double dt) const override;

private:
  NewtonMethod _newtonMethod;
};

// For X: (Equation 33) For V: Truncated Gaussian (TG) (Equation 13)
class TGPathSimulator2D : public PathSimulator2D
{
public:
  TGPathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model,
                    size_t randomSeed, double gamma1 = 0.5,
                    double gamma2 = 0.5);
  virtual ~TGPathSimulator2D() override;

  /**
   * Perform one time step using TG scheme for V and Equation 33 for X
   * @param timeIndex Current time index in _timeSteps
   * @param assetPrice Current asset price S(t)
   * @param variance Current variance V(t)
   * @return Pair of (S(t+Δ), V(t+Δ))
   */
  std::pair<double, double> nextStep(size_t timeIndex, double assetPrice,
                                     double variance) const override;

protected:
  double _gamma1; // γ₁: weight for V(t) in integral approximation
  double _gamma2; // γ₂: weight for V(t+Δt) in integral approximation
};

// For X: (Equation 33) For V: Quadratic-Exponential (QE) (Equation 23/26)
class QEPathSimulator2D : public PathSimulator2D
{
public:
  QEPathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model,
                    size_t randomSeed,
                    double psi_c = 1.5, // for the switching rule
                    double gamma1 = 0.5, double gamma2 = 0.5);
  virtual ~QEPathSimulator2D() override;

  /**
   * Perform one time step using QE scheme for V and Equation 33 for X
   * @param timeIndex Current time index in _timeSteps
   * @param assetPrice Current asset price S(t)
   * @param variance Current variance V(t)
   * @return Pair of (S(t+Δ), V(t+Δ))
   */
  std::pair<double, double> nextStep(size_t timeIndex, double assetPrice,
                                     double variance) const override;

protected:
  double _psi_c;  // ψ_c: Threshold parameter for switching between
                  // quadratic/exponential
  double _gamma1; // γ₁: weight for V(t) in integral approximation
  double _gamma2; // γ₂: weight for V(t+Δt) in integral approximation
};

// ============================================================================
// HestonSLVPathSimulator2D
// ============================================================================

/**
 * Heston Stochastic Local Volatility (SLV) Simulator
 * Reference: van der Stoep et al. (2013) "The Heston SLV Model"
 *
 * Inherits from PathSimulator2D to reuse:
 * - _modelPtr (Model2D*) - cast to HestonModel* when needed
 * - _randomEngine - for random generation
 * - _timeSteps (reference) - caller must ensure lifetime
 * - generateStandardNormalSLV(), generateUniform() - random utilities
 * - stepVarianceQE() - QE variance discretization
 *
 * SLV-specific additions:
 * - _volSurfacePtr - for Dupire local volatility σ_LV(S,t)
 * - Binning algorithm (Algorithm 1) for E[V|S]
 * - Leverage function L2(t,S) = σ^2_LV / E[V|S]
 * - Batch simulation (all paths simultaneously)
 */
class HestonSLVPathSimulator2D : public PathSimulator2D
{
public:
  HestonSLVPathSimulator2D(
      const HestonModel &model, const VolatilitySurface &volSurface,
      const std::vector<double> &timeSteps, // Caller must ensure lifetime
      size_t numPaths,
      size_t numBins = 20, // Paper recommends 20 bins
      size_t randomSeed = 42);

  ~HestonSLVPathSimulator2D() override;

  HestonSLVPathSimulator2D(const HestonSLVPathSimulator2D &) = delete;
  HestonSLVPathSimulator2D &
  operator=(const HestonSLVPathSimulator2D &) = delete;

  /**
   * nextStep() not used - SLV uses batch simulation
   * @throws std::runtime_error always
   */
  std::pair<double, double> nextStep(size_t timeIndex, double assetPrice,
                                     double variance) const override;

  /**
   * Simulate ALL paths simultaneously (batch mode)
   * Implements Algorithm 1 + Section 3.3 (QE scheme)
   * @return vector of (S_T, V_T) pairs for all paths at terminal time
   */
  std::vector<std::pair<double, double>> simulateAllPaths() const;
  std::vector<std::pair<double, double>> simulateAllPathsParallel() const;
  std::vector<std::pair<double, double>> simulateAllPathsOptimized() const;

  /**
   * Simulate and return full path history
   * @return [path_index][time_index] -> (S, V)
   */
  std::vector<std::vector<std::pair<double, double>>> simulateAllPathsFull() const;

private:
  // =========================================================================
  // Binning Data Structure (Section 3.1)
  // =========================================================================
  struct BinData
  {
    std::vector<size_t> pathIndices; // Paths assigned to this bin (J_{i,k})
    double lowerBound;               // b_k
    double upperBound;               // b_{k+1}
    double midpoint;                 // (b_k + b_{k+1}) / 2 - for interpolation
    double conditionalExpectation;   // E[V | S in bin_k]
  };

  // =========================================================================
  // Binning Algorithm 1
  // =========================================================================

  /**
   * Compute bins and conditional expectations E[V|S] at current time step
   * @param spotValues Current spot values for all paths
   * @param varianceValues Current variance values for all paths
   * @return Vector of bins with computed conditional expectations
   */
  std::vector<BinData> computeBins(const std::vector<double> &spotValues, const std::vector<double> &varianceValues) const;

  std::vector<BinData> computeBinsVectorized(const std::vector<double> &spotValues, const std::vector<double> &varianceValues) const;

  /**
   * Find bin index for a given spot value (binary search)
   * @param spot Spot value to locate
   * @param bins Computed bins
   * @return Bin index k such that spot in [b_k, b_{k+1})
   */
  size_t findBinIndex(double spot, const std::vector<BinData> &bins) const;

  // =========================================================================
  // Leverage Function (Eq. 2.9)
  // =========================================================================

  /**
   * Compute leverage function squared: L2(t,S) = σ^2_LV(t,S) / E[V|S]
   * This is the core of the SLV model (Equation 2.9)
   * Uses linear interpolation for E[V|S] (Section 3.2)
   *
   * @param spot Current spot S
   * @param time Current time t
   * @param bins Pre-computed bins with conditional expectations
   * @return L2(t, S)
   */
  double leverageSquared(double spot, double time,
                         const std::vector<BinData> &bins) const;

  /**
   * Linear interpolation of E[V|S] between bin midpoints
   * @param spot Spot value for interpolation
   * @param bins Computed bins with conditional expectations
   * @return Interpolated E[V|S]
   */
  double
  interpolateConditionalExpectation(double spot,
                                    const std::vector<BinData> &bins) const;

  // =========================================================================
  // SLV-Specific Discretization (Section 3.3)
  // =========================================================================

  /**
   * Step variance using QE scheme with explicit uniform input
   * Wrapper around base class stepVarianceQE for SLV batch simulation
   * @param V_t Current variance
   * @param dt Time step
   * @param Z_V Standard normal for variance
   * @param U_V Uniform(0,1) for QE switching
   * @return V_{t+dt}
   */
  // ? Reuse this from before no?
  double stepVarianceQE_SLV(double V_t, double dt, double Z_V,
                            double U_V) const;

  /**
   * Step log-price using SLV scheme (Equation 3.18)
   * X_{i+1} = X_i + rΔt - ½L²V_iΔt
   *         + (ρ/γ)L(V_{i+1} - V_i - κθΔt + κV_iΔt)
   *         + √(1-ρ²)·√(L²V_iΔt)·Z
   *
   * @param X_t Current log-price ln(S)
   * @param V_t Current variance
   * @param V_next Next variance (already sampled)
   * @param leverageSq L²(t, S) from leverage function
   * @param dt Time step
   * @param Z Standard normal for X
   * @return X_{t+Δt}
   */
  double stepLogPriceSLV(double X_t, double V_t, double V_next,
                         double leverageSq, double dt, double Z) const;

  // =========================================================================
  // SLV-Specific Member Variables
  // =========================================================================
  //
  // INHERITED from PathSimulator2D (DO NOT REDECLARE):
  // - _modelPtr (Model2D*) - use getHestonModel() to cast
  // - _randomEngine - use generateStandardNormalSLV(), generateUniform()
  // - _timeSteps (reference)
  // =========================================================================

  const VolatilitySurface
      *_volSurfacePtr; // Cloned vol surface for Dupire σ_LV(S,t)
  size_t _numPaths;    // N: number of MC paths
  size_t _numBins;     // l: number of bins for E[V|S]
  double _psiC;        // ψ_c: QE switching threshold (default 1.5)

  // Helper to cast _modelPtr to HestonModel*
  const HestonModel *getHestonModel() const;
  /** What does this do?
   *   * Allows PathSimulator2D to work with any Model2D subclass
   * (polymorphism). But Model2D only provide generic methods like drift2D(),
   * diffusion2D(), correlation(), initValue() - getters
   *     HestonSLVPathSimulator2D requires Heston-specific params which are not
   * in the Model2D - i.e. v0(), kappa(), vbar(), ... Solution: getHestonModel()
   * - safely converts the base pointer to the derived type. Helper to cast
   * _modelPtr
   */
};

#endif // CPPFM_PATHSIMULATOR2D_H