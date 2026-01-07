//
// Created by Patrik  Liba on 27/11/2025.
//

#include "PathSimulator2D.h"
#include "Utils.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <tbb/parallel_sort.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// PathSimulator2D Base Class Implementation
// ============================================================================

PathSimulator2D::PathSimulator2D(const std::vector<double> &timeSteps,
                                 const Model2D &model, size_t randomSeed)
    : _timeSteps(timeSteps), _modelPtr(model.clone()), _randomSeed(randomSeed),
      _randomEngine(randomSeed), _antitheticMode(false), _cacheIndex(0)
{
  // _modelPtr should point to a COPY of model (use clone() method)
  // This ensures PathSimulator2D owns its model and model can be destroyed
  // safely

  // Call a method that does SANITY CHECKS on timeSteps
  // timeSteps init value is zero
  // TimeSteps is strictly increasing sequence
  if (!timeStepsSanityCheck())
    throw std::runtime_error("The Time Steps are not correct!"); // exception

  // Reserve space for random cache (estimate: 2 randoms per time step)
  _randomCache.reserve(timeSteps.size() * 4);
}

PathSimulator2D::~PathSimulator2D()
{
  // Delete the model pointer - the 'new' called by clone() is deleted here
  if (_modelPtr)
  { // Check if Null before deleting pointer
    delete _modelPtr;
    _modelPtr = nullptr; // Prevent double deletion
  }
}

bool PathSimulator2D::timeStepsSanityCheck() const
{
  // required more logic
  if (_timeSteps.empty())
    return false;
  if (_timeSteps[0] < 0.0)
    return false;
  for (size_t i = 1; i < _timeSteps.size(); ++i)
  {
    if (!std::isfinite(_timeSteps[i]))
      return false; // finiteness
    if (_timeSteps[i] <= _timeSteps[i - 1])
      return false; // increasing
  }
  return true;
}

/**
 *  PATH-WISE IMPLEMENTATION
 */
// (X, V) put them together
std::pair<std::vector<double>, std::vector<double>>
PathSimulator2D::paths() const
{
  /**
   * Path generation with PATH-LEVEL antithetic sampling
   *
   * For proper antithetic sampling:
   * - Path 1: Generate with fresh randoms, cache them
   * - Path 2: Replay cached randoms negated
   * - Path 3: Generate fresh again
   * - etc.
   */

  // Setup for current path
  if (_antitheticMode)
  {
    // Antithetic mode: replay cached randoms (negated)
    _cacheIndex = 0;
  }
  else
  {
    // Fresh mode: clear cache for new path
    _randomCache.clear();
  }

  // Initialisation stage : Path[0] = S0, Variance[0] = V0
  std::vector<double> assetPath;
  std::vector<double> variancePath;

  assetPath.push_back(_modelPtr->initValue()); // S_0

  // For 2D models like Heston, we need V_0
  // Since Model2D interface doesn't have auxInitValue(), we downcast to
  // HestonModel
  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error(
        "PathSimulator2D currently only supports HestonModel");
  }
  variancePath.push_back(hestonPtr->v0());

  // Iteration stage : Path[i] to Path[i+1]
  for (size_t idx = 0; idx < _timeSteps.size() - 1; ++idx)
  {
    auto [S_next, V_next] = nextStep(idx, assetPath[idx], variancePath[idx]);
    assetPath.push_back(S_next);
    variancePath.push_back(V_next);
  }

  // Toggle antithetic mode for next path
  _antitheticMode = !_antitheticMode;

  return {assetPath, variancePath};
}

std::pair<double, double> PathSimulator2D::generateCorrelatedNormals() const
{
  /**
   * Generate correlated standard normals with PATH-LEVEL antithetic sampling
   *
   * NOTE: This method is used ONLY by Euler scheme (Equations 6-7).
   *
   * Algorithm:
   * 1. If antithetic mode: replay cached Z_V, Z_perp (negated)
   * 2. Otherwise: generate fresh and cache them
   * 3. Apply Cholesky: Z_X = ρ·Z_V + √(1-ρ²)·Z_perp
   */

  double Z_V, Z_perp;

  if (_antitheticMode)
  {
    // Replay cached values (negated)
    Z_V = -_randomCache[_cacheIndex++];
    Z_perp = -_randomCache[_cacheIndex++];
  }
  else
  {
    // Generate fresh normals
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double U1 = uniform(_randomEngine);
    double U2 = uniform(_randomEngine);

    Z_V = Utils::inverseNormalCDF(U1);
    Z_perp = Utils::inverseNormalCDF(U2);

    // Cache for antithetic path
    _randomCache.push_back(Z_V);
    _randomCache.push_back(Z_perp);
  }

  // Apply correlation structure
  double rho = _modelPtr->correlation();
  double Z_X = rho * Z_V + std::sqrt(1.0 - rho * rho) * Z_perp;

  return {Z_V, Z_X};
}

double PathSimulator2D::generateUniform() const
{
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  return uniform(_randomEngine);
}

double PathSimulator2D::generateStandardNormalSLV() const
{ // ! because SLV doesn't use antithetic sampling at path-level
  return Utils::inverseNormalCDF(generateUniform());
}

double PathSimulator2D::generateStandardNormal() const
{
  /**
   * Generate a single standard normal with PATH-LEVEL antithetic sampling
   *
   * Used by TG, QE, and BK schemes for independent normal generation
   *
   * If antithetic mode: replay cached value (negated)
   * Otherwise: generate fresh and cache it
   */

  if (_antitheticMode)
  {
    // Replay cached value (negated)
    return -_randomCache[_cacheIndex++];
  }
  else
  {
    // Generate fresh normal
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double U = uniform(_randomEngine);
    double Z = Utils::inverseNormalCDF(U);

    // Cache for antithetic path
    _randomCache.push_back(Z);

    return Z;
  }
}

void PathSimulator2D::initThreadLocalRNGs() const
{
#ifdef _OPENMP
  int numThreads = omp_get_max_threads();
#else
  int numThreads = 1;
#endif

  _threadRNGs.resize(numThreads);

  for (int t = 0; t < numThreads; ++t)
  {
    // PCG32: use stream parameter for independent sequences per thread
    // Each thread gets the same seed but a different stream ID (0, 1, 2, ...)
    // PCG's stream feature guarantees statistically independent sequences
    _threadRNGs[t] = pcg32(_randomSeed, t);
  }
}

double PathSimulator2D::generateStandardNormalParallel(int threadId) const
{
  // thread-safe RNG
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  double U = uniform(_threadRNGs[threadId]); // thread's private RNG
  return Utils::inverseNormalCDF(U);
}

// ============================================================================
// Step Log Price Eq. 33 - Correlation-Preserving Discretization for X
// ============================================================================

double PathSimulator2D::stepLogPriceEq33(double X_t, double V_t, double V_next,
                                         double dt, double Z, double gamma1,
                                         double gamma2) const
{
  /**
   * Correlation-Preserving Discretization for X (Andersen Eq. 33)
   *
   * Notation: X(t) = ln(S(t)) is the log-price
   *
   * X̂(t + Δ) = X̂(t) + K₀ + K₁V̂(t) + K₂V̂(t+Δ) + √(K₃V̂(t) + K₄V̂(t+Δ)) · Z
   *
   * where Z is a standard Gaussian independent of V̂, and:
   *   K₀ = -(ρκθ/ε)Δ
   *   K₁ = γ₁Δ(κρ/ε - 1/2) - ρ/ε
   *   K₂ = γ₂Δ(κρ/ε - 1/2) + ρ/ε
   *   K₃ = γ₁Δ(1 - ρ²)
   *   K₄ = γ₂Δ(1 - ρ²)
   *
   * This scheme preserves the correlation ρ/ε · V(t+Δ) between X(t+Δ) and
   * V(t+Δ)
   */

  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("stepLogPriceEq33 requires HestonModel");
  }

  double r = hestonPtr->riskFreeRate(); // r: risk-free rate
  double kappa = hestonPtr->kappa();
  double vbar = hestonPtr->vbar();       // θ
  double sigma_v = hestonPtr->sigma_v(); // ε
  double rho = hestonPtr->rho();         // ρ

  // ========================================================================
  // Compute K₀, K₁, K₂, K₃, K₄ coefficients (Andersen paper)
  // ========================================================================

  double K0 =
      -(rho * kappa * vbar / sigma_v) * dt; // interest rate could go here.

  double K1 = gamma1 * dt * (kappa * rho / sigma_v - 0.5) - rho / sigma_v;

  double K2 = gamma2 * dt * (kappa * rho / sigma_v - 0.5) + rho / sigma_v;

  double K3 = gamma1 * dt * (1.0 - rho * rho);

  double K4 = gamma2 * dt * (1.0 - rho * rho);

  // ========================================================================
  // Apply Equation 33
  // ========================================================================

  double drift = r * dt;

  double variance_coeff = K3 * V_t + K4 * V_next;
  double term1 = std::sqrt(std::max(variance_coeff, 0.0)) *
                 Z; // Ensure non-negative argument

  double X_next = X_t + drift + K0 + K1 * V_t + K2 * V_next + term1;

  return X_next; // X(t+Δ) = ln(S(t+Δ))
}

// ============================================================================
// Shared Variance Discretization Methods
// ============================================================================

double PathSimulator2D::stepVarianceTG(double V_t, double dt,
                                       double Z_V) const
{
  /**
   * Truncated Gaussian (TG) Scheme for Variance (Andersen Eq. 13)
   *
   * CIR process: dV = κ(θ - V)dt + σ_v√V dW_V
   *
   * TG scheme: V̂(t+Δ) = [V̂(t) + κ(θ - V̂(t))Δ + σ_v√(V̂(t))√Δ·Z_V]⁺
   *
   * where [x]⁺ = max(x, 0) is the positive part (truncation)
   *
   * Used by: TGPathSimulator2D, BKTGPathSimulator2D
   */

  const HestonModel *hestonPtr = static_cast<const HestonModel *>(_modelPtr);
  double kappa = hestonPtr->kappa();
  double vbar = hestonPtr->vbar();
  double sigma_v = hestonPtr->sigma_v();

  double sqrtDt = std::sqrt(dt);

  // Apply truncation to current variance
  double V_plus = std::max(V_t, 0.0);

  // Drift term: κ(θ - V⁺)Δ
  double drift = kappa * (vbar - V_plus) * dt;

  // Diffusion term: σ_v√(V⁺)√Δ·Z_V
  double diffusion = sigma_v * std::sqrt(V_plus) * sqrtDt * Z_V;

  // Update with truncation
  double V_next = std::max(V_t + drift + diffusion, 0.0); // Eq. 13

  return V_next;
}

double PathSimulator2D::stepVarianceQE(double V_t, double dt, double Z_V,
                                       double psi_c) const
{
  /**
   * Quadratic Exponential (QE) Scheme for Variance
   *
   * CIR process: dV = κ(θ - V)dt + σ_v√V dW_V
   *
   * QE scheme switches between quadratic and exponential approximations
   * based on ψ = s²/m² where:
   *   m = E[V(t+Δ) | V(t)]
   *   s² = Var[V(t+Δ) | V(t)]
   *
   * If ψ ≤ ψ_c: Quadratic scheme
   * If ψ > ψ_c: Exponential scheme
   *
   * Used by: QEPathSimulator2D, BKQEPathSimulator2D
   */

  const HestonModel *hestonPtr = static_cast<const HestonModel *>(_modelPtr);
  double kappa = hestonPtr->kappa();
  double vbar = hestonPtr->vbar();
  double sigma_v = hestonPtr->sigma_v();

  // ========================================================================
  // Compute moments of V(t+Δ) conditional on V(t)
  // ========================================================================

  double exp_kappa_dt = std::exp(-kappa * dt);

  // Mean: m = E[V(t+Δ) | V(t)] = θ + (V(t) - θ)e^(-κΔt)
  double m = vbar + (V_t - vbar) * exp_kappa_dt;

  // Variance: s² = Var[V(t+Δ) | V(t)]
  //           = V(t)·σ²/κ·(e^(-κΔt) - e^(-2κΔt)) + θ·σ²/(2κ)·(1 - e^(-κΔt))²
  double s2 = V_t * sigma_v * sigma_v / kappa *
                  (exp_kappa_dt - exp_kappa_dt * exp_kappa_dt) +
              vbar * sigma_v * sigma_v / (2.0 * kappa) * (1.0 - exp_kappa_dt) *
                  (1.0 - exp_kappa_dt);

  // Scaled variance: ψ = s²/m²
  // Guard against division by zero when mean variance is very small
  const double MIN_MEAN_VARIANCE = 1e-10;
  double m_safe = std::max(m, MIN_MEAN_VARIANCE);
  double psi = s2 / (m_safe * m_safe); // Eq. 19

  // ========================================================================
  // 3.2.3 Switching rule
  // ========================================================================

  // CRITICAL: Generate ALL randoms upfront to ensure consistent consumption
  // This prevents cache misalignment when branch decisions differ between
  // fresh and antithetic paths due to different variance values.
  // NOTE: For parallel execution, the caller must provide Z_V which already
  // encodes the extra randomness needed for the exponential branch.
  // We use Z_V transformed to uniform for the switching logic.
  double U_V = Utils::stdNormCdf(Z_V); // Use Z_V for both quadratic and exponential

  if (psi <= psi_c)
  {
    // ? How does this mix-in with Ant-sampling?
    // --------------------------------------------------------------------
    // Quadratic scheme (ψ ≤ ψ_c): V ~ a(b + Z_V)²
    // --------------------------------------------------------------------
    // Note: U_V is generated but not used (ensures consistent random
    // consumption)

    double psi_inv = 2.0 / psi;
    double temp = std::max(0.0, psi_inv - 1.0);                     // Numerical guard
    double b_squared = temp + std::sqrt(psi_inv) * std::sqrt(temp); // Eq. 27
    double a = m_safe / (1.0 + b_squared);                          // Eq. 28
    double b = std::sqrt(b_squared);

    double V_next = a * (b + Z_V) * (b + Z_V);
    return std::max(V_next, 0.0); // Ensure non-negativity
  }
  else
  {
    // --------------------------------------------------------------------
    // Exponential scheme (ψ > ψ_c): V ~ exponential/point mass mixture
    // --------------------------------------------------------------------
    double p = (psi - 1.0) / (psi + 1.0);
    double beta = (1.0 - p) / m_safe;

    if (U_V <= p)
    { // Eq. 25
      // Point mass at 0
      return 0.0;
    }
    else
    {
      // Exponential: V = (1/β)·ln((1-p)/(1-U))
      double V_next = std::log((1.0 - p) / (1.0 - U_V)) / beta;
      return V_next;
    }
  }
}

// ============================================================================
// Euler Path Simulator 2D - Full Truncation Scheme for Heston Model
// ============================================================================

EulerPathSimulator2D::EulerPathSimulator2D(const std::vector<double> &timeSteps,
                                           const Model2D &model,
                                           size_t randomSeed)
    : PathSimulator2D(timeSteps, model, randomSeed) {}

std::pair<double, double>
EulerPathSimulator2D::nextStep(size_t timeIndex, double assetPrice,
                               double variance) const
{
  /**
   * Full Truncation Euler Scheme (Equations 6-7)
   * Uses CORRELATED Brownian motions for V and X directly.
   */

  // Get time step size Δ = t_{i+1} - t_i
  double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];
  double sqrtDt = std::sqrt(dt);

  // Get Heston model parameters
  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("EulerPathSimulator2D requires HestonModel");
  }
  double r = hestonPtr->riskFreeRate();  // r: risk-free rate
  double kappa = hestonPtr->kappa();     // κ: mean reversion speed
  double vbar = hestonPtr->vbar();       // v̄: long-term variance mean
  double sigma_v = hestonPtr->sigma_v(); // σ_v: volatility of volatility

  // Generate correlated standard normals (Z_V, Z_X) with correlation ρ
  // Both Z_V and Z_X are used in this scheme
  auto [Z_V, Z_X] = generateCorrelatedNormals();

  // FULL TRUNCATION SCHEME: Apply x⁺ = max(x, 0) to variance
  double V_plus = std::max(variance, 0.0);

  // VARIANCE UPDATE (Equation 7):
  // V̂(t + Δ) = V̂(t) + κ(θ - V̂(t)⁺)Δ + ε√(V̂(t)⁺)Z_V√Δ
  double term1V = kappa * (vbar - V_plus) * dt;
  double term2V = sigma_v * std::sqrt(V_plus) * Z_V * sqrtDt;
  double V_next =
      std::max(variance + term1V + term2V,
               0.0); // FULL TRUNCATION SCHEME: Apply x+ = max(x, 0) to variance

  // LOG-PRICE UPDATE (Equation 6):
  // Notation: X(t) = ln(S(t)) is the log-price
  // X̂(t + Δ) = X̂(t) - (1/2)V̂(t)⁺Δ + √(V̂(t)⁺)·Z_X·√Δ
  // Therefore: S(t+Δ) = S(t)·exp(-(1/2)V⁺(t)Δ + √(V⁺(t))·Z_X·√Δ)   //
  // exponentiate
  double drift = r * dt;
  double term1X = -0.5 * V_plus * dt;
  double term2X = std::sqrt(V_plus) * Z_X * sqrtDt;
  double S_next =
      assetPrice * std::exp(drift + term1X + term2X); // S(t+Δ) = S(t)·exp(ΔX)

  return {S_next, V_next};
}

// ============================================================================
// BK Scheme Base - Equation 11 (Exact SDE Representation)
// ============================================================================

BKSchemeBase::BKSchemeBase(const std::vector<double> &timeSteps,
                           const Model2D &model, size_t randomSeed)
    : PathSimulator2D(timeSteps, model, randomSeed) {}

std::pair<double, double> BKSchemeBase::nextStep(size_t timeIndex,
                                                 double assetPrice,
                                                 double variance) const
{
  /**
   * Common nextStep for all BK variants
   *
   * Algorithm:
   * 1. Generate V(t+Δ) and ∫V using variant-specific method
   * 2. Generate independent normal Z for asset price
   * 3. Apply Equation 11 (same for all variants)
   */

  double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

  // ========================================================================
  // STEP 1: Get V_next and integratedV from subclass
  // ========================================================================
  auto [V_next, integratedV] = generateIntegratedVariance(variance, dt);

  // ========================================================================
  // STEP 2: Generate independent Z for asset price
  // ========================================================================
  double Z = generateStandardNormal();

  // ========================================================================
  // STEP 3: Apply Equation 11 (common to all BK variants)
  // ========================================================================
  double X_t = std::log(assetPrice);
  double X_next = stepLogPriceEq11(X_t, variance, V_next, integratedV, dt, Z);
  double S_next = std::exp(X_next);

  return {S_next, V_next};
}

double BKSchemeBase::stepLogPriceEq11(double X_t, double V_t, double V_next,
                                      double integratedV, double dt,
                                      double Z) const
{
  /**
   * Broadie-Kaya Equation 11 (Exact SDE Representation)
   *
   * From Andersen (2008) Equation 11:
   * X(t+Δ) = X(t) + r·Δt + (ρ/ε)(V(t+Δ) - V(t) - κθΔ)
   *                      + (κρ/ε - 1/2)∫V(s)ds
   *                      + √(1-ρ²)·√(∫V)·Z
   *
   * Shared by ALL BK variants (exact and approximate)
   */

  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("stepLogPriceEq11 requires HestonModel");
  }

  double r = hestonPtr->riskFreeRate();
  double kappa = hestonPtr->kappa();
  double vbar = hestonPtr->vbar();
  double sigma_v = hestonPtr->sigma_v();
  double rho = hestonPtr->rho();

  // BK Equation 11 decomposition:
  double drift = r * dt;
  double term1 = (rho / sigma_v) * (V_next - V_t - kappa * vbar * dt);
  double term2 = (kappa * rho / sigma_v - 0.5) * integratedV;
  double term3 = std::sqrt(1.0 - rho * rho) * std::sqrt(integratedV) * Z;

  double X_next = X_t + drift + term1 + term2 + term3;

  return X_next;
}

// ============================================================================
// BK Approximate Scheme - Equations 32 & 33
// ============================================================================

BKApproximateScheme::BKApproximateScheme(const std::vector<double> &timeSteps,
                                         const Model2D &model,
                                         size_t randomSeed, double gamma1,
                                         double gamma2)
    : BKSchemeBase(timeSteps, model, randomSeed), _gamma1(gamma1),
      _gamma2(gamma2) {}

std::pair<double, double>
BKApproximateScheme::generateIntegratedVariance(double V_t, double dt) const
{
  /**
   * Implements approximation (32): ∫V ≈ γ₁V(t) + γ₂V(t+Δ)
   *
   * Uses trapezoidal rule by default (γ₁ = γ₂ = 0.5)
   * This leads to efficient computation via Equation 33
   */

  // Generate V(t+Δ) using subclass-specific method (TG or QE)
  double V_next = generateNextVariance(V_t, dt);

  // Apply approximation (32)
  double integratedV = _gamma1 * V_t + _gamma2 * V_next;
  integratedV *= dt; // Scale by time step

  return {V_next, integratedV};
}

// ============================================================================
// BKExact - Exact CIR + Exact ∫V (Original Broadie-Kaya)
// ============================================================================

BKExactPathSimulator2D::BKExactPathSimulator2D(
    const std::vector<double> &timeSteps, const Model2D &model,
    size_t randomSeed, NewtonMethod newtonMethod)
    : BKSchemeBase(timeSteps, model, randomSeed), _newtonMethod(newtonMethod) {}

std::pair<double, double>
BKExactPathSimulator2D::generateIntegratedVariance(double V_t,
                                                   double dt) const
{
  /**
   * BKExact: Original Broadie-Kaya - bypasses approximation (32)
   *
   * 1. Sample V(t+Δ) exactly from non-central χ² (CIR distribution)
   * 2. Sample ∫[t,t+Δ] V(s)ds exactly using COS method + Newton inversion
   *
   * Most accurate but slowest method
   */

  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("BKExact requires HestonModel");
  }

  double kappa = hestonPtr->kappa();
  double vbar = hestonPtr->vbar();
  double sigma_v = hestonPtr->sigma_v();

  // ========================================================================
  // STEP 1: Sample V(t+Δ) exactly using CIR distribution (with antithetic
  // sampling)
  // ========================================================================
  // Generate uniform via antithetic mechanism: Z → Φ(Z) → Uniform(0,1)
  double Z_cir = generateStandardNormal();
  double u_cir = Utils::stdNormCdf(Z_cir); // Forward CDF: Normal → Uniform

  // Use existing CIR sampler from Utils with pre-generated uniform
  double V_next =
      CIRSampler::sampleCIR(kappa, vbar, sigma_v, V_t, 0.0, dt, u_cir);

  // ========================================================================
  // STEP 2: Sample integrated variance ∫[t,t+Δ] V(s)ds exactly
  // ========================================================================

  // Define the characteristic function for integrated variance
  auto chf = [&](double omega) -> std::complex<double>
  {
    return ChFIntegratedVariance::compute(omega, kappa, vbar, sigma_v, V_t,
                                          V_next, dt);
  };

  // Truncation bounds for COS method [a, b]
  double V_avg = 0.5 * (V_t + V_next);
  double expected_intV = V_avg * dt;
  double std_intV = sigma_v * std::sqrt(V_avg * dt * dt * dt / 3.0);

  double a = std::max(0.0, expected_intV - 10.0 * std_intV);
  double b = expected_intV + 10.0 * std_intV;
  size_t N = 128;

  // ========================================================================
  // Generate uniform for CDF inversion via antithetic mechanism
  // ========================================================================
  // CRITICAL: Use generateStandardNormal() to ensure antithetic sampling works!
  // The uniform is obtained via Φ(Z) where Φ is the standard normal CDF.
  // For antithetic paths: Z₂ = -Z₁ → U₂ = Φ(-Z₁) = 1 - Φ(Z₁) = 1 - U₁

  double Z_u = generateStandardNormal();
  double u = Utils::stdNormCdf(Z_u); // Use existing forward CDF

  // Invert CDF to get sample of integrated variance
  double integratedV;
  if (_newtonMethod == NewtonMethod::Original)
  {
    integratedV = Transforms::invertCDF(a, b, N, chf, u).first;
  }
  else
  {
    integratedV = Transforms::invertCDF_Optimized(a, b, N, chf, u).first;
  }

  integratedV = std::max(0.0, integratedV);

  return {V_next, integratedV};
}

// ============================================================================
// TG Path Simulator 2D (Truncated Gaussian Scheme)
// ============================================================================

TGPathSimulator2D::TGPathSimulator2D(const std::vector<double> &timeSteps,
                                     const Model2D &model, size_t randomSeed,
                                     double gamma1, double gamma2)
    : PathSimulator2D(timeSteps, model, randomSeed), _gamma1(gamma1),
      _gamma2(gamma2)
{
  // Validate γ₁ + γ₂ should typically equal 1 for proper approximation
  if (std::abs(_gamma1 + _gamma2 - 1.0) > 1e-10)
  {
    // Warning: non-standard choice, but allow flexibility
  }
}

TGPathSimulator2D::~TGPathSimulator2D()
{
  // Base class destructor handles cleanup
}

std::pair<double, double> TGPathSimulator2D::nextStep(size_t timeIndex,
                                                      double assetPrice,
                                                      double variance) const
{
  /**
   * TG Scheme:
   *
   * Step 1: Generate V̂(t+Δ) using TG scheme (Eq. 13)
   * Step 2: Draw uniform U ~ Uniform(0,1), independent of V̂(t+Δ)
   * Step 3: Set Z = Φ^(-1)(U) using Moro's algorithm
   * Step 4: Compute X̂(t+Δ) from Equation 33
   *
   * Key insight: Correlation between X and V is preserved through the
   * K coefficients in Eq. 33, NOT through correlated random numbers!
   */

  // Get time step size Δ = t_{i+1} - t_i
  double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

  // Get Heston model parameters
  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("TGPathSimulator2D requires HestonModel");
  }

  // ========================================================================
  // STEP 1: Generate V(t+Δ) using TG scheme (Andersen Eq. 13)
  // ========================================================================

  // Generate standard normal for variance using Moro's inverse CDF
  // This supports antithetic sampling automatically
  double Z_V = generateStandardNormal();

  // Apply TG scheme to get V(t+Δ) using shared implementation
  double V_next = PathSimulator2D::stepVarianceTG(variance, dt, Z_V);

  // ========================================================================
  // STEP 2-3: Generate INDEPENDENT Z for X discretization (Eq. 33)
  // ========================================================================

  // Generate independent standard normal for asset price
  // This also supports antithetic sampling
  double Z = generateStandardNormal();

  // ========================================================================
  // STEP 4: Compute X(t+Δ) using correlation-preserving scheme (Eq. 33)
  // ========================================================================
  // Notation: X(t) = ln(S(t)) is the log-price
  // The correlation between X and V is preserved through K₁ and K₂ coefficients

  double X_t = std::log(assetPrice); // X(t) = ln(S(t))
  double X_next = PathSimulator2D::stepLogPriceEq33(X_t, variance, V_next, dt, Z, _gamma1, _gamma2);
  double S_next = std::exp(X_next); // S(t+Δ) = exp(X(t+Δ))

  return {S_next, V_next};
}

// ============================================================================
// QE Path Simulator 2D
// ============================================================================

QEPathSimulator2D::QEPathSimulator2D(const std::vector<double> &timeSteps,
                                     const Model2D &model, size_t randomSeed,
                                     double psi_c, double gamma1, double gamma2)
    : PathSimulator2D(timeSteps, model, randomSeed), _psi_c(psi_c),
      _gamma1(gamma1), _gamma2(gamma2)
{
  // Validate γ₁ + γ₂ should typically equal 1 for proper approximation
  if (std::abs(_gamma1 + _gamma2 - 1.0) > 1e-10)
  {
    // Warning: non-standard choice, but allow flexibility
  }
}

QEPathSimulator2D::~QEPathSimulator2D()
{
  // Base class destructor handles destruction
}

std::pair<double, double> QEPathSimulator2D::nextStep(size_t timeIndex, double assetPrice, double variance) const
{
  /**
   * QE Scheme:
   *
   * Step 1: Generate V̂(t+Δ) using QE scheme (Eq. 23/26)
   * Step 2: Draw uniform U ~ Uniform(0,1), independent of V̂(t+Δ)
   * Step 3: Set Z = Φ^(-1)(U) using Moro's algorithm
   * Step 4: Compute X̂(t+Δ) from Equation 33
   *
   * Key insight: Correlation between X and V is preserved through the
   * K coefficients in Eq. 33, NOT through correlated random numbers!
   */

  // Get time step size Δ = t_{i+1} - t_i
  double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

  // Get Heston model parameters
  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("QEPathSimulator2D requires HestonModel");
  }

  // ========================================================================
  // STEP 1: Generate V(t+Δ) using QE scheme (Andersen Eq. 23/26)
  // ========================================================================

  // Generate standard normal for variance using Moro's inverse CDF
  // This supports antithetic sampling automatically
  double Z_V = generateStandardNormal();

  // Apply QE scheme to get V(t+Δ) using shared implementation
  double V_next = PathSimulator2D::stepVarianceQE(variance, dt, Z_V, _psi_c);

  // ========================================================================
  // STEP 2-3: Generate INDEPENDENT Z for X discretization (Eq. 33)
  // ========================================================================

  // Generate independent standard normal for asset price
  // This also supports antithetic sampling
  double Z = generateStandardNormal();

  // ========================================================================
  // STEP 4: Compute X(t+Δ) using correlation-preserving scheme (Eq. 33)
  // ========================================================================
  // Notation: X(t) = ln(S(t)) is the log-price
  // The correlation between X and V is preserved through K1 and K2 coefficients

  double X_t = std::log(assetPrice); // X(t) = ln(S(t))
  double X_next = PathSimulator2D::stepLogPriceEq33(X_t, variance, V_next, dt, Z, _gamma1, _gamma2);
  double S_next = std::exp(X_next); // S(t+Δ) = exp(X(t+Δ))

  return {S_next, V_next};
}

// ============================================================================
// BKTG - TG Variance + Approximation (32)
// ============================================================================

BKTGPathSimulator2D::BKTGPathSimulator2D(const std::vector<double> &timeSteps,
                                         const Model2D &model,
                                         size_t randomSeed, double gamma1,
                                         double gamma2)
    : BKApproximateScheme(timeSteps, model, randomSeed, gamma1, gamma2) {}

double BKTGPathSimulator2D::generateNextVariance(double V_t, double dt) const
{
  /**
   * BKTG: Differs from BKQE ONLY in variance generation
   * Uses Truncated Gaussian (TG) scheme via shared implementation
   */
  double Z_V = generateStandardNormal();
  return PathSimulator2D::stepVarianceTG(V_t, dt, Z_V);
}

// ============================================================================
// BKQE - QE Variance + Approximation (32)
// ============================================================================

BKQEPathSimulator2D::BKQEPathSimulator2D(const std::vector<double> &timeSteps,
                                         const Model2D &model,
                                         size_t randomSeed, double psi_c,
                                         double gamma1, double gamma2)
    : BKApproximateScheme(timeSteps, model, randomSeed, gamma1, gamma2),
      _psi_c(psi_c) {}

double BKQEPathSimulator2D::generateNextVariance(double V_t, double dt) const
{
  /**
   * BKQE: Differs from BKTG ONLY in variance generation
   * Uses Quadratic-Exponential (QE) scheme via shared implementation
   */
  double Z_V = generateStandardNormal();
  return PathSimulator2D::stepVarianceQE(V_t, dt, Z_V, _psi_c);
}

// ============================================================================
// Heston Stochastic Local Volatility (SLV) Simulator
// Reference: van der Stoep et al. (2013) "The Heston SLV Model"
// ============================================================================

HestonSLVPathSimulator2D::HestonSLVPathSimulator2D(
    const HestonModel &model,
    const std::vector<double> &timeSteps, // ensure lifetime
    size_t numPaths, size_t numBins, size_t randomSeed)
    : PathSimulator2D(timeSteps, model, randomSeed),
      _hestonLocalVol(std::make_unique<HestonLocalVol>(
          model.kappa(), model.vbar(), model.sigma_v(), model.rho(),
          model.v0(), model.riskFreeRate(), model.initValue())),
      _numPaths(numPaths), _numBins(numBins), _psiC(1.5)
{
  // Validate SLV-specific parameters
  if (_numPaths == 0)
  {
    throw std::runtime_error("Number of paths must be positive");
  }
  if (_numBins == 0 || _numBins > _numPaths)
  {
    throw std::runtime_error("Number of bins must be in (0, numPaths]");
  }
  // timeSteps validation done by base class

  // Initialize spot grid and precompute local vol grid using analytical Dupire
  // This is done at construction for fast lookups during simulation
  initializeLeverageGrid();
}

HestonSLVPathSimulator2D::~HestonSLVPathSimulator2D()
{
  // _hestonLocalVol is automatically cleaned up by unique_ptr
}

size_t HestonSLVPathSimulator2D::findBinIndex(
    double spot, const std::vector<BinData> &bins) const
{

  // Binary search for bin containing spot
  size_t lo = 0, hi = bins.size() - 1;

  while (lo < hi)
  {
    size_t mid = (lo + hi) / 2;
    if (spot < bins[mid].upperBound)
    {
      hi = mid;
    }
    else
    {
      lo = mid + 1;
    }
  }
  return lo;
}

double HestonSLVPathSimulator2D::interpolateConditionalExpectation(double spot, const std::vector<BinData> &bins) const
{
  size_t k = findBinIndex(spot, bins);

  // Linear interpolation between bin midpoints
  if (spot <= bins[0].midpoint)
  {
    return bins[0].conditionalExpectation;
  }
  if (spot >= bins.back().midpoint)
  {
    return bins.back().conditionalExpectation;
  }

  // Find adjacent bins for interpolating
  size_t k_left = (spot < bins[k].midpoint && k > 0) ? k - 1 : k;
  size_t k_right = k_left + 1;
  if (k_right >= bins.size())
    k_right = bins.size() - 1;

  double t = (spot - bins[k_left].midpoint) / (bins[k_right].midpoint - bins[k_left].midpoint);

  return (1 - t) * bins[k_left].conditionalExpectation + t * bins[k_right].conditionalExpectation;
}

// ! Equation 2.9
double HestonSLVPathSimulator2D::leverageSquared(double spot, double time, const std::vector<BinData> &bins) const
{
  // L2(t, S) = σ^2_LV(t, S) / E[V|S]
  // Uses precomputed local vol grid from COS-based analytical Dupire

  // Find time index for grid lookup
  size_t timeIdx = 0;
  for (size_t i = 0; i < _timeSteps.size(); ++i) {
    if (std::abs(_timeSteps[i] - time) < 1e-10) {
      timeIdx = i;
      break;
    }
  }

  double sigma_LV = getLocalVolFromGrid(spot, timeIdx);
  double E_V_S = interpolateConditionalExpectation(spot, bins);

  const double minExpectation = 1e-8;
  E_V_S = std::max(minExpectation, E_V_S);

  return (sigma_LV * sigma_LV) / E_V_S;
}

double HestonSLVPathSimulator2D::stepLogPriceSLV(double X_t, double V_t, double V_next, double leverageSq, double dt,
                                                 double Z) const
{
  const HestonModel *heston = getHestonModel();
  double r = heston->riskFreeRate();
  double kappa = heston->kappa();
  double vbar = heston->vbar();
  double sigma_v = heston->sigma_v();
  double rho = heston->rho();

  // Equation 3.18
  double L = std::sqrt(leverageSq);
  double L_V_dt = leverageSq * V_t * dt;

  double drift = r * dt - 0.5 * L_V_dt;
  double term1 = (rho / sigma_v) * L * (V_next - kappa * vbar * dt + V_t * (kappa * dt - 1));
  double term2 = std::sqrt(1 - rho * rho) * std::sqrt(L_V_dt) * Z;

  return X_t + drift + term1 + term2;
}

const HestonModel *HestonSLVPathSimulator2D::getHestonModel() const
{
  const HestonModel *hestonPtr = dynamic_cast<const HestonModel *>(_modelPtr);
  if (!hestonPtr)
  {
    throw std::runtime_error("This simulator requires HestonModel2D");
  }
  return hestonPtr;
}

// HestonSLVPathSimulator2D::nextStep - not used, throws error
std::pair<double, double> HestonSLVPathSimulator2D::nextStep(
    size_t /*timeIndex*/, double /*assetPrice*/, double /*variance*/) const
{
  throw std::runtime_error(
      "HestonSLVPathSimulator2D::nextStep() not supported - use simulateAllPaths() instead");
}

// ============================================================================
// HestonSLVPathSimulator2D::simulateAllPaths - SLV Monte Carlo simulation - for European options
// ============================================================================
// TODO
std::vector<std::pair<double, double>> HestonSLVPathSimulator2D::simulateAllPaths() const
{
  /**
   * HESTON SLV MONTE CARLO SIMULATION
   * =================================
   * Algorithm 1
   * 1. Generate paths under Heston dynamics
   * 2. At each time step, apply leverage function L(S,t)
   * 3. L(S,t) calibrates local vol to match market implied vol surface
   *
   * WIP
   */

  const HestonModel *heston = getHestonModel();
  double S0 = heston->initValue(); // S_0 from base class ModelBase
  double V0 = heston->v0();        // Initial variance from HestonModel

  // Current state vectors for all paths
  std::vector<double> spots(_numPaths, S0);
  std::vector<double> variances(_numPaths, V0);

  // Temporary storate for next step
  std::vector<double> nextSpots(_numPaths);
  std::vector<double> nextVariances(_numPaths);

  // Time step - all paths together
  for (size_t i = 0; i < _timeSteps.size() - 1; ++i)
  {
    double t = _timeSteps[i];
    double dt = _timeSteps[i + 1] - t;

    // 1 evolve variances for ALL paths
    for (size_t p = 0; p < _numPaths; ++p)
    {
      double Z_V = generateStandardNormalSLV();
      nextVariances[p] = stepVarianceQE(variances[p], dt, Z_V, _psiC); // QE scheme for variance
    }

    // 2 compute bins for current (S, V)-pair
    std::vector<BinData> bins = computeBins(spots, variances);

    // 3 evolve the log-spot for ALL paths using leverage
    for (size_t p = 0; p < _numPaths; ++p)
    {
      double X_t = std::log(spots[p]);
      double leverageSq = leverageSquared(spots[p], t, bins);
      double Z = generateStandardNormalSLV(); // * the non-antithetic one
      double X_next = stepLogPriceSLV(X_t, variances[p], nextVariances[p], leverageSq, dt, Z);

      nextSpots[p] = std::exp(X_next);
    }
    // 4 update state for the next iteration; swap currs to nexts
    std::swap(spots, nextSpots);
    std::swap(variances, nextVariances);
  }
  // pack terminal values
  std::vector<std::pair<double, double>> terminalValues(_numPaths);
  for (size_t p = 0; p < _numPaths; ++p)
  {
    terminalValues[p] = {spots[p], variances[p]};
  }

  return terminalValues;
}

std::vector<std::pair<double, double>> HestonSLVPathSimulator2D::simulateAllPathsParallel() const
{
  const HestonModel *heston = getHestonModel();
  double S0 = heston->initValue();
  double V0 = heston->v0();

  initThreadLocalRNGs();

  std::vector<double> spots(_numPaths, S0);
  std::vector<double> variances(_numPaths, V0);
  std::vector<double> nextSpots(_numPaths);
  std::vector<double> nextVariances(_numPaths);

  // only compute bins on-the-fly if not calibrated
  std::vector<BinData> currentBins;

  for (size_t i = 0; i < _timeSteps.size() - 1; ++i)
  {
    double t = _timeSteps[i];
    double dt = _timeSteps[i + 1] - t;

    if (!_leverageCalibrated) {
      currentBins = computeBinsVectorized(spots, variances);
    }

#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p)
    {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;
#endif
      double Z_V = generateStandardNormalParallel(tid);
      nextVariances[p] = stepVarianceQE(variances[p], dt, Z_V, _psiC);
    }

#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p)
    {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;
#endif
      double Z = generateStandardNormalParallel(tid);
      double X_t = std::log(spots[p]);

      double leverageSq;
      if (_leverageCalibrated) {
        leverageSq = getLeverageSquaredFromGrid(spots[p], i);
      } else {
        leverageSq = leverageSquared(spots[p], t, currentBins);
      }

      double X_next = stepLogPriceSLV(X_t, variances[p], nextVariances[p], leverageSq, dt, Z);
      nextSpots[p] = std::exp(X_next);
    }

    std::swap(spots, nextSpots);
    std::swap(variances, nextVariances);
  }

  std::vector<std::pair<double, double>> terminalValues(_numPaths);

#pragma omp parallel for
  for (size_t p = 0; p < _numPaths; ++p)
  {
    terminalValues[p] = {spots[p], variances[p]};
  }

  return terminalValues;
}

// ============================================================================
// HestonSLVPathSimulator2D::simulateAllPathsFull - Full path history - for Asian, barrier, lookback options
// ============================================================================
std::vector<std::vector<std::pair<double, double>>> HestonSLVPathSimulator2D::simulateAllPathsFull() const
{

  // Returns full path history: [path_index][time_index] -> (S, V)

  const HestonModel *heston = getHestonModel();
  double S0 = heston->initValue(); // S_0 from base class ModelBase
  double V0 = heston->v0();        // Initial variance from HestonModel

  size_t numSteps = _timeSteps.size();

  // 1 allocate space to store full path: [path_index][time_index] -> (S, V)
  std::vector<std::vector<std::pair<double, double>>> allPaths(_numPaths, std::vector<std::pair<double, double>>(numSteps));

  // current state vectors
  std::vector<double> spots(_numPaths, S0);
  std::vector<double> variances(_numPaths, V0);

  // store initial values
  for (size_t p = 0; p < _numPaths; ++p)
  {
    allPaths[p][0] = {S0, V0};
  }

  // temporary storage of nexts
  std::vector<double> nextSpots(_numPaths);
  std::vector<double> nextVariances(_numPaths);

  // time step - all paths together
  for (size_t i = 0; i < numSteps - 1; ++i)
  {
    double t = _timeSteps[i];
    double dt = _timeSteps[i + 1] - t;

    // 1 evolve variances for ALL paths
    for (size_t p = 0; p < _numPaths; ++p)
    {
      double Z_V = generateStandardNormalSLV();
      nextVariances[p] = stepVarianceQE(variances[p], dt, Z_V, _psiC); // QE scheme for variance
    }

    // 2 compute bins for current (S, V)-pair
    std::vector<BinData> bins = computeBins(spots, variances); // ? change back to parallel

    // 3 evolve the log-spot for ALL paths using leverage
    for (size_t p = 0; p < _numPaths; ++p)
    {
      double X_t = std::log(spots[p]);
      double leverageSq = leverageSquared(spots[p], t, bins);
      double Z = generateStandardNormalSLV(); // * the non-antithetic one
      double X_next = stepLogPriceSLV(X_t, variances[p], nextVariances[p], leverageSq, dt, Z);
      nextSpots[p] = std::exp(X_next);

      // store in path history
      allPaths[p][i + 1] = {nextSpots[p], nextVariances[p]};
    }

    // 4 update state for the next iteration; swap currs to nexts
    std::swap(spots, nextSpots);
    std::swap(variances, nextVariances);
  }

  return allPaths;
}

std::vector<std::vector<std::pair<double, double>>> HestonSLVPathSimulator2D::simulateAllPathsFullParallel() const
{
  const HestonModel *heston = getHestonModel();
  double S0 = heston->initValue();
  double V0 = heston->v0();

  size_t numSteps = _timeSteps.size();

  initThreadLocalRNGs();

  std::vector<std::vector<std::pair<double, double>>> allPaths(_numPaths, std::vector<std::pair<double, double>>(numSteps));

  std::vector<double> spots(_numPaths, S0);
  std::vector<double> variances(_numPaths, V0);
  std::vector<double> nextSpots(_numPaths);
  std::vector<double> nextVariances(_numPaths);

  std::vector<BinData> currentBins;

#pragma omp parallel for schedule(static)
  for (size_t p = 0; p < _numPaths; ++p)
  {
    allPaths[p][0] = {S0, V0};
  }

  for (size_t i = 0; i < numSteps - 1; ++i)
  {
    double t = _timeSteps[i];
    double dt = _timeSteps[i + 1] - t;

    if (!_leverageCalibrated) {
      currentBins = computeBinsVectorized(spots, variances);
    }

#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p)
    {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;
#endif
      double Z_V = generateStandardNormalParallel(tid);
      nextVariances[p] = stepVarianceQE(variances[p], dt, Z_V, _psiC);
    }

#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p)
    {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;
#endif
      double Z = generateStandardNormalParallel(tid);
      double X_t = std::log(spots[p]);

      double leverageSq;
      if (_leverageCalibrated) {
        leverageSq = getLeverageSquaredFromGrid(spots[p], i);
      } else {
        leverageSq = leverageSquared(spots[p], t, currentBins);
      }

      double X_next = stepLogPriceSLV(X_t, variances[p], nextVariances[p], leverageSq, dt, Z);

      nextSpots[p] = std::exp(X_next);
      allPaths[p][i + 1] = {nextSpots[p], nextVariances[p]};
    }
    std::swap(spots, nextSpots);
    std::swap(variances, nextVariances);
  }

  return allPaths;
}

std::vector<HestonSLVPathSimulator2D::BinData> HestonSLVPathSimulator2D::computeBins(
    const std::vector<double> &spotValues, const std::vector<double> &varianceValues) const
{

  size_t n = spotValues.size();

  // Create indices
  std::vector<size_t> sortedIndices(n);
  std::iota(sortedIndices.begin(), sortedIndices.end(), 0);

  // Sort by spot values
  std::sort(sortedIndices.begin(), sortedIndices.end(),
            [&spotValues](size_t a, size_t b)
            {
              return spotValues[a] < spotValues[b];
            });

  // Divide into bins
  std::vector<BinData> bins(_numBins);
  size_t pathsPerBin = n / _numBins;
  size_t remainder = n % _numBins;

  size_t currentIndex = 0;
  for (size_t k = 0; k < _numBins; ++k)
  {
    size_t binSize = pathsPerBin + (k < remainder ? 1 : 0);

    double sumVariance = 0.0;
    double sumSpot = 0.0;
    double minSpot = std::numeric_limits<double>::infinity();
    double maxSpot = -std::numeric_limits<double>::infinity();

    for (size_t j = 0; j < binSize; ++j)
    {
      size_t idx = sortedIndices[currentIndex + j];
      double spot = spotValues[idx];
      sumVariance += varianceValues[idx];
      sumSpot += spot;
      minSpot = std::min(minSpot, spot);
      maxSpot = std::max(maxSpot, spot);
    }

    bins[k].lowerBound = minSpot;
    bins[k].upperBound = maxSpot;
    // Sample mean per van der Stoep Algorithm 1
    bins[k].midpoint = sumSpot / binSize;
    bins[k].conditionalExpectation = sumVariance / binSize; // E[V|S]

    currentIndex += binSize;
  }

  // Last bin upper bound should be +infinity for safety
  bins.back().upperBound = std::numeric_limits<double>::infinity();

  return bins;
}

std::vector<HestonSLVPathSimulator2D::BinData> HestonSLVPathSimulator2D::computeBinsVectorized(const std::vector<double> &spotValues, const std::vector<double> &varianceValues) const
{
  size_t n = spotValues.size();

  // create indices
  std::vector<size_t> sortedIndices(n);
  std::iota(sortedIndices.begin(), sortedIndices.end(), 0);

  // Parallel sort - uses all CPU cores via TBB
  tbb::parallel_sort(sortedIndices.begin(),
                     sortedIndices.end(),
                     [&spotValues](size_t a, size_t b)
                     {
                       return spotValues[a] < spotValues[b];
                     });

  std::vector<BinData> bins(_numBins);
  size_t pathsPerBin = n / _numBins;
  size_t remainder = n % _numBins;

  // precompute the start and end position of the bin
  std::vector<size_t> binStarts(_numBins + 1);
  binStarts[0] = 0;

  for (size_t k = 0; k < _numBins; ++k)
  {
    // first 'remainder' bins get one extra path to handle leftovers
    size_t binSize = pathsPerBin + (k < remainder ? 1 : 0);
    binStarts[k + 1] = binStarts[k] + binSize;
  }
// each bin is independent -> parallelize it
#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < _numBins; ++k)
  {
    // get bin's range in the sorted array
    size_t start = binStarts[k];
    size_t end = binStarts[k + 1];
    size_t binSize = end - start;

    // accumulators for this bin
    double sumVariance = 0.0;
    double sumSpot = 0.0;
    double minSpot = std::numeric_limits<double>::infinity();
    double maxSpot = -std::numeric_limits<double>::infinity();

    // sum over all paths in this bin
    for (size_t j = start; j < end; ++j)
    {
      size_t idx = sortedIndices[j];
      double spot = spotValues[idx];

      sumVariance += varianceValues[idx];
      sumSpot += spot;
      minSpot = std::min(minSpot, spot);
      maxSpot = std::max(maxSpot, spot);
    }
    bins[k].lowerBound = minSpot;
    bins[k].upperBound = maxSpot;
    // Sample mean per van der Stoep Algorithm 1
    bins[k].midpoint = sumSpot / static_cast<double>(binSize);
    bins[k].conditionalExpectation = sumVariance / static_cast<double>(binSize);
  }
  bins.back().upperBound = std::numeric_limits<double>::infinity();

  return bins;
}

// Leverage calibration
void HestonSLVPathSimulator2D::initializeLeverageGrid()
{
  const HestonModel* heston = getHestonModel();
  double S0 = heston->initValue();

  _nSpotGrid = 100;   // Good resolution for leverage interpolation
  double spotMin = 0.3 * S0;
  double spotMax = 2.5 * S0;

  // log-uniform spot grid
  _spotGridPoints.resize(_nSpotGrid);
  double logMin = std::log(spotMin);
  double logMax = std::log(spotMax);
  for (size_t i = 0; i < _nSpotGrid; ++i) {
    double frac = static_cast<double>(i) / static_cast<double>(_nSpotGrid - 1);
    _spotGridPoints[i] = std::exp(logMin + frac * (logMax - logMin));
  }

  // L^2 = 1 everywhere (pure Heston)
  size_t numTimeSteps = _timeSteps.size();
  _leverageGrid.resize(numTimeSteps * _nSpotGrid, 1.0);

  // Precompute local vol grid using COS-based analytical Dupire
  precomputeLocalVolGrid();

  _leverageCalibrated = false;
  _calibrationErrors.clear();
}

void HestonSLVPathSimulator2D::resetCalibration()
{
  _leverageCalibrated = false;
  _calibrationErrors.clear();
  std::fill(_leverageGrid.begin(), _leverageGrid.end(), 1.0);
}

void HestonSLVPathSimulator2D::precomputeLocalVolGrid()
{
  /**
   * Precompute local vol on the same grid as leverage using vectorized COS
   *
   * Performance optimization:
   * - Uses batch localVolsAtTime() which vectorizes COS pricing
   * - Parallelizes across time slices using OpenMP
   * - Reduces CF evaluations from O(N_time × N_spot × N_terms × 13)
   *   to O(N_time × 5 × N_terms) - roughly 100x faster!
   */
  size_t numTimeSteps = _timeSteps.size();
  _localVolGrid.resize(numTimeSteps * _nSpotGrid);

  double sqrtV0 = std::sqrt(getHestonModel()->v0());

  // Parallel loop over time slices
  // Each time slice is independent - perfect for parallelization
  #pragma omp parallel for schedule(dynamic)
  for (size_t t_idx = 0; t_idx < numTimeSteps; ++t_idx) {
    double t = _timeSteps[t_idx];

    // Handle t=0 case (use sqrt(v0) as default)
    if (t < 1e-6) {
      for (size_t s_idx = 0; s_idx < _nSpotGrid; ++s_idx) {
        _localVolGrid[t_idx * _nSpotGrid + s_idx] = sqrtV0;
      }
      continue;
    }

    // Batch compute local vols for all spots at this time using vectorized COS
    std::vector<double> lvs = _hestonLocalVol->localVolsAtTime(_spotGridPoints, t);

    // Copy to grid
    for (size_t s_idx = 0; s_idx < _nSpotGrid; ++s_idx) {
      _localVolGrid[t_idx * _nSpotGrid + s_idx] = lvs[s_idx];
    }
  }
}

double HestonSLVPathSimulator2D::getLocalVolFromGrid(double spot, size_t timeIdx) const
{
  // Linear interpolation in spot dimension
  size_t s_lo = findSpotGridIndex(spot);
  size_t s_hi = std::min(s_lo + 1, _nSpotGrid - 1);

  double spot_lo = _spotGridPoints[s_lo];
  double spot_hi = _spotGridPoints[s_hi];

  double lv_lo = _localVolGrid[timeIdx * _nSpotGrid + s_lo];
  double lv_hi = _localVolGrid[timeIdx * _nSpotGrid + s_hi];

  if (s_lo == s_hi || std::abs(spot_hi - spot_lo) < 1e-12) {
    return lv_lo;
  }

  double t = (spot - spot_lo) / (spot_hi - spot_lo);
  return (1.0 - t) * lv_lo + t * lv_hi;
}

size_t HestonSLVPathSimulator2D::findSpotGridIndex(double spot) const
{
  auto it = std::lower_bound(_spotGridPoints.begin(), _spotGridPoints.end(), spot);

  if (it == _spotGridPoints.begin()) {
    return 0;
  }
  if (it == _spotGridPoints.end()) {
    return _nSpotGrid - 2;
  }

  return static_cast<size_t>(std::distance(_spotGridPoints.begin(), it)) - 1;
}

double HestonSLVPathSimulator2D::getLeverageSquaredFromGrid(double spot, size_t timeIdx) const
{
  if (!_leverageCalibrated || _leverageGrid.empty()) {
    return 1.0;
  }

  size_t t_idx = std::min(timeIdx, _timeSteps.size() - 1);

  size_t s_lo = findSpotGridIndex(spot);
  size_t s_hi = std::min(s_lo + 1, _nSpotGrid - 1);

  double S_lo = _spotGridPoints[s_lo];
  double S_hi = _spotGridPoints[s_hi];
  double alpha = 0.0;
  if (S_hi > S_lo) {
    alpha = (spot - S_lo) / (S_hi - S_lo);
    alpha = std::clamp(alpha, 0.0, 1.0);
  }

  double L2_lo = _leverageGrid[t_idx * _nSpotGrid + s_lo];
  double L2_hi = _leverageGrid[t_idx * _nSpotGrid + s_hi];

  return (1.0 - alpha) * L2_lo + alpha * L2_hi;
}

std::vector<std::vector<HestonSLVPathSimulator2D::BinData>>
HestonSLVPathSimulator2D::simulateAndCollectAllBins() const
{
  const HestonModel* heston = getHestonModel();
  double S0 = heston->initValue();
  double V0 = heston->v0();

  size_t numTimeSteps = _timeSteps.size();
  std::vector<std::vector<BinData>> allTimeBins(numTimeSteps);

  std::vector<double> spots(_numPaths, S0);
  std::vector<double> variances(_numPaths, V0);
  std::vector<double> nextSpots(_numPaths);
  std::vector<double> nextVariances(_numPaths);

  allTimeBins[0] = computeBins(spots, variances);

  initThreadLocalRNGs();

  for (size_t i = 0; i < numTimeSteps - 1; ++i) {
    double t = _timeSteps[i];
    double dt = _timeSteps[i + 1] - t;

    const std::vector<BinData>& currentBins = allTimeBins[i];

    #pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p) {
      #ifdef _OPENMP
        int tid = omp_get_thread_num();
      #else
        int tid = 0;
      #endif
      double Z_V = generateStandardNormalParallel(tid);
      nextVariances[p] = stepVarianceQE(variances[p], dt, Z_V, _psiC);
    }

    #pragma omp parallel for schedule(static)
    for (size_t p = 0; p < _numPaths; ++p) {
      #ifdef _OPENMP
        int tid = omp_get_thread_num();
      #else
        int tid = 0;
      #endif

      double X_t = std::log(spots[p]);

      double leverageSq;
      if (_leverageCalibrated) {
        leverageSq = getLeverageSquaredFromGrid(spots[p], i);
      } else {
        leverageSq = leverageSquared(spots[p], t, currentBins);
      }

      double Z = generateStandardNormalParallel(tid);
      double X_next = stepLogPriceSLV(X_t, variances[p], nextVariances[p], leverageSq, dt, Z);
      nextSpots[p] = std::exp(X_next);
    }

    std::swap(spots, nextSpots);
    std::swap(variances, nextVariances);

    allTimeBins[i + 1] = computeBinsVectorized(spots, variances);
  }

  return allTimeBins;
}

double HestonSLVPathSimulator2D::updateLeverageGrid(
    const std::vector<std::vector<BinData>>& allTimeBins,
    double damping)
{
  double maxRelChange = 0.0;
  size_t numTimeSteps = _timeSteps.size();

  for (size_t t_idx = 1; t_idx < numTimeSteps; ++t_idx) {
    double t = _timeSteps[t_idx];
    const std::vector<BinData>& bins = allTimeBins[t_idx];

    for (size_t s_idx = 0; s_idx < _nSpotGrid; ++s_idx) {
      double spot = _spotGridPoints[s_idx];

      // Use precomputed local vol grid
      double sigma_LV = _localVolGrid[t_idx * _nSpotGrid + s_idx];

      double E_V_S = interpolateConditionalExpectation(spot, bins);
      E_V_S = std::max(1e-8, E_V_S);

      double L2_computed = (sigma_LV * sigma_LV) / E_V_S;
      L2_computed = std::max(L2_computed, 0.0);  // Ensure non-negative

      size_t gridIdx = t_idx * _nSpotGrid + s_idx;
      double L2_old = _leverageGrid[gridIdx];

      // Apply damping directly to L² for linear averaging in variance space
      // This avoids the non-linear cross-terms from (αL₁ + (1-α)L₂)²
      double L2_new = damping * L2_computed + (1.0 - damping) * L2_old;

      // Compute relative change using L (not L²) for convergence metric
      double L_old = std::sqrt(std::max(L2_old, 0.0));
      double L_new = std::sqrt(std::max(L2_new, 0.0));

      double relChange = (L_old > 1e-10)
                         ? std::abs(L_new - L_old) / L_old
                         : std::abs(L_new - L_old);
      maxRelChange = std::max(maxRelChange, relChange);

      _leverageGrid[gridIdx] = L2_new;
    }
  }

  // Copy leverage from t_idx=1 to t_idx=0 (Dupire formula unreliable at t->0)
  // This ensures the first simulation step uses calibrated leverage
  if (numTimeSteps > 1) {
    for (size_t s_idx = 0; s_idx < _nSpotGrid; ++s_idx) {
      _leverageGrid[0 * _nSpotGrid + s_idx] = _leverageGrid[1 * _nSpotGrid + s_idx];
    }
  }

  return maxRelChange;
}

size_t HestonSLVPathSimulator2D::calibrateLeverage(size_t maxIterations, double tol, double damping)
{
  initializeLeverageGrid();

  // CRITICAL: Keep _leverageCalibrated = false during the iteration loop
  // so that simulateAndCollectAllBins() always computes fresh leverage from formula.
  // Only set it true AFTER calibration completes for subsequent pricing calls.
  _leverageCalibrated = false;

  size_t iter = 0;
  for (; iter < maxIterations; ++iter) {
    auto allTimeBins = simulateAndCollectAllBins();

    double maxChange = updateLeverageGrid(allTimeBins, damping);
    _calibrationErrors.push_back(maxChange);

    if (maxChange < tol) {
      _leverageCalibrated = true;  // Only now mark as calibrated
      return iter + 1;
    }
  }

  _leverageCalibrated = true;  // Mark as calibrated even if max iterations reached
  return iter;
}

