//
// Created by Patrik  Liba on 27/11/2025.
//

#include "PathSimulator2D.h"
#include "Utils.h"
#include <stdexcept>
#include <cmath>

// ============================================================================
// PathSimulator2D Base Class Implementation
// ============================================================================

PathSimulator2D::PathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model, size_t randomSeed)
    : _timeSteps(timeSteps), _modelPtr(model.clone()), _randomSeed(randomSeed), _randomEngine(randomSeed)
{
    // _modelPtr should point to a COPY of model (use clone() method)
    // This ensures PathSimulator2D owns its model and model can be destroyed safely
    
    // Call a method that does SANITY CHECKS on timeSteps
    // TimeSteps init value is zero
    // TimeSteps is strictly increasing sequence
    if (!timeStepsSanityCheck())
        throw std::runtime_error("The Time Steps are not correct!"); // exception
} 

PathSimulator2D::~PathSimulator2D()
{
    // Delete the model pointer - the 'new' called by clone() is deleted here
    if (_modelPtr) {  // Check if Null before deleting pointer
        delete _modelPtr;
        _modelPtr = nullptr;  // Prevent double deletion
    }
}

bool PathSimulator2D::timeStepsSanityCheck() const {
    // required more logic
    if (_timeSteps.empty()) return false;
    if (_timeSteps[0] < 0.0) return false;
    for (size_t i = 1; i < _timeSteps.size(); ++i)
    {
        if (!std::isfinite(_timeSteps[i])) return false; // finiteness
        if (_timeSteps[i] <= _timeSteps[i-1]) return false; // increasing
    }
    return true;
}

std::pair<std::vector<double>, std::vector<double>> PathSimulator2D::paths() const
{
    // Initialisation stage : Path[0] = S0, Variance[0] = V0
    std::vector<double> assetPath;
    std::vector<double> variancePath;
    
    assetPath.push_back(_modelPtr->initValue());  // S_0
    
    // For 2D models like Heston, we need V_0
    // Since Model2D interface doesn't have auxInitValue(), we downcast to HestonModel
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("PathSimulator2D currently only supports HestonModel");
    }
    variancePath.push_back(hestonPtr->v0());
    
    // Iteration stage : Path[i] to Path[i+1]
    for (size_t idx = 0; idx < _timeSteps.size() - 1; ++idx) {
        auto [S_next, V_next] = nextStep(idx, assetPath[idx], variancePath[idx]);
        assetPath.push_back(S_next);
        variancePath.push_back(V_next);
    }
    
    return {assetPath, variancePath};
}

std::pair<double, double> PathSimulator2D::generateCorrelatedNormals() const
{
    /**
     * Generate correlated standard normals using Moro's inverse CDF
     * 
     * NOTE: This method is used ONLY by Euler scheme (Equations 6-7).
     *       TG and QE schemes do NOT use this method - they generate
     *       independent uniforms separately for V and X.
     * 
     * Algorithm (standard Cholesky decomposition):
     * 1. Generate two INDEPENDENT uniforms U1, U2 ~ Uniform(0, 1)
     * 2. Transform to standard normals using Moro's algorithm: 
     *    Z_V = Φ^(-1)(U1), Z_perp = Φ^(-1)(U2)
     * 3. Apply correlation via Cholesky decomposition: 
     *    Z_X = ρ·Z_V + √(1-ρ²)·Z_perp
     * 
     * Returns: (Z_V, Z_X) where:
     * - Z_V: standard normal for variance process
     * - Z_X: standard normal for asset price (correlated with Z_V by ρ)
     */
    
    // Step 1: Generate two independent uniforms
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double U1 = uniform(_randomEngine);
    double U2 = uniform(_randomEngine);
    
    // Step 2: Inverse CDF transform using Moro's algorithm (Φ^(-1))
    double Z_V = Utils::inverseNormalCDF(U1);       // For variance process
    double Z_perp = Utils::inverseNormalCDF(U2);    // Independent normal
    
    // Step 3: Apply correlation structure - Cholesky decomposition
    // Z_X = ρ·Z_V + √(1-ρ²)·Z_perp
    // This ensures Cov(Z_V, Z_X) = ρ
    double rho = _modelPtr->correlation();
    double Z_X = rho * Z_V + std::sqrt(1.0 - rho * rho) * Z_perp;
    
    return {Z_V, Z_X};  
}

// ============================================================================
// Step Log Price Eq. 33 - Correlation-Preserving Discretization for ln X
// ============================================================================

double PathSimulator2D::stepLogPriceEq33(double X_t, double V_t, double V_next, double dt, double Z,
                                          double gamma1, double gamma2) const
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
     * This scheme preserves the correlation ρ/ε · V(t+Δ) between X(t+Δ) and V(t+Δ)
     */
    
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("stepLogPriceEq33 requires HestonModel");
    }

    double r = hestonPtr->riskFreeRate();    // r: risk-free rate
    double kappa = hestonPtr->kappa();
    double vbar = hestonPtr->vbar();        // θ
    double sigma_v = hestonPtr->sigma_v();  // ε
    double rho = hestonPtr->rho();          // ρ
    
    // ========================================================================
    // Compute K₀, K₁, K₂, K₃, K₄ coefficients (Andersen paper)
    // ========================================================================
    
    double K0 = - (rho * kappa * vbar / sigma_v) * dt; // interest rate could go here.
    
    double K1 = gamma1 * dt * (kappa * rho / sigma_v - 0.5) - rho / sigma_v;
    
    double K2 = gamma2 * dt * (kappa * rho / sigma_v - 0.5) + rho / sigma_v;
    
    double K3 = gamma1 * dt * (1.0 - rho * rho);
    
    double K4 = gamma2 * dt * (1.0 - rho * rho);
    
    // ========================================================================
    // Apply Equation 33
    // ========================================================================
    
    double drift = r * dt;
    
    double variance_coeff = K3 * V_t + K4 * V_next;
    double term1 = std::sqrt(std::max(variance_coeff, 0.0)) * Z; // Ensure non-negative argument
    
    double X_next = X_t + drift + K0 + K1 * V_t + K2 * V_next + term1;
    
    return X_next;     // X(t+Δ) = ln(S(t+Δ))
}

// ============================================================================
// Euler Path Simulator 2D - Full Truncation Scheme for Heston Model
// ============================================================================

EulerPathSimulator2D::EulerPathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model, size_t randomSeed) 
    : PathSimulator2D(timeSteps, model, randomSeed)
{
}

std::pair<double, double> EulerPathSimulator2D::nextStep(size_t timeIndex, double assetPrice, double variance) const
{
    /**
     * Full Truncation Euler Scheme (Equations 6-7)
     * Uses CORRELATED Brownian motions for V and X directly.
     */
    
    // Get time step size Δ = t_{i+1} - t_i
    double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];
    double sqrtDt = std::sqrt(dt);
    
    // Get Heston model parameters
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("EulerPathSimulator2D requires HestonModel");
    }
    double r = hestonPtr->riskFreeRate();    // r: risk-free rate
    double kappa = hestonPtr->kappa();      // κ: mean reversion speed
    double vbar = hestonPtr->vbar();        // v̄: long-term variance mean
    double sigma_v = hestonPtr->sigma_v();  // σ_v: volatility of volatility
    
    // Generate correlated standard normals (Z_V, Z_X) with correlation ρ
    // Both Z_V and Z_X are used in this scheme
    auto [Z_V, Z_X] = generateCorrelatedNormals();
    
    // FULL TRUNCATION SCHEME: Apply x⁺ = max(x, 0) to variance
    double V_plus = std::max(variance, 0.0);
    
    // VARIANCE UPDATE (Equation 7):
    // V̂(t + Δ) = V̂(t) + κ(θ - V̂(t)⁺)Δ + ε√(V̂(t)⁺)Z_V√Δ
    double term1V = kappa * (vbar - V_plus) * dt;
    double term2V = sigma_v * std::sqrt(V_plus) * Z_V * sqrtDt;
    double V_next = std::max(variance + term1V + term2V, 0.0);     // FULL TRUNCATION SCHEME: Apply x+ = max(x, 0) to variance
    
    // LOG-PRICE UPDATE (Equation 6):
    // Notation: X(t) = ln(S(t)) is the log-price
    // X̂(t + Δ) = X̂(t) - (1/2)V̂(t)⁺Δ + √(V̂(t)⁺)·Z_X·√Δ
    // Therefore: S(t+Δ) = S(t)·exp(-(1/2)V⁺(t)Δ + √(V⁺(t))·Z_X·√Δ)   // exponentiate
    double drift = r * dt;
    double term1X = -0.5 * V_plus * dt;
    double term2X = std::sqrt(V_plus) * Z_X * sqrtDt;
    double S_next = assetPrice * std::exp(drift + term1X + term2X);   // S(t+Δ) = S(t)·exp(ΔX)
    
    return {S_next, V_next};
}

// ============================================================================
// BK Path Simulator 2D
// ============================================================================

BKPathSimulator2D::BKPathSimulator2D(const std::vector<double> &timeSteps, const Model2D &model, size_t randomSeed, NewtonMethod newtonMethod)
    : PathSimulator2D(timeSteps, model, randomSeed), _newtonMethod(newtonMethod)
{
}

std::pair<double, double> BKPathSimulator2D::nextStep(size_t timeIndex, double assetPrice, double variance) const
{   /// CAUTION: Notation X(t) = ln(S(t)) is the log-price; comments have ln(X(t)) instead of X(t)
    /**
     * Broadie-Kaya Exact Simulation Scheme for Heston Model
     * 
     * Following Andersen (2008) description of BK scheme:
     * 
     * Algorithm (from paper):
     * 1. Sample V(t+Δ) from non-central chi-squared distribution (exact)
     * 2. Conditional on V(t+Δ) and V(t), draw a sample ∫[t,t+Δ] V(u)du 
     *    using characteristic function + numerical Fourier inversion (COS method)
     * 3. Conditional on V(t+Δ) and ∫V(u)du, use equation (11) to draw ln X(t+Δ)
     *    from a Gaussian distribution
     * 
     * Key insight: ln X(t+Δ) | V(t+Δ), ∫V is Gaussian with easily computable moments
     */
    
    // Get time step size Δ = t_{i+1} - t_i
    double dt = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];
    double t_current = _timeSteps[timeIndex];
    double t_next = _timeSteps[timeIndex + 1];
    
    // Get Heston model parameters
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("BKPathSimulator2D requires HestonModel");
    }

    double kappa = hestonPtr->kappa();       // κ: mean reversion speed
    double vbar = hestonPtr->vbar();         // θ: long-term variance mean
    double sigma_v = hestonPtr->sigma_v();   // σ_v = ε: volatility of volatility
    double rho = hestonPtr->rho();           // ρ: correlation
    double r = hestonPtr->riskFreeRate();    // r: risk-free rate

   // ========================================================================
   // STEP 1: Sample V(t+Δ) exactly using CIR distribution (non-central χ²)
   // ========================================================================
   // V(t+Δ) | V(t) ~ (1/2c) · χ²(δ, λ) where:
   //   δ = 4κθ/ε²  (degrees of freedom)
   //   λ = 2c·exp(-κΔ)·V(t)  (noncentrality parameter)
   //   c = 2κ/(ε²(1-exp(-κΔ)))
   
   double V_t = variance;  // Current variance V(t)
   double V_next = CIRSampler::sampleCIR(kappa, vbar, sigma_v, V_t, t_current, t_next, _randomEngine);
   
   // ========================================================================
   // STEP 2: Sample integrated variance ∫[t,t+Δ] V(u)du exactly
   // ========================================================================
   // Conditional on V(t) and V(t+Δ), sample ∫V using its characteristic function
   // The ChF is known analytically (Broadie-Kaya) and we use COS + Newton inversion
   
   // Define the characteristic function as a lambda (captures all necessary parameters)
   auto chf = [&](double omega) -> std::complex<double> {
       return ChFIntegratedVariance::compute(omega, kappa, vbar, sigma_v, V_t, V_next, dt);
   };
   
   // Truncation bounds for COS method [a, b]
   // We need to choose [a, b] such that most of the probability mass is captured
   // Heuristic: use E[∫V] ± L·SD[∫V] where L = 10 is a safety factor
   //   E[∫V] ≈ (V(t) + V(t+Δ))/2 · Δ  (trapezoidal approximation)
   //   SD[∫V] ≈ ε·√(V_avg·Δ³/3)  (rough estimate from CIR variance)
   
   double V_avg = 0.5 * (V_t + V_next);
   double expected_intV = V_avg * dt;
   double std_intV = sigma_v * std::sqrt(V_avg * dt * dt * dt / 3.0);
   
   double a = std::max(0.0, expected_intV - 10.0 * std_intV);  // Lower bound (≥ 0)
   double b = expected_intV + 10.0 * std_intV;                      // Upper bound
   
   // Number of terms in COS expansion (trade-off: accuracy vs speed)
   size_t N = 128;
   
   // Generate uniform random variable for CDF inversion
   std::uniform_real_distribution<double> uniform(0.0, 1.0);
   double u = uniform(_randomEngine);
   
   // Invert CDF to get sample of integrated variance
   double integratedV;
   if (_newtonMethod == NewtonMethod::Original) {
       integratedV = COSMethod::invertCDF(a, b, N, chf, u).first;
   } else {
       // Optimized version with coefficient caching (default)
       integratedV = COSMethod::invertCDF_Optimized(a, b, N, chf, u).first;  
   }
   
   // Safeguard: ensure integrated variance is positive
   integratedV = std::max(0.0, integratedV);
   
   // ========================================================================
   // STEP 3: Sample ln X(t+Δ) from Gaussian distribution using Equation (11)
   // ========================================================================
   // Conditional on V(t+Δ) and ∫V(u)du, ln X(t+Δ) is Gaussian
   //
   // From Andersen's Equation (11):
   // ln X(t+Δ) = ln X(t) + (ρ/ε)(V(t+Δ) - V(t) - κθΔ) 
   //                     + (κρ/ε - 1/2)∫V(u)du 
   //                     + √(1-ρ²)∫√V(u)dW(u)
   //
   // Rearranging:
   // ln X(t+Δ) = ln X(t) + (ρ/ε)(V(t+Δ) - V(t) - κθΔ + κ∫V) - (1/2)∫V
   //                     + √((1-ρ²)∫V)·Z
   //
   // where Z ~ N(0,1) is INDEPENDENT of W_V
   
   // Generate ONE independent standard normal using Moro's algorithm
   // NOTE: This Z is INDEPENDENT - we do NOT use generateCorrelatedNormals()
   //       because V is sampled exactly (no Brownian motion for V in BK scheme)
   std::uniform_real_distribution<double> uniform_Z(0.0, 1.0);
   double U_Z = uniform_Z(_randomEngine);
   double Z = Utils::inverseNormalCDF(U_Z);  // Independent N(0,1) for diffusion term
   
   double X_t = std::log(assetPrice);  // Current log-price
   
   // Drift term: r·Δ - (1/2)∫V
   // double drift = r * dt - 0.5 * integratedV;
   // double drift = X_t + (rho / sigma_v) * (V_next - V_t - kappa * sigma_v * dt);
   // double correlation_term = (kappa * rho / sigma_v - 1/2) * (V_next - V_t - kappa * vbar * dt + kappa * integratedV);
   // Correlation term: (ρ/ε)·(V(t+Δ) - V(t) - κθΔ + κ∫V)
   // This captures the correlation between S and V through equation (10)
   // double correlation_term = (rho / sigma_v) * (V_next - V_t - kappa * vbar * dt + kappa * integratedV);
   
   // Independent Brownian motion contribution: √((1-ρ²)∫V)·Z
   // This is the √(1-ρ²)∫√V(u)dW(u) term from equation (11)
   // double diffusion = std::sqrt(std::max(0.0, (1.0 - rho * rho) * integratedV)) * Z;
   
   // Update log-price using Equation (11)
   // double X_next = X_t + drift + correlation_term + diffusion;
   
   /// CAUTION: Add drift term for option pricing; remove drift for the "historical simulation" wrt Andersen.
    double drift = r * dt;
    double term1 = (rho/sigma_v) * (V_next - V_t - kappa * vbar * dt);
    double term2 = (kappa * rho / sigma_v - 0.5) * integratedV;
    double term3 = std::sqrt(1 - rho * rho) * std::sqrt(integratedV) * Z;
    double X_next = X_t + drift + term1 + term2 + term3;
    // Convert back to price
    double S_next = std::exp(X_next);
    
    return {S_next, V_next};
}

// ============================================================================
// TG Path Simulator 2D (Truncated Gaussian Scheme)
// ============================================================================

TGPathSimulator2D::TGPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, 
                                     size_t randomSeed, double gamma1, double gamma2)
    : PathSimulator2D(timeSteps, model, randomSeed), _gamma1(gamma1), _gamma2(gamma2)
{
    // Validate γ₁ + γ₂ should typically equal 1 for proper approximation
    if (std::abs(_gamma1 + _gamma2 - 1.0) > 1e-10) {
        // Warning: non-standard choice, but allow flexibility
    }
}

TGPathSimulator2D::~TGPathSimulator2D()
{
    // Base class destructor handles cleanup
}

std::pair<double, double> TGPathSimulator2D::nextStep(size_t timeIndex, double assetPrice, double variance) const
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
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("TGPathSimulator2D requires HestonModel");
    }
    
    // ========================================================================
    // STEP 1: Generate V(t+Δ) using TG scheme (Andersen Eq. 13)
    // ========================================================================
    
    // Generate uniform for variance: U_V ~ Uniform(0,1)
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double U_V = uniform(_randomEngine);
    
    // Transform to standard normal using Moro: Z_V = Φ^(-1)(U_V)
    double Z_V = Utils::inverseNormalCDF(U_V);
    
    // Apply TG scheme to get V(t+Δ)
    double V_next = stepVarianceTG(variance, dt, Z_V);

    // ========================================================================
    // STEP 2-3: Generate INDEPENDENT Z for X discretization (Eq. 33)
    // ========================================================================
    
    // Draw independent uniform U ~ Uniform(0,1), independent of U_V
    double U = uniform(_randomEngine);
    
    // Inverse CDF: Z = Φ^(-1)(U) using Moro's algorithm
    double Z = Utils::inverseNormalCDF(U);
    
    // ========================================================================
    // STEP 4: Compute X(t+Δ) using correlation-preserving scheme (Eq. 33)
    // ========================================================================
    // Notation: X(t) = ln(S(t)) is the log-price
    // The correlation between X and V is preserved through K₁ and K₂ coefficients
    
    double X_t = std::log(assetPrice);                // X(t) = ln(S(t))
    double X_next = PathSimulator2D::stepLogPriceEq33(X_t, variance, V_next, dt, Z, _gamma1, _gamma2);
    double S_next = std::exp(X_next);                 // S(t+Δ) = exp(X(t+Δ))
    
    return {S_next, V_next};
}

double TGPathSimulator2D::stepVarianceTG(double V_t, double dt, double Z_V) const
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
     * Q: For what γ's this scheme coincides with Euler with full truncation, similar to Eq. 6-7? γ1 = 1.0, γ2 = 0.0 - gives Euler-like scheme - test
     * Q: Can we do better? set γ1 = γ2 = 0.5 - test
     * Q: What is the best way to choose γ's? - Moment Matching is proposed see QE below - test; compare all 3 cases.
     */
    
    const HestonModel* hestonPtr = static_cast<const HestonModel*>(_modelPtr);
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
    double V_next = std::max(V_t + drift + diffusion, 0.0);    // Eq. 13
    
    return V_next;
}

// ============================================================================
// QE Path Simulator 2D 
// ============================================================================

QEPathSimulator2D::QEPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, 
                                     size_t randomSeed, double psi_c, double gamma1, double gamma2)
    : PathSimulator2D(timeSteps, model, randomSeed), _psi_c(psi_c), _gamma1(gamma1), _gamma2(gamma2)
{
    // Validate γ₁ + γ₂ should typically equal 1 for proper approximation
    if (std::abs(_gamma1 + _gamma2 - 1.0) > 1e-10) {
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
    const HestonModel* hestonPtr = dynamic_cast<const HestonModel*>(_modelPtr);
    if (!hestonPtr) {
        throw std::runtime_error("QEPathSimulator2D requires HestonModel");
    }
    
    // ========================================================================
    // STEP 1: Generate V(t+Δ) using QE scheme (Andersen Eq. 23/26)
    // ========================================================================
    
    // Generate uniform for variance: U_V ~ Uniform(0,1)
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double U_V = uniform(_randomEngine);
    
    // Transform to standard normal using Moro: Z_V = Φ^(-1)(U_V)
    double Z_V = Utils::inverseNormalCDF(U_V);
    
    // Apply QE scheme to get V(t+Δ)
    double V_next = stepVarianceQE(variance, dt, Z_V);
    
    // ========================================================================
    // STEP 2-3: Generate INDEPENDENT Z for X discretization (Eq. 33)
    // ========================================================================
    
    // Draw independent uniform U ~ Uniform(0,1), independent of U_V
    double U = uniform(_randomEngine);
    
    // Inverse CDF: Z = Φ^(-1)(U) using Moro's algorithm
    double Z = Utils::inverseNormalCDF(U);
    
    // ========================================================================
    // STEP 4: Compute X(t+Δ) using correlation-preserving scheme (Eq. 33)
    // ========================================================================
    // Notation: X(t) = ln(S(t)) is the log-price
    // The correlation between X and V is preserved through K1 and K2 coefficients
    
    double X_t = std::log(assetPrice);                // X(t) = ln(S(t))
    double X_next = PathSimulator2D::stepLogPriceEq33(X_t, variance, V_next, dt, Z, _gamma1, _gamma2);
    double S_next = std::exp(X_next);                 // S(t+Δ) = exp(X(t+Δ))
    
    return {S_next, V_next};
}

double QEPathSimulator2D::stepVarianceQE(double V_t, double dt, double Z_V) const
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
     */
    
    const HestonModel* hestonPtr = static_cast<const HestonModel*>(_modelPtr);
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
    double s2 = V_t * sigma_v * sigma_v / kappa * (exp_kappa_dt - exp_kappa_dt * exp_kappa_dt)
              + vbar * sigma_v * sigma_v / (2.0 * kappa) * (1.0 - exp_kappa_dt) * (1.0 - exp_kappa_dt);       
    
    // Scaled variance: ψ = s²/m²
    double psi = s2 / (m * m);   // Eq. 19
    
    // ========================================================================
    // 3.2.3 Switching rule 
    // ========================================================================
    
    if (psi <= _psi_c) {
        // --------------------------------------------------------------------
        // Quadratic scheme (ψ ≤ ψ_c): V ~ a(b + Z_V)²
        // --------------------------------------------------------------------
        double b_squared = 2.0 / psi - 1.0 + std::sqrt(2.0 / psi) * std::sqrt(2.0 / psi - 1.0); // Eq. 27
        double a = m / (1.0 + b_squared);   // Eq. 28
        double b = std::sqrt(b_squared);    
        
        double V_next = a * (b + Z_V) * (b + Z_V);
        return std::max(V_next, 0.0); // Ensure non-negativity
    } else {
        // --------------------------------------------------------------------
        // Exponential scheme (ψ > ψ_c): V ~ exponential/point mass mixture
        // --------------------------------------------------------------------
        double p = (psi - 1.0) / (psi + 1.0);
        double beta = (1.0 - p) / m;
        
        // Generate uniform for mixture
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double U_V = uniform(_randomEngine);
        
        if (U_V <= p) { // Eq. 25
            // Point mass at 0
            return 0.0;
        } else {
            // Exponential: V = (1/β)·ln((1-p)/(1-U))
            double V_next = std::log((1.0 - p) / (1.0 - U_V)) / beta;
            return V_next;
        }
    }
}


