//
// Created by Patrik  Liba on 27/11/2025.
//

#ifndef CPPFM_PATHSIMULATOR2D_H
#define CPPFM_PATHSIMULATOR2D_H

#include "Model.h"
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
 */

/** 
 * Monte Carlo Simulation Methods for Heston Model
 *   |   
 *   ├─ Scheme 1: Euler-Euler [DONE]
 *   ├─ V: Euler with truncation (or reflection) (eq. 6-7)
 *   └─ X: Euler (eq. 6)
 *
 *   ├─ Scheme 2: IM-IJK (Kahl-Jäckel) 
 *   ├─ V: Implicit Milstein (eq. 9)
 *   └─ X: IJK (eq. 8)
 *
 *   ├─ Scheme 3: Broadie-Kaya (BK) [DONE]
 *   ├─ V: Exact non-central χ² sampling
 *   └─ X: Exact via COS and Newton (for CDF inversion) 
 *
 *   ├─ Scheme 4: TG (Andersen) [DONE]
 *   ├─ V: Truncated Gaussian (eq. 13)
 *   └─ X: Correlation-preserving (eq. 33)
 *
 *   └─ Scheme 5: QE (Andersen) [DONE]
 *   ├─ V: Quadratic-Exponential (eq. 23/26)
 *   └─ X: Correlation-preserving (eq. 33)
 *
 * 
 
 */


// Abstract class for different discretization schemes for 2D models (e.g., Heston)
class PathSimulator2D
{
public: 
    PathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, size_t randomSeed);
    virtual ~PathSimulator2D();
    virtual std::pair<double, double> nextStep(size_t timeIndex, double assetPrice, double variance) const = 0;
    
    // Generate full paths for both S and V
    std::pair<std::vector<double>, std::vector<double>> paths() const;

    /// TODO: Can tidy with all those params declared here?

protected:
    bool timeStepsSanityCheck() const;
    
    // Generate correlated normal random variables with correlation from model
    std::pair<double, double> generateCorrelatedNormals() const;
    
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
    double stepLogPriceEq33(double X_t, double V_t, double V_next, double dt, double Z,
                           double gamma1, double gamma2) const;

    const std::vector<double>& _timeSteps;
    const Model2D* _modelPtr; // Pointer to base class Model2D - for polymorphic 2D models
    size_t _randomSeed;
    mutable std::default_random_engine _randomEngine; 
};

/** EULER as a pair for (X, V)
 * Full Truncation Scheme for Heston Model (Equations 6-7)
 * 
 * Notation: X(t) = ln(S(t)) is the log-price
 *
 * X̂(t + Δ) = X̂(t) - (1/2)V̂(t)⁺Δ + √(V̂(t)⁺)·Z_X·√Δ
 * V̂(t + Δ) = V̂(t) + κ(θ - V̂(t)⁺)Δ + ε·√(V̂(t)⁺)·Z_V·√Δ
 *
 * where x⁺ = max(x, 0), and Z_X, Z_V are correlated standard normals with correlation ρ.
 *
 * This scheme prevents negative variance by using V⁺ in all terms.
 * When V goes negative, it reverts back deterministically via drift term κθΔ.
 */
class EulerPathSimulator2D : public PathSimulator2D
{
public:
    EulerPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, size_t randomSeed);
    /**
     * Perform one time step using Full Truncation discretization
     * @param timeIndex Current time index in _timeSteps
     * @param assetPrice Current asset price S(t)
     * @param variance Current variance V(t)
     * @return Pair of (S(t+Δ), V(t+Δ))
     */
     std::pair<double, double> nextStep(size_t timeIndex, double assetPrice, double variance) const override;

};       
  

class MilsteinPathSimulator2D : public PathSimulator2D
{

};

// BK as a pair for (X, V)
/// DONE: Implement "true" BK scheme matching Andersen's paper
/// DONE: Implement BK with a choice between `your Newton` and `Newton Optimized` method - test compare performance

enum class NewtonMethod { Original,     // my implementation; 
                          Optimized     // Optimized version - cached coefficients (as default)
};

class BKPathSimulator2D : public PathSimulator2D
{
public:
    BKPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, size_t randomSeed,
                      NewtonMethod newtonMethod = NewtonMethod::Optimized);
    std::pair<double, double> nextStep(size_t timeIndex, double assetPrice, double variance) const override;

protected:
    NewtonMethod _newtonMethod;
};

// For X: (Equation 33) For V: Truncated Gaussian (TG) (Equation 13)
class TGPathSimulator2D : public PathSimulator2D
{
public:
    TGPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, size_t randomSeed,
                      double gamma1 = 0.5, 
                      double gamma2 = 0.5);
    virtual ~TGPathSimulator2D() override;
    
    /**
     * Perform one time step using TG scheme for V and Equation 33 for X
     * @param timeIndex Current time index in _timeSteps
     * @param assetPrice Current asset price S(t)
     * @param variance Current variance V(t)
     * @return Pair of (S(t+Δ), V(t+Δ))
     */
    std::pair<double, double> nextStep(size_t timeIndex, double assetPrice, double variance) const override;
    
protected:
    double _gamma1;         // γ₁: weight for V(t) in integral approximation
    double _gamma2;         // γ₂: weight for V(t+Δt) in integral approximation
    
private:
    /**
     * Compute V(t+Δ) using Truncated Gaussian scheme (Andersen Eq. 13)
     * @param V_t Current variance V(t)
     * @param dt Time step Δt
     * @param Z_V Standard normal for variance
     * @return V(t+Δ)
     */
    double stepVarianceTG(double V_t, double dt, double Z_V) const;
};

// For X: (Equation 33) For V: Quadratic-Exponential (QE) (Equation 23/26)
class QEPathSimulator2D : public PathSimulator2D
{
public:
    QEPathSimulator2D(const std::vector<double>& timeSteps, const Model2D& model, size_t randomSeed,
                      double psi_c = 1.5,     // for the switching rule
                      double gamma1 = 0.5, 
                      double gamma2 = 0.5);
    virtual ~QEPathSimulator2D() override;
    
    /**
     * Perform one time step using QE scheme for V and Equation 33 for X
     * @param timeIndex Current time index in _timeSteps
     * @param assetPrice Current asset price S(t)
     * @param variance Current variance V(t)
     * @return Pair of (S(t+Δ), V(t+Δ))
     */
    std::pair<double, double> nextStep(size_t timeIndex, double assetPrice, double variance) const override;
    
protected:
    double _psi_c;          // ψ_c: Threshold parameter for switching between quadratic/exponential
    double _gamma1;         // γ₁: weight for V(t) in integral approximation
    double _gamma2;         // γ₂: weight for V(t+Δt) in integral approximation
    
private:
    /**
     * Compute V(t+Δ) using Quadratic-Exponential scheme (Andersen Eq. 23/26)
     * @param V_t Current variance V(t)
     * @param dt Time step Δt
     * @param Z_V Standard normal for variance
     * @return V(t+Δ)
     */
    double stepVarianceQE(double V_t, double dt, double Z_V) const;
};





#endif //CPPFM_PATHSIMULATOR2D_H