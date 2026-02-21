# C++ Option Pricing Library

## Contents

| Component | Description |
|-----------|-------------|
| Models | Black-Scholes, Dupire local vol, Heston stochastic vol |
| Monte Carlo | 4 Heston schemes from Andersen (2008), SLV |
| Fourier | COS method, CDF/PDF recovery, Newton inversion |
| Finite Difference | Theta-method PDE solver, automatic Greeks |
| Market Data | Vol surface builder, Dupire local vol, discount curves |
| Calibration | SVI/SSVI smile fitting, Levenberg-Marquardt, arbitrage checks |
| Linear Algebra | Cholesky, Householder QR, Golub-Kahan SVD, least-squares |

## Models

```cpp
// constant vol
BlackScholesModel model(S0, r, sigma);

// local vol from IV surface
DupireModel model(S0, r, volSurface);

// stochastic vol
HestonModel model(S0, V0, r, kappa, theta, sigma, rho);
```

## Pricing Methods

**Monte Carlo** with 4 Heston discretization schemes:

- Euler (full truncation, Eq. 6-7)
- Truncated Gaussian (Eq. 13), Quadratic-Exponential (Eq. 23/26)
- Broadie-Kaya family (TG, QE, Exact) - exact log-price distribution via Eq. 11
- Heston SLV (van der Stoep et al.) - stochastic-local vol hybrid with leverage function

QE is the main choice - handles Feller violation, much faster than exact sampling.

**Finite Difference**: theta-method solver (explicit/implicit/Crank-Nicolson), Dirichlet boundaries, automatic Greeks via bump-and-revalue (delta, gamma, vega, rho, vanna, volga).

**Fourier**: COS method for fast European pricing, CDF/PDF recovery from characteristic function, Newton-based CDF inversion for sampling.

## Volatility Surface

```cpp
auto surface = VolatilitySurfaceBuilder()
    .setSpot(100.0)
    .setRate(0.05)
    .addSlice(0.25, {80, 90, 100, 110, 120}, {0.25, 0.22, 0.20, 0.22, 0.25})
    .addSlice(0.5, {80, 90, 100, 110, 120}, {0.24, 0.21, 0.19, 0.21, 0.24})
    .build();

// Dupire local vol
double localVol = surface.localVol(S, t);
```

Linear or cubic spline interpolation in strike, linear in time.

## Calibration

Market IV quotes to smooth parametric surface via SVI/SSVI.

```cpp
// fit single smile slice with SVI (5 params: a, b, rho, m, sigma)
SviParams guess(0.04, 0.4, -0.5, 0.0, 0.1);
auto result = calibrateSmile(&guess, strikes, marketVols, forward, T);

// SSVI surface - 3 global params (rho, eta, gamma), per-maturity theta
// enforces no-butterfly-arbitrage (Gatheral & Jacquier 2014, Theorem 4.3)
auto ssviResult = calibrateSsvi(strikesPerMaturity, volsPerMaturity, forwards, maturities);
auto surface = ssviResult.buildSurface(strikeGrid, forwards, discountCurve);
```

Optimizer pipeline: LM minimizes `||w_model(k) - w_market||` in total variance space, with parameter bounds and analytical Jacobians. Roger Lee moment formula enforced on wings.

## Linear Algebra & Optimization

Hand-rolled numerical linear algebra (no Eigen dependency):

- **MatrixOps** - matrix-vector/matrix-matrix multiply, cache-blocked transpose, `A^T*A` without forming `A^T`
- **Cholesky** - SPD decomposition with forward/back substitution
- **QR** - Householder reflections, implicit Q storage (no explicit Q formed)
- **SVD** - Golub-Kahan bidiagonalization + implicit QR shifts, truncated solve for rank-deficient systems
- **LeastSquares** - QR or SVD backend, wraps the above

Optimizers:

- **Levenberg-Marquardt** - nonlinear least-squares with bound constraints, numerical or analytical Jacobian
- **L-BFGS** - quasi-Newton for smooth unconstrained problems (Nocedal & Wright Alg. 7.4)
- **Gradient Descent** - steepest descent with Armijo backtracking

## Building

Needs C++20, CMake 3.18+, Boost, OpenMP, TBB.

```bash
# macOS
brew install boost libomp tbb

# build
cmake -S . -B build && cmake --build build

# test
cd build && ctest
```

## Code Structure

```
include/cppfm/
├── models/
│   └── Model.h              ModelBase, Model1D, Model2D, BS, Dupire, Heston
├── simulators/
│   ├── PathSimulator.h      Euler, Milstein (1D)
│   └── PathSimulator2D.h    All Heston schemes (2D)
├── pricers/
│   ├── Pricer.h             MonteCarloPricer, FDPricer, COSPricer
│   └── BlackScholesFormulas.h  Closed-form BS
├── pde/
│   ├── Grid.h               Uniform/log-spaced grids
│   ├── PDE.h                PDE coefficient definition
│   ├── Solver.h             ThetaMethodSolver
│   ├── BoundaryConditions.h Dirichlet BC
│   └── PDEGreeksCalculator.h Bump-and-revalue Greeks
├── market/
│   ├── VolatilitySurface.h  IV surface + Dupire local vol
│   ├── SviCalibration.h     SVI/SSVI params, smile fitting, arbitrage checks
│   ├── DiscountCurve.h      Flat/forward rate curves
│   └── FinancialInstrument.h European/American options
├── calibration/
│   ├── LinearAlgebra.h      MatrixOps, Cholesky, QR, SVD, LeastSquares
│   └── Optimizer.h          LM, L-BFGS, GradientDescent
└── utils/
    ├── Utils.h              COS method, CIR sampler, HestonCF, Thomas algorithm
    └── InterpolationSchemes.h Linear, cubic spline

src/                         Implementation (mirrors include/)
tests/                       Test files
python/                      pybind11 bindings + Jupyter notebooks
docs/                        ONBOARDING.md - further details about the codebase
```

## Tests

```bash
./build/test_heston_scheme   # scheme convergence
./build/test_cos_method      # Fourier pricing accuracy
./build/test_pde_solver      # FD stability and Greeks
./build/test_heston_slv      # SLV calibration
./build/test_linear_algebra  # decompositions, least-squares
./build/test_optimizer       # LM, L-BFGS convergence
./build/test_svi_calibration # SVI/SSVI smile fitting
```

## Papers

- Andersen (2008) - Efficient Simulation of the Heston Process
- Broadie & Kaya (2006) - Exact Simulation of Stochastic Volatility
- Fang & Oosterlee (2008) - A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions
- van der Stoep, Grzelak & Oosterlee (2013) - The Heston Stochastic-Local Volatility Model
- Albrecher et al. (2007) - The Little Heston Trap
