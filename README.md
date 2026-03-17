# C++ Option Pricing Library

## Contents

| Component | Description |
|-----------|-------------|
| Models | Black-Scholes, Dupire local vol, Heston stochastic vol |
| Monte Carlo | 4 Heston schemes from Andersen (2008), SLV |
| Fourier | COS method with cumulant-based truncation, CDF inversion |
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
- Broadie-Kaya family (TG, QE, Exact) - exact log-price distribution via Eq. 11, COS-based integrated variance sampling with cumulant bounds (N=48)
- Heston SLV (van der Stoep et al.) - stochastic-local vol hybrid with leverage function

QE is the main choice - handles Feller violation, much faster than exact sampling.

**Finite Difference**: theta-method solver (explicit/implicit/Crank-Nicolson), Dirichlet boundaries, automatic Greeks via bump-and-revalue (delta, gamma, vega, rho, vanna, volga).

**Fourier**: COS method for European pricing with cumulant-based truncation bounds (Fang & Oosterlee 2008 Eq. 23). Analytical Heston cumulants from Le Floc'h (2018). CDF/PDF recovery from characteristic function, Newton-based CDF inversion with coefficient caching for BK exact sampling.

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
├── cos/
│   └── COS.h                COS transforms, HestonCF + cumulants, ChFIntegratedVariance, truncation bounds
└── utils/
    ├── Utils.h              CIR sampler, Thomas algorithm
    └── InterpolationSchemes.h Linear, cubic spline

src/                         Implementation (mirrors include/)
tests/                       Test files
python/                      pybind11 bindings + Jupyter notebooks
docs/                        ONBOARDING.md - further details about the codebase
```

## Tests

```bash
cd build && ctest            # all 22 test suites
./build/test_heston_scheme   # scheme convergence
./build/test_heston_slv      # SLV calibration + control variate
./build/test_cos_method      # Fourier pricing accuracy
./build/test_pde_greeks      # full Greeks suite vs analytical
./build/test_svi_calibration # SVI/SSVI smile fitting + arb checks
```

### Heston MC Discretization Schemes

7 schemes from Andersen (2008), compared head-to-head on 50k paths.

| Scheme | Variance Discretization | Log-Price | Bias vs BK Exact | Feller Robust |
|--------|------------------------|-----------|-------------------|---------------|
| Euler (full truncation) | Eq. 6-7 | Euler | ~25-30% | no |
| Truncated Gaussian | Eq. 13 | Eq. 33 | ~5-10% | partial |
| Quadratic-Exponential | Eq. 23/26 | Eq. 33 | <1% | yes |
| BK/TG | CIR exact | Eq. 11 (COS) | <1% | yes |
| BK/QE | CIR exact | Eq. 11 (COS) | <1% | yes |
| BK Exact | CIR exact + Newton CDF | Eq. 11 (COS) | reference | yes |
| SLV (van der Stoep) | QE + leverage L(S,t) | adjusted drift | N/A | yes |

Feller condition test runs 4 parameter regimes (violated, satisfied, boundary, equality) across all schemes - QE stays stable where Euler blows up.

### Heston Stochastic-Local Volatility

Full SLV pipeline validated against van der Stoep et al. (2013).

| Test | Description | Result |
|------|-------------|--------|
| Van der Stoep Table 1 params | COS surface -> Dupire -> leverage calibration -> MC repricing | Max error 39 bp |
| Heston-consistent surface | Market = Heston IVs, so L(S,t) should be ~1 | Max error 20 bp |
| Control variate (Heston as CV) | Same Brownians, COS analytical baseline | **120x variance reduction** |
| CV across strikes (K=0.90-1.10) | Per-strike SE reduction | 80x-135x, avg 110x |
| External surface fitting | U-shaped smile (non-Heston) via leverage | Fits arbitrary smiles |
| Mixing factor eta | eta=0 (pure LV) to eta=1 (full SLV) | eta=1 within 44 bp of Heston |

### Calibration

| Test | Method | Result |
|------|--------|--------|
| SVI round-trip | 5 params -> 9 strikes -> LM -> recovery | RMSE < 1e-6, all params within 1e-6 |
| SSVI round-trip | 3 global params -> 4 maturities x 9 strikes -> LM | RMSE < 1e-4 (sub-bp) |
| SSVI butterfly arb | g(k) >= 0 across 41 points (Gatheral & Jacquier 2014) | arb-free verified |
| SSVI calendar arb | total variance non-decreasing in theta | monotonicity enforced (PAVA) |
| SSVI analytical Jacobian | dw/d{rho,eta,gamma} vs finite differences | all within 1e-4 |
| Heston round-trip | 4 mat x 7 strikes -> LM (tol=1e-10) -> recovery | RMSE < 1e-4 (sub-bp) |
| Heston grid search | 3^5 = 243 starting points -> parallel -> LM refinement | >1.5x speedup, same optimum |

### COS Method & Cumulants

| Test | Description | Result |
|------|-------------|--------|
| PDF recovery vs N | Replicates F&O 2008 Table 1, N={4..256} | Exponential convergence |
| Heston cumulants | Le Floc'h (2018) analytical vs numerical FD | c1, c2 within 1e-6 |
| Cumulant-based bounds | Auto [a,b] from c1, c2, c4 (F&O Eq. 23) | Matches sigma-hint pricing within 0.02 |
| CDF inversion benchmark | Newton iterations + timing for N={64..512} | Coefficient caching + bisection init |

### PDE Solver & Greeks

| Test | Description | Result |
|------|-------------|--------|
| 8 Greeks vs BS analytical | delta, gamma, vega, theta, rho, vanna, volga (call + put) | All 16 within 1e-3 |
| Crank-Nicolson convergence | L2 error ratio on grid doubling | >3.0 (second-order) |
| Heat equation (3 schemes) | Explicit, Implicit, CN vs analytical | CN error < 0.001 |
| BS European call (PDE vs formula) | Price + delta accuracy | <0.5% relative error |

### Dupire Local Volatility

| Test | Description | Result |
|------|-------------|--------|
| Flat vol = BS (FD + MC) | Dupire with constant surface vs BS analytical | FD < 0.10, MC < 0.30 |
| Smile affects prices | Non-flat surface diverges from BS | verified for OTM puts |
| FD grid convergence | Error shrinks with refinement (100->200->400) | second-order convergence |
| Put-call parity | C - P = S - Ke^{-rT} via FD and MC | FD < 0.10, MC < 0.30 |
| Milstein vs Euler MC | Both 1D discretizations agree | within 0.50 |

## Papers

- Andersen (2008) - Efficient Simulation of the Heston Process
- Broadie & Kaya (2006) - Exact Simulation of Stochastic Volatility
- Fang & Oosterlee (2008) - A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions
- Le Floc'h (2018) - Analytical Heston Cumulants (arXiv:2005.13248, corrects F&O Table 2)
- van der Stoep, Grzelak & Oosterlee (2013) - The Heston Stochastic-Local Volatility Model
- Albrecher et al. (2007) - The Little Heston Trap
