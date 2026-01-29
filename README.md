# C++ Option Pricing Library

## Contents

| Component | Description |
|-----------|-------------|
| Models | Black-Scholes, Dupire local vol, Heston stochastic vol |
| Monte Carlo | 4 Heston schemes from Andersen (2008), SLV |
| Fourier | COS method, CDF/PDF recovery, Newton inversion |
| Finite Difference | Theta-method PDE solver, automatic Greeks |
| Market Data | Vol surface builder, Dupire local vol, discount curves |

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
│   ├── DiscountCurve.h      Flat/forward rate curves
│   └── FinancialInstrument.h European/American options
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
```

## Papers

- Andersen (2008) - Efficient Simulation of the Heston Process
- Broadie & Kaya (2006) - Exact Simulation of Stochastic Volatility
- Fang & Oosterlee (2008) - A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions
- van der Stoep, Grzelak & Oosterlee (2013) - The Heston Stochastic-Local Volatility Model
- Albrecher et al. (2007) - The Little Heston Trap
