# Onboarding

Quick guide to get productive with the C++ option pricing library.

## Quick Start

### Dependencies

- C++20
- CMake 3.18+
- Google Test 1.14.0 (fetched automatically)
- Boost, libomp, TBB (Homebrew on macOS)

```bash
# macOS install
brew install boost libomp tbb
```

Paths expected at `/opt/homebrew/opt/{boost,libomp,tbb}`.

### Build

```bash
cmake -S . -B build && cmake --build build
```

### Run Tests

```bash
cd build && ctest
# or single test
./build/test_heston_model
```

## Architecture

```
ModelBase
├── Model1D
│   ├── BlackScholesModel  (constant vol)
│   └── DupireModel        (local vol from IV surface)
└── Model2D
    └── HestonModel        (stochastic vol)

PathSimulator (1D)
├── EulerPathSimulator
└── MilsteinPathSimulator

PathSimulator2D (Heston schemes - Andersen 2008)
├── EulerPathSimulator2D     (Eq. 6-7, full truncation)
├── TGPathSimulator2D        (Eq. 13 + 33)
├── QEPathSimulator2D        (Eq. 23/26 + 33)
├── BKTGPathSimulator2D      (TG variance + Eq. 11)
├── BKQEPathSimulator2D      (QE variance + Eq. 11)
├── BKExactPathSimulator2D   (exact CIR + exact ∫V)
└── HestonSLVPathSimulator2D (van der Stoep 2013)

Pricers
├── BlackScholesPricer  (closed-form)
├── MonteCarloPricer    (any PathSimulator)
├── FDPricer            (theta-method PDE)
└── COSPricer           (Fang & Oosterlee 2008)
```


## File Map

| File | Purpose |
|------|---------|
| `Model.h/cpp` | Model hierarchy (BS, Dupire, Heston) |
| `PathSimulator.h/cpp` | 1D Euler/Milstein |
| `PathSimulator2D.h/cpp` | Heston schemes (TG, QE, BK family, SLV) |
| `Pricer.h/cpp` | BS, MC, FD pricers |
| `VolatilitySurface.h/cpp` | IV surface + Dupire local vol |
| `DiscountCurve.h/cpp` | Flat/forward rate curves |
| `Utils.h/cpp` | COS method, CIR sampler, HestonLocalVol |
| `PDEs/` | Grid, PDE, ThetaMethodSolver, BoundaryConditions |
| `BlackScholesFormulas.h/cpp` | Closed-form BS formulas |
| `InterpolationSchemes.h/cpp` | Linear, cubic spline interpolation |

## Heston Discretization Schemes

The library's main complexity is in `PathSimulator2D`. All schemes reference Andersen (2008).

**Variance (V) discretization:**
- Euler: Full truncation (Eq. 6-7) - use V⁺ = max(V, 0) everywhere
- TG: Truncated Gaussian (Eq. 13) - moment matching
- QE: Quadratic-Exponential (Eq. 23/26) - switches based on ψ threshold

**Log-price (X) discretization:**
- Euler: Standard (Eq. 6)
- Eq. 33: Correlation-preserving, uses γ₁V(t) + γ₂V(t+Δ) approximation for ∫V
- Eq. 11: Exact SDE representation (Broadie-Kaya)

**Scheme combinations (PathSimulator2D.h:49-86):**
```
X scheme    V scheme       Class
────────────────────────────────────────
Euler       Euler          EulerPathSimulator2D
Eq.33       TG (Eq.13)     TGPathSimulator2D
Eq.33       QE (Eq.23/26)  QEPathSimulator2D
Eq.11       TG + approx    BKTGPathSimulator2D
Eq.11       QE + approx    BKQEPathSimulator2D
Eq.11       Exact CIR      BKExactPathSimulator2D
```

QE is the workhorse - handles Feller violation gracefully while being much faster than BKExact.

## Tests

Location: `/tests/`
Framework: Google Test

Good starting points:
- `test_heston_model.cpp` - Feller condition impact, scheme comparison
- `test_black_scholes.cpp` - basic BS pricer tests
- `test_pde_solver.cpp` - FD convergence tests
- `test_cos_method.cpp` - Fourier pricing tests

Run specific test:
```bash
./build/test_heston_model
```
