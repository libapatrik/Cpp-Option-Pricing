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

## Directory Structure

```
CppFM/
├── include/cppfm/           # Public headers
│   ├── models/              # Model hierarchy (BS, Dupire, Heston)
│   ├── simulators/          # Path simulation (1D Euler/Milstein, 2D Heston)
│   ├── pricers/             # Closed-form, MC, FD pricers
│   ├── pde/                 # Finite difference framework
│   ├── market/              # Vol surfaces, discount curves, instruments
│   ├── montecarlo/          # MC engine + variance reduction (WIP)
│   └── utils/               # COS method, interpolation, tridiagonal solver
├── src/                     # Implementation files (mirrors include/)
│   └── bindings/            # Python bindings (pybind11)
├── tests/                   # Google Test suite
├── python/examples/         # Jupyter notebooks
└── docs/                    # ONBOARDING.md, OPTIMIZATION.md
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

| Module | Headers | Implementation |
|--------|---------|----------------|
| Models | `include/cppfm/models/Model.h` | `src/models/Model.cpp` |
| Path Sim 1D | `include/cppfm/simulators/PathSimulator.h` | `src/simulators/PathSimulator.cpp` |
| Path Sim 2D | `include/cppfm/simulators/PathSimulator2D.h` | `src/simulators/PathSimulator2D.cpp` |
| Pricers | `include/cppfm/pricers/Pricer.h` | `src/pricers/Pricer.cpp` |
| Vol Surface | `include/cppfm/market/VolatilitySurface.h` | `src/market/VolatilitySurface.cpp` |
| Discount Curve | `include/cppfm/market/DiscountCurve.h` | `src/market/DiscountCurve.cpp` |
| Utils/COS | `include/cppfm/utils/Utils.h` | `src/utils/Utils.cpp` |
| PDE | `include/cppfm/pde/*.h` | `src/pde/*.cpp` |
| Interpolation | `include/cppfm/utils/InterpolationSchemes.h` | `src/utils/InterpolationSchemes.cpp` |

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

**Scheme combinations:**
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

Location: `tests/`
Framework: Google Test

Good starting points:
- `test_heston_model.cpp` - Feller condition impact, scheme comparison
- `test_heston_scheme.cpp` - scheme convergence and accuracy
- `test_black_scholes.cpp` - basic BS pricer tests
- `test_pde_solver.cpp` - FD convergence tests
- `test_cos_method.cpp` - Fourier pricing tests
- `test_heston_slv.cpp` - SLV calibration tests

Run specific test:
```bash
./build/test_heston_model
```

## Python Bindings

Python module exposes core functionality via pybind11.

Build:
```bash
cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON
cmake --build build
```

Example notebooks in `python/examples/`.
