# Onboarding

## Repository Structure

```
Cpp-Option-Pricing/
├── include/cppfm/        # Public headers
│   ├── models/           # BS, Dupire, Heston
│   ├── simulators/       # Path simulation (1D, 2D Heston)
│   ├── pricers/          # MC, FD, COS, closed-form
│   ├── pde/              # Finite difference framework
│   ├── market/           # Vol surfaces, discount curves
│   └── utils/            # COS internals, interpolation, tridiagonal
├── src/                  # Implementation (mirrors include/)
│   └── bindings/         # pybind11 Python bindings
├── tests/                # Google Test suite
├── python/               # Python package + Jupyter examples
├── .github/workflows/    # CI
├── CMakeLists.txt
├── pyproject.toml        # pip install config
└── README.md
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

PathSimulator2D (Heston - Andersen 2008)
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

| Module | Header | Implementation |
|--------|--------|----------------|
| Models | `include/cppfm/models/Model.h` | `src/models/Model.cpp` |
| Path Sim 1D | `include/cppfm/simulators/PathSimulator.h` | `src/simulators/PathSimulator.cpp` |
| Path Sim 2D | `include/cppfm/simulators/PathSimulator2D.h` | `src/simulators/PathSimulator2D.cpp` |
| Pricers | `include/cppfm/pricers/Pricer.h` | `src/pricers/Pricer.cpp` |
| Vol Surface | `include/cppfm/market/VolatilitySurface.h` | `src/market/VolatilitySurface.cpp` |
| Discount Curve | `include/cppfm/market/DiscountCurve.h` | `src/market/DiscountCurve.cpp` |
| Utils/COS | `include/cppfm/utils/Utils.h` | `src/utils/Utils.cpp` |
| PDE | `include/cppfm/pde/*.h` | `src/pde/*.cpp` |

## Heston Discretization

All schemes reference Andersen (2008).

**Variance (V):**
- Euler: Full truncation (Eq. 6-7)
- TG: Truncated Gaussian (Eq. 13)
- QE: Quadratic-Exponential (Eq. 23/26)

**Log-price (X):**
- Euler: Standard (Eq. 6)
- Eq. 33: Correlation-preserving
- Eq. 11: Exact (Broadie-Kaya)

**Combinations:**
```
X        V              Class
─────────────────────────────────────
Euler    Euler          EulerPathSimulator2D
Eq.33    TG             TGPathSimulator2D
Eq.33    QE             QEPathSimulator2D
Eq.11    TG + approx    BKTGPathSimulator2D
Eq.11    QE + approx    BKQEPathSimulator2D
Eq.11    Exact CIR      BKExactPathSimulator2D
```

QE handles Feller violation gracefully, much faster than BKExact.

## Python
Then from any Python project use:
```python
import cppfm

# Heston pricing
price = cppfm.heston_call_price(S=100, K=100, T=1.0, r=0.05,
                                 V0=0.04, kappa=2.0, theta=0.04,
                                 sigma=0.3, rho=-0.7)

# Black-Scholes
call = cppfm.bs_call_price(S=100, K=100, T=1.0, r=0.05, sigma=0.2)
```

Example notebooks in `python/examples/`.
