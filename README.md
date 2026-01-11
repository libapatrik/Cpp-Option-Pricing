# Option Pricing Library

## Models

- **Black-Scholes**: constant vol, analytical Greeks
- **Dupire**: local vol $\sigma(S,t)$ from IV surface
- **Heston**: stochastic vol with mean reversion

## Path Simulators

1D schemes (Black-Scholes, Dupire):
- Euler, Milstein

2D Heston schemes (Andersen 2008):
- Euler (full truncation, Eq. 6-7)
- Truncated Gaussian (Eq. 13)
- Quadratic-Exponential (Eq. 23/26)
- Broadie-Kaya variants: BK-TG, BK-QE, BK-Exact
- Eq. 33 log-price discretization

Heston SLV (van der Stoep et al. 2013):
- Leverage function $L^2(S,t) = \sigma^2_{LV} / \mathbb{E}[V|S]$
- Binning for conditional expectation
- Iterative calibration

## Pricers

- **BlackScholesPricer**: closed-form
- **MonteCarloPricer**: any PathSimulator, antithetic sampling
- **FDPricer**: theta-method PDE solver
- **COSPricer**: Fourier-cosine expansion (Fang & Oosterlee 2008)

## Market Data

- **DiscountCurve**: flat curve, instantaneous forward rates
- **VolatilitySurface**: strike/maturity grid, linear/cubic smile interpolation
- **VolatilitySurfaceBuilder**: fluent API

## PDE Framework

- **Grid**: uniform or log-spaced
- **ThetaMethodSolver**: explicit ($\theta=0$), implicit ($\theta=1$), Crank-Nicolson ($\theta=0.5$)
- **BoundaryConditions**: Dirichlet
- **PDEGreeksCalculator**: bump-and-revalue for vega, rho, vanna, volga

## Utilities

- Thomas algorithm (tridiagonal solver)
- COS method: CDF/PDF recovery, CDF inversion via Newton
- CIR exact sampling (non-central $\chi^2$)
- ChF of integrated variance (Broadie & Kaya 2006)
- Heston ChF with Little Heston Trap (Albrecher 2007)
- HestonLocalVol: analytical Dupire via COS
- Interpolation: linear, natural cubic spline
- PCG32 RNG
