"""
CppFM Python Bindings
=====================

Python interface to the CppFM quantitative finance library.

Provides access to:
- Discount curves (FlatDiscountCurve)
- Volatility surfaces with Dupire local vol
- Heston stochastic volatility model
- Heston SLV simulator (van der Stoep et al. 2013)
- Black-Scholes formulas

Quick Start
-----------
>>> import cppfm
>>> import numpy as np

>>> # Create discount curve (5% flat rate)
>>> curve = cppfm.FlatDiscountCurve(0.05)

>>> # Create volatility surface from market data
>>> strikes = [80, 90, 100, 110, 120]
>>> maturities = [0.25, 0.5, 1.0]
>>> ivs = [[0.25, 0.22, 0.20, 0.22, 0.25],
...        [0.24, 0.21, 0.19, 0.21, 0.24],
...        [0.23, 0.20, 0.18, 0.20, 0.23]]
>>> vol_surface = cppfm.VolatilitySurface(strikes, maturities, ivs, curve)

>>> # Create Heston model
>>> heston = cppfm.HestonModel(
...     spot=100, discount_curve=curve,
...     v0=0.04, kappa=2.0, vbar=0.04, sigma_v=0.3, rho=-0.7
... )

>>> # Run SLV simulation
>>> time_steps = list(np.linspace(0, 1, 252))  # daily steps
>>> sim = cppfm.HestonSLVSimulator(
...     heston, vol_surface, time_steps, num_paths=10000
... )
>>> paths = sim.simulate_full()  # Full path history for hedging

>>> # Analyze terminal distribution
>>> terminal_spots = np.array([p[-1][0] for p in paths])
>>> print(f"Mean S_T: {terminal_spots.mean():.2f}")
"""

# Try to import the compiled module
try:
    from cppfm import *
except ImportError:
    # If not yet built, provide helpful message
    import sys
    import os

    _msg = """
    CppFM Python bindings not yet built.

    To build:
        cd /path/to/CppFM
        mkdir build_python && cd build_python
        cmake -DBUILD_PYTHON_BINDINGS=ON ..
        make -j

    Then add to Python path:
        import sys
        sys.path.insert(0, '/path/to/CppFM/build_python/python')
        import cppfm
    """

    print(_msg, file=sys.stderr)
    raise
