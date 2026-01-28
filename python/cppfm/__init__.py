"""
CppFM: Quantitative Finance Library

Provides:
- Discount curves (FlatDiscountCurve)
- Volatility surfaces with Dupire local vol
- Heston stochastic volatility model
- Heston SLV simulator
- Black-Scholes formulas and Greeks
"""

from cppfm._core import (
    DiscountCurve,
    FlatDiscountCurve,
    VolatilitySurface,
    VolatilitySurfaceBuilder,
    SmileInterpolation,
    MaturityInterpolation,
    HestonModel,
    HestonSLVSimulator,
    # HestonQESimulator,  # temporarily disabled
    OptionType,
    bs_call_price,
    bs_put_price,
    bs_delta,
    bs_gamma,
    bs_vega,
    bs_theta,
    bs_implied_volatility,
    heston_call_price,
)

from cppfm._core import __version__

__all__ = [
    "DiscountCurve",
    "FlatDiscountCurve",
    "VolatilitySurface",
    "VolatilitySurfaceBuilder",
    "SmileInterpolation",
    "MaturityInterpolation",
    "HestonModel",
    "HestonSLVSimulator",
    # "HestonQESimulator",  # temporarily disabled
    "OptionType",
    "bs_call_price",
    "bs_put_price",
    "bs_delta",
    "bs_gamma",
    "bs_vega",
    "bs_theta",
    "bs_implied_volatility",
    "heston_call_price",
    "__version__",
]
