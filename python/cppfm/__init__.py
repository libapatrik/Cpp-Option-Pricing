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
    HestonParams,
    HestonSliceData,
    HestonCalibrationResult,
    LMOptions,
    calibrate_heston,
    heston_grid_search,
    SviParams,
    SmileFitResult,
    calibrate_svi,
    SsviParams,
    SsviCalibrationResult,
    SsviCalibrationOptions,
    calibrate_ssvi,
    SsviSurface,
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
    "HestonParams",
    "HestonSliceData",
    "HestonCalibrationResult",
    "LMOptions",
    "calibrate_heston",
    "heston_grid_search",
    "SviParams",
    "SmileFitResult",
    "calibrate_svi",
    "SsviParams",
    "SsviCalibrationResult",
    "SsviCalibrationOptions",
    "calibrate_ssvi",
    "SsviSurface",
    "__version__",
]
