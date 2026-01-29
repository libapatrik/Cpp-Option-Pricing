// pybind11 module entry point
#include "bindings_common.h"

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        CppFM Python Bindings
        ---------------------

        Python interface to the CppFM quantitative finance library.
        Provides access to:
        - Discount curves
        - Volatility surfaces (with Dupire local vol)
        - Heston model
        - Heston SLV simulator
        - Black-Scholes formulas

        Example:
            import cppfm

            # Create discount curve
            curve = cppfm.FlatDiscountCurve(0.05)

            # Create volatility surface
            vol_surface = cppfm.VolatilitySurface(
                strikes, maturities, volatilities, curve
            )

            # Create Heston model
            heston = cppfm.HestonModel(
                spot=100, discount_curve=curve, v0=0.04,
                kappa=2.0, vbar=0.04, sigma_v=0.3, rho=-0.7
            )

            # Run SLV simulation
            sim = cppfm.HestonSLVSimulator(
                heston, time_steps, num_paths=10000
            )
            paths = sim.simulate_full()  # Full path history
    )pbdoc";

    bind_market(m);
    bind_models(m);
    bind_simulators(m);
    bind_pricers(m);

    m.attr("__version__") = "0.1.0";
}
