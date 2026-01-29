// bindings for DiscountCurve, VolatilitySurface, VolatilitySurfaceBuilder
#include "bindings_common.h"
#include <cppfm/market/DiscountCurve.h>
#include <cppfm/market/VolatilitySurface.h>

// trampoline for abstract DiscountCurve
class PyDiscountCurve : public DiscountCurve
{
public:
    using DiscountCurve::DiscountCurve;

    double discount(double time) const override
    {
        PYBIND11_OVERRIDE_PURE(double, DiscountCurve, discount, time);
    }

    DiscountCurve *clone() const override
    {
        PYBIND11_OVERRIDE_PURE(DiscountCurve *, DiscountCurve, clone);
    }

    double rate() const override
    {
        PYBIND11_OVERRIDE_PURE(double, DiscountCurve, rate);
    }
};

void bind_market(py::module_ &m)
{
    // DiscountCurve
    py::class_<DiscountCurve, PyDiscountCurve>(m, "DiscountCurve",
                                               "Abstract base class for discount curves")
        .def("discount", &DiscountCurve::discount, py::arg("time"),
             "Get discount factor B(0,t) for given time t")
        .def("rate", &DiscountCurve::rate,
             "Get the rate parameter")
        .def("instantaneous_rate", &DiscountCurve::instantaneousRate,
             py::arg("time"), py::arg("eps") = 1e-6,
             "Get instantaneous forward rate r(t) = -d/dt[log(B(t))]");

    py::class_<FlatDiscountCurve, DiscountCurve>(m, "FlatDiscountCurve",
                                                 "Flat (constant) interest rate discount curve: B(t) = exp(-r*t)")
        .def(py::init<double>(), py::arg("rate"),
             "Create flat discount curve with constant rate r")
        .def("discount", &FlatDiscountCurve::discount, py::arg("time"))
        .def("rate", &FlatDiscountCurve::rate)
        .def("__repr__", [](const FlatDiscountCurve &c)
             { return "<FlatDiscountCurve rate=" + std::to_string(c.rate()) + ">"; });

    // VolatilitySurface enums
    py::enum_<VolatilitySurface::SmileInterpolationType>(m, "SmileInterpolation",
                                                         "Interpolation method across strikes (smile dimension)")
        .value("Linear", VolatilitySurface::SmileInterpolationType::Linear,
               "Linear interpolation in volatility space")
        .value("CubicSpline", VolatilitySurface::SmileInterpolationType::CubicSpline,
               "Cubic spline interpolation (smoother, recommended)")
        .export_values();

    py::enum_<VolatilitySurface::MaturityInterpolationType>(m, "MaturityInterpolation",
                                                            "Interpolation method across maturities (term structure)")
        .value("Bilinear", VolatilitySurface::MaturityInterpolationType::Bilinear,
               "Simple linear interpolation in volatility space")
        .value("ForwardMoneyness", VolatilitySurface::MaturityInterpolationType::ForwardMoneyness,
               "Variance interpolation with constant forward moneyness (recommended)")
        .export_values();

    // VolatilitySurface
    py::class_<VolatilitySurface>(m, "VolatilitySurface",
                                  R"pbdoc(
            Implied volatility surface with Dupire local volatility support.

            Implements equation (1.4) from the documentation for IV surface
            and provides local volatility via Dupire formula.

            Args:
                strikes: List of strike prices
                maturities: List of maturities in years
                volatilities: 2D list of IVs [maturity_idx][strike_idx]
                discount_curve: DiscountCurve for forward calculations
                smile_interp: SmileInterpolation type (default: CubicSpline)
                maturity_interp: MaturityInterpolation type (default: ForwardMoneyness)
        )pbdoc")
        .def(py::init<const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<std::vector<double>> &,
                      const DiscountCurve &,
                      VolatilitySurface::SmileInterpolationType,
                      VolatilitySurface::MaturityInterpolationType>(),
             py::arg("strikes"),
             py::arg("maturities"),
             py::arg("volatilities"),
             py::arg("discount_curve"),
             py::arg("smile_interp") = VolatilitySurface::SmileInterpolationType::CubicSpline,
             py::arg("maturity_interp") = VolatilitySurface::MaturityInterpolationType::ForwardMoneyness)
        .def("implied_volatility", &VolatilitySurface::impliedVolatility,
             py::arg("strike"), py::arg("maturity"),
             "Get implied volatility sigma*(K, T)")
        .def("local_volatility", &VolatilitySurface::localVolatility,
             py::arg("spot"), py::arg("time"),
             "Get Dupire local volatility sigma_LV(S, t)")
        .def("get_bounds", &VolatilitySurface::getBounds,
             "Get surface bounds: ((min_K, max_K), (min_T, max_T))")
        .def_property_readonly("strikes", &VolatilitySurface::strikes,
                               "Get strike grid")
        .def_property_readonly("maturities", &VolatilitySurface::maturities,
                               "Get maturity grid")
        .def_property_readonly("volatilities", &VolatilitySurface::volatilities,
                               "Get volatility matrix");

    // VolatilitySurfaceBuilder
    py::class_<VolatilitySurfaceBuilder>(m, "VolatilitySurfaceBuilder",
                                         R"pbdoc(
            Builder for constructing VolatilitySurface incrementally.

            Example:
                builder = cppfm.VolatilitySurfaceBuilder()
                builder.add_strike(90).add_strike(100).add_strike(110)
                builder.add_maturity(0.25).add_maturity(0.5)
                builder.set_volatility(90, 0.25, 0.22)
                builder.set_volatility(100, 0.25, 0.20)
                # ... etc
                builder.set_discount_curve(curve)
                surface = builder.build()
        )pbdoc")
        .def(py::init<>())
        .def("add_strike", &VolatilitySurfaceBuilder::addStrike,
             py::arg("strike"), py::return_value_policy::reference,
             "Add a strike to the grid")
        .def("add_maturity", &VolatilitySurfaceBuilder::addMaturity,
             py::arg("maturity"), py::return_value_policy::reference,
             "Add a maturity to the grid")
        .def("set_volatility", &VolatilitySurfaceBuilder::setVolatility,
             py::arg("strike"), py::arg("maturity"), py::arg("volatility"),
             py::return_value_policy::reference,
             "Set implied volatility at (strike, maturity)")
        .def("set_smile_interpolation", &VolatilitySurfaceBuilder::setSmileInterpolationType,
             py::arg("interp_type"), py::return_value_policy::reference,
             "Set smile interpolation method")
        .def("set_maturity_interpolation", &VolatilitySurfaceBuilder::setMaturityInterpolationType,
             py::arg("interp_type"), py::return_value_policy::reference,
             "Set maturity interpolation method")
        .def("set_discount_curve", &VolatilitySurfaceBuilder::setDiscountCurve,
             py::arg("discount_curve"), py::return_value_policy::reference,
             "Set the discount curve")
        .def("build", [](VolatilitySurfaceBuilder &builder)
             { return builder.build(); }, py::return_value_policy::take_ownership, "Build and return the VolatilitySurface");
}
