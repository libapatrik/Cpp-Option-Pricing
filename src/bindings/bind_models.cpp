// bindings for HestonModel
#include "bindings_common.h"
#include <cppfm/market/DiscountCurve.h>
#include <cppfm/models/Model.h>

void bind_models(py::module_ &m)
{
    py::class_<HestonModel>(m, "HestonModel",
                            R"pbdoc(
            Heston Stochastic Volatility Model.

            2D SDE system:
                dS = r*S*dt + sqrt(V)*S*dW^S
                dV = kappa*(vbar - V)*dt + sigma_v*sqrt(V)*dW^V
                d<W^S, W^V> = rho*dt

            Args:
                spot: Initial asset price S_0
                discount_curve: Risk-free rate curve
                v0: Initial variance V_0 > 0
                kappa: Mean reversion speed kappa > 0
                vbar: Long-term variance theta > 0
                sigma_v: Vol-of-vol sigma_v > 0
                rho: Correlation rho in [-1, 1]
        )pbdoc")
        .def(py::init<double, const DiscountCurve &, double, double, double, double, double>(),
             py::arg("spot"),
             py::arg("discount_curve"),
             py::arg("v0"),
             py::arg("kappa"),
             py::arg("vbar"),
             py::arg("sigma_v"),
             py::arg("rho"))
        .def_property_readonly("spot", [](const HestonModel &m)
                               { return m.initValue(); }, "Initial spot price S_0")
        .def_property_readonly("v0", &HestonModel::v0, "Initial variance V_0")
        .def_property_readonly("kappa", &HestonModel::kappa, "Mean reversion speed kappa")
        .def_property_readonly("vbar", &HestonModel::vbar, "Long-term variance theta (vbar)")
        .def_property_readonly("sigma_v", &HestonModel::sigma_v, "Volatility of volatility sigma_v")
        .def_property_readonly("rho", &HestonModel::rho, "Correlation rho")
        .def_property_readonly("risk_free_rate", &HestonModel::riskFreeRate, "Risk-free rate r")
        .def("correlation", &HestonModel::correlation, "Get correlation coefficient rho")
        .def("satisfies_feller", &HestonModel::satisfiesFellerCondition, "Check Feller condition: 2*kappa*vbar >= sigma_v^2")
        .def("__repr__", [](const HestonModel &m)
             { return "<HestonModel S0=" + std::to_string(m.initValue()) +
                      " v0=" + std::to_string(m.v0()) +
                      " kappa=" + std::to_string(m.kappa()) +
                      " vbar=" + std::to_string(m.vbar()) +
                      " sigma_v=" + std::to_string(m.sigma_v()) +
                      " rho=" + std::to_string(m.rho()) + ">"; });
}
