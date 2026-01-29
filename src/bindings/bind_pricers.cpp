// bindings for BlackScholesFormulas and heston_call_price
#include "bindings_common.h"
#include <cppfm/pricers/BlackScholesFormulas.h>
#include <cppfm/market/FinancialInstrument.h>
#include <cppfm/utils/Utils.h>

void bind_pricers(py::module_ &m)
{
    // Option type enum
    py::enum_<Option::Type>(m, "OptionType", "Option type enumeration")
        .value("Call", Option::Type::Call)
        .value("Put", Option::Type::Put)
        .export_values();

    // explicit function pointer casts for overloaded functions
    m.def("bs_call_price",
          static_cast<double (*)(double, double, double, double, double)>(
              &BlackScholesFormulas::callPrice),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"),
          "Black-Scholes call option price");

    m.def("bs_put_price",
          static_cast<double (*)(double, double, double, double, double)>(
              &BlackScholesFormulas::putPrice),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"),
          "Black-Scholes put option price");

    m.def("bs_delta",
          static_cast<double (*)(double, double, double, double, double, Option::Type)>(
              &BlackScholesFormulas::delta),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"), py::arg("option_type"),
          "Option delta");

    m.def("bs_gamma",
          static_cast<double (*)(double, double, double, double, double)>(
              &BlackScholesFormulas::gamma),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"),
          "Option gamma");

    m.def("bs_vega",
          static_cast<double (*)(double, double, double, double, double)>(
              &BlackScholesFormulas::vega),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"),
          "Option vega");

    m.def("bs_theta",
          static_cast<double (*)(double, double, double, double, double, Option::Type)>(
              &BlackScholesFormulas::theta),
          py::arg("spot"), py::arg("strike"), py::arg("rate"),
          py::arg("volatility"), py::arg("maturity"), py::arg("option_type"),
          "Option theta");

    m.def("bs_implied_volatility",
          static_cast<double (*)(double, double, double, double, double, Option::Type, double, size_t, double)>(
              &BlackScholesFormulas::impliedVolatility),
          py::arg("spot"), py::arg("strike"), py::arg("rate"), py::arg("maturity"),
          py::arg("market_price"), py::arg("option_type"),
          py::arg("initial_guess") = 0.2, py::arg("max_iter") = 100, py::arg("tol") = 1e-8,
          "Compute implied volatility from option price");

    // Heston pricing via COS method
    m.def("heston_call_price", [](double S0, double K, double r, double T,
                                   double kappa, double vbar, double sigma,
                                   double rho, double v0) {
        HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
        return COSPricer::callPrice(S0, K, r, T, cf);
    }, py::arg("S0"), py::arg("K"), py::arg("r"), py::arg("T"),
       py::arg("kappa"), py::arg("vbar"), py::arg("sigma"),
       py::arg("rho"), py::arg("v0"),
       "Heston call price using COS method");
}
