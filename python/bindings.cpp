/**
 * Python bindings for CppFM library using pybind11
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // std::vector, std::pair automatic conversion
#include <pybind11/operators.h> // For operator overloading

// CppFM headers
#include "../DiscountCurve.h"
#include "../VolatilitySurface.h"
#include "../Model.h"
#include "../PathSimulator2D.h"
#include "../BlackScholesFormulas.h"
#include "../FinancialInstrument.h"
// Utils.h causes crashes - skip for now

namespace py = pybind11;

// =============================================================================
// Trampoline class for DiscountCurve (abstract base class)
// =============================================================================
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

// =============================================================================
// Main module definition
// =============================================================================
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

    // =========================================================================
    // DiscountCurve classes
    // =========================================================================

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

    // =========================================================================
    // VolatilitySurface enums
    // =========================================================================

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

    // =========================================================================
    // VolatilitySurface
    // =========================================================================

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

    // =========================================================================
    // VolatilitySurfaceBuilder
    // =========================================================================

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

    // =========================================================================
    // HestonModel
    // =========================================================================

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

    // =========================================================================
    // HestonSLVPathSimulator2D
    // =========================================================================
    //
    // remember that C++ class stores timeSteps by reference, so we need a
    // wrapper class to keep the vector alive for the lifetime of the simulator.
    //

    // Wrapper class that owns both the time_steps vector and the simulator
    struct HestonSLVSimulatorWrapper
    {
        std::vector<double> time_steps;
        std::unique_ptr<HestonSLVPathSimulator2D> simulator;
        std::unique_ptr<VolatilitySurface> owned_vol_surface; // owned copy for external surface mode

        // Constructor 1: HestonModel only (uses internal HestonLocalVol)
        HestonSLVSimulatorWrapper(const HestonModel &model,
                                  std::vector<double> ts,
                                  size_t num_paths,
                                  size_t num_bins,
                                  size_t seed)
            : time_steps(std::move(ts)), owned_vol_surface(nullptr)
        {
            simulator = std::make_unique<HestonSLVPathSimulator2D>(
                model, time_steps, num_paths, num_bins, seed);
        }

        // Constructor 2: HestonModel + external VolatilitySurface
        HestonSLVSimulatorWrapper(const HestonModel &model,
                                  const VolatilitySurface &volSurface,
                                  std::vector<double> ts,
                                  size_t num_paths,
                                  size_t num_bins,
                                  size_t seed)
            : time_steps(std::move(ts))
        {
            // Clone the surface so wrapper owns it (Python object may go out of scope)
            owned_vol_surface = volSurface.clone();
            simulator = std::make_unique<HestonSLVPathSimulator2D>(
                model, owned_vol_surface.get(), time_steps, num_paths, num_bins, seed);
        }

        std::vector<std::pair<double, double>> simulate(bool parallel = true)
        {
            if (parallel)
            {
                return simulator->simulateAllPathsParallel();
            }
            return simulator->simulateAllPaths();
        }

        std::vector<std::vector<std::pair<double, double>>> simulate_full(bool parallel = true)
        {
            if (parallel)
            {
                return simulator->simulateAllPathsFullParallel();
            }
            return simulator->simulateAllPathsFull();
        }

        size_t calibrate_leverage(size_t max_iters = 10, double tol = 1e-3, double damping = 0.5)
        {
            return simulator->calibrateLeverage(max_iters, tol, damping);
        }

        bool is_calibrated() const
        {
            return simulator->isCalibrated();
        }

        std::vector<double> get_calibration_errors() const
        {
            return simulator->calibrationErrors();
        }

        void reset_calibration()
        {
            simulator->resetCalibration();
        }
    };

    // =========================================================================
    // Pure Heston QE Simulator - TEMPORARILY DISABLED
    // =========================================================================
    /*
    struct HestonQESimulatorWrapper { ... };
    py::class_<HestonQESimulatorWrapper>(m, "HestonQESimulator", ...);
    */

    // =========================================================================
    // HestonSLVSimulator
    // =========================================================================

    py::class_<HestonSLVSimulatorWrapper>(m, "HestonSLVSimulator",
                                          R"pbdoc(
            Heston Stochastic Local Volatility (SLV) Monte Carlo Simulator.

            Implements the van der Stoep et al. (2013) methodology for
            calibrating SLV models. Uses:
            - QE (Quadratic-Exponential) scheme for variance
            - Binning algorithm for conditional expectation E[V|S]
            - Leverage function L^2(t,S) = sigma_LV^2 / E[V|S]

            Two modes of operation:
            1. Without vol_surface: Uses internal HestonLocalVol (COS-based Dupire)
               - For validating SLV implementation (L≈1 everywhere)
            2. With vol_surface: Uses external VolatilitySurface for local vol
               - For fitting arbitrary market smiles

            Args:
                model: HestonModel instance
                vol_surface: Optional VolatilitySurface for external local vol
                time_steps: List of time points [0, t1, t2, ..., T]
                num_paths: Number of Monte Carlo paths
                num_bins: Number of bins for E[V|S] (default: 20)
                seed: Random seed (default: 42)

            Example (internal Dupire):
                sim = cppfm.HestonSLVSimulator(
                    heston, time_steps=list(np.linspace(0, 1, 252)),
                    num_paths=50000
                )

            Example (external surface):
                sim = cppfm.HestonSLVSimulator(
                    heston, vol_surface, time_steps=list(np.linspace(0, 1, 252)),
                    num_paths=50000
                )
        )pbdoc")
        // Constructor 1: HestonModel only (internal Dupire)
        .def(py::init<const HestonModel &,
                      std::vector<double>,
                      size_t,
                      size_t,
                      size_t>(),
             py::arg("model"),
             py::arg("time_steps"),
             py::arg("num_paths"),
             py::arg("num_bins") = 20,
             py::arg("seed") = 42)
        // Constructor 2: HestonModel + external VolatilitySurface
        .def(py::init<const HestonModel &,
                      const VolatilitySurface &,
                      std::vector<double>,
                      size_t,
                      size_t,
                      size_t>(),
             py::arg("model"),
             py::arg("vol_surface"),
             py::arg("time_steps"),
             py::arg("num_paths"),
             py::arg("num_bins") = 20,
             py::arg("seed") = 42)
        .def("simulate", &HestonSLVSimulatorWrapper::simulate,
             py::arg("parallel") = true,
             "Simulate and return terminal values (S_T, V_T)")
        .def("simulate_full", &HestonSLVSimulatorWrapper::simulate_full,
             py::arg("parallel") = true,
             "Simulate and return full path history")
        .def("calibrate_leverage", &HestonSLVSimulatorWrapper::calibrate_leverage,
             py::arg("max_iters") = 10,
             py::arg("tol") = 1e-3,
             py::arg("damping") = 0.5,
             "Run iterative leverage calibration (Van der Stoep)")
        .def("is_calibrated", &HestonSLVSimulatorWrapper::is_calibrated,
             "Check if leverage has been calibrated")
        .def("get_calibration_errors", &HestonSLVSimulatorWrapper::get_calibration_errors,
             "Get convergence history")
        .def("reset_calibration", &HestonSLVSimulatorWrapper::reset_calibration,
             "Reset calibration to single-pass mode");

    // =========================================================================
    // Option type enum (must be defined before BlackScholes that uses it)
    // =========================================================================

    py::enum_<Option::Type>(m, "OptionType", "Option type enumeration")
        .value("Call", Option::Type::Call)
        .value("Put", Option::Type::Put)
        .export_values();

    // =========================================================================
    // BlackScholesFormulas (static utility functions)
    // =========================================================================

    // Use explicit function pointer casts for overloaded functions
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

    // Heston pricing via COS method (HestonCF and COSPricer from Utils.h via PathSimulator2D.h)
    m.def("heston_call_price", [](double S0, double K, double r, double T,
                                   double kappa, double vbar, double sigma,
                                   double rho, double v0) {
        HestonCF cf(kappa, vbar, sigma, rho, v0, r, T);
        return COSPricer::callPrice(S0, K, r, T, cf);
    }, py::arg("S0"), py::arg("K"), py::arg("r"), py::arg("T"),
       py::arg("kappa"), py::arg("vbar"), py::arg("sigma"),
       py::arg("rho"), py::arg("v0"),
       "Heston call price using COS method");

    // =========================================================================
    // Version info
    // =========================================================================

    m.attr("__version__") = "0.1.0";
}
