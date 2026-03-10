// bindings for calibration: Heston + SVI/SSVI
#include "bindings_common.h"
#include <cppfm/calibration/HestonCalibration.h>
#include <cppfm/calibration/Optimizer.h>
#include <cppfm/market/SviCalibration.h>

void bind_calibration(py::module_ &m)
{
    // HestonParams
    py::class_<HestonParams>(m, "HestonParams",
        "Heston model parameters: v0, kappa, vbar, sigma_v, rho")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double>(),
             py::arg("v0") = 0.04, py::arg("kappa") = 1.5,
             py::arg("vbar") = 0.04, py::arg("sigma_v") = 0.3,
             py::arg("rho") = -0.7)
        .def_readwrite("v0", &HestonParams::v0)
        .def_readwrite("kappa", &HestonParams::kappa)
        .def_readwrite("vbar", &HestonParams::vbar)
        .def_readwrite("sigma_v", &HestonParams::sigma_v)
        .def_readwrite("rho", &HestonParams::rho)
        .def("to_vector", &HestonParams::toVector)
        .def("from_vector", &HestonParams::fromVector)
        .def("satisfies_feller", &HestonParams::satisfiesFellerCondition)
        .def_static("lower_bounds", &HestonParams::lowerBounds)
        .def_static("upper_bounds", &HestonParams::upperBounds)
        .def("__repr__", [](const HestonParams &p) {
            return "HestonParams(v0=" + std::to_string(p.v0)
                + ", kappa=" + std::to_string(p.kappa)
                + ", vbar=" + std::to_string(p.vbar)
                + ", sigma_v=" + std::to_string(p.sigma_v)
                + ", rho=" + std::to_string(p.rho) + ")";
        });

    // HestonSliceData
    py::class_<HestonSliceData>(m, "HestonSliceData",
        "Per-maturity market data: T, strikes, marketIVs")
        .def(py::init<>())
        .def_readwrite("T", &HestonSliceData::T)
        .def_readwrite("strikes", &HestonSliceData::strikes)
        .def_readwrite("market_ivs", &HestonSliceData::marketIVs);

    // HestonCalibrationResult
    py::class_<HestonCalibrationResult>(m, "HestonCalibrationResult")
        .def_readonly("params", &HestonCalibrationResult::params)
        .def_readonly("rmse", &HestonCalibrationResult::rmse)
        .def_readonly("slice_rmse", &HestonCalibrationResult::sliceRmse)
        .def_readonly("iterations", &HestonCalibrationResult::iterations)
        .def_readonly("converged", &HestonCalibrationResult::converged)
        .def_readonly("message", &HestonCalibrationResult::message);

    // LMOptions (already exists in C++, expose for calibration tuning)
    py::class_<LMOptions>(m, "LMOptions",
        "Levenberg-Marquardt optimizer options")
        .def(py::init<>())
        .def_readwrite("tol", &LMOptions::tol)
        .def_readwrite("grad_tol", &LMOptions::gradTol)
        .def_readwrite("max_iter", &LMOptions::maxIter)
        .def_readwrite("lambda0", &LMOptions::lambda0)
        .def_readwrite("verbose", &LMOptions::verbose);

    // calibrateHeston
    m.def("calibrate_heston", &calibrateHeston,
          py::arg("slices"), py::arg("S0"), py::arg("r"),
          py::arg("guess") = HestonParams{},
          py::arg("opts") = LMOptions{},
          "Calibrate Heston model to market implied vols across multiple maturities");

    // hestonGridSearch
    m.def("heston_grid_search", &hestonGridSearch,
          py::arg("slices"), py::arg("S0"), py::arg("r"),
          py::arg("grid_points") = 3,
          "Coarse grid search for Heston starting point");

    m.def("heston_grid_search_parallel", &hestonGridSearchParallel,
          py::arg("slices"), py::arg("S0"), py::arg("r"),
          py::arg("grid_points") = 3,
          "Parallel grid search for Heston starting point");

    // ========================================================================
    // SVI / SSVI
    // ========================================================================

    // SviParams
    py::class_<SviParams>(m, "SviParams",
        "SVI smile parameters: a, b, rho, m, sigma")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double>(),
             py::arg("a") = 0.04, py::arg("b") = 0.01,
             py::arg("rho") = -0.3, py::arg("m") = 0.0,
             py::arg("sigma") = 0.2)
        .def_readwrite("a", &SviParams::a)
        .def_readwrite("b", &SviParams::b)
        .def_readwrite("rho", &SviParams::rho)
        .def_readwrite("m", &SviParams::m)
        .def_readwrite("sigma", &SviParams::sigma)
        .def("total_variance", &SviParams::totalVariance, py::arg("k"))
        .def("implied_vol", &SviParams::impliedVol, py::arg("k"), py::arg("T"))
        .def("dw", &SviParams::dw, py::arg("k"))
        .def("d2w", &SviParams::d2w, py::arg("k"))
        .def("g_function", &SviParams::gFunction, py::arg("k"))
        .def("satisfies_constraints", &SviParams::satisfiesConstraints)
        .def("to_vector", &SviParams::toVector)
        .def("from_vector", &SviParams::fromVector);

    // SmileFitResult
    py::class_<SmileFitResult>(m, "SmileFitResult")
        .def_readonly("rmse", &SmileFitResult::rmse)
        .def_readonly("converged", &SmileFitResult::converged)
        .def_property_readonly("params", [](const SmileFitResult &r) -> const SmileParams& {
            return *r.params;
        }, py::return_value_policy::reference_internal);

    // calibrateSmile
    m.def("calibrate_svi", [](const std::vector<double> &strikes,
                              const std::vector<double> &vols,
                              double forward, double T,
                              const LMOptions &opts) {
        SviParams guess;
        return calibrateSmile(guess, strikes, vols, forward, T, opts);
    },  py::arg("strikes"), py::arg("vols"),
        py::arg("forward"), py::arg("T"),
        py::arg("opts") = LMOptions{},
        "Calibrate SVI smile to market vols for a single maturity");

    // SsviParams
    py::class_<SsviParams>(m, "SsviParams",
        "SSVI surface parameters: rho, eta, gamma (+ per-slice theta)")
        .def(py::init<>())
        .def(py::init<double, double, double>(),
             py::arg("rho") = -0.5, py::arg("eta") = 0.1,
             py::arg("gamma") = 0.5)
        .def_readwrite("rho", &SsviParams::rho)
        .def_readwrite("eta", &SsviParams::eta)
        .def_readwrite("gamma", &SsviParams::gamma)
        .def_readwrite("theta", &SsviParams::theta)
        .def("phi", py::overload_cast<>(&SsviParams::phi, py::const_),
             "phi(theta) using current theta")
        .def("phi", py::overload_cast<double>(&SsviParams::phi, py::const_),
             py::arg("theta"), "phi for arbitrary theta")
        .def("total_variance", &SsviParams::totalVariance, py::arg("k"))
        .def("implied_vol", &SsviParams::impliedVol, py::arg("k"), py::arg("T"))
        .def("dw", &SsviParams::dw, py::arg("k"))
        .def("d2w", &SsviParams::d2w, py::arg("k"))
        .def("g_function", &SsviParams::gFunction, py::arg("k"))
        .def("satisfies_constraints", &SsviParams::satisfiesConstraints)
        .def("satisfies_butterfly_constraint", &SsviParams::satisfiesButterflyConstraint,
             py::arg("theta"))
        .def("to_vector", &SsviParams::toVector)
        .def("from_vector", &SsviParams::fromVector);

    // SsviCalibrationResult
    py::class_<SsviCalibrationResult>(m, "SsviCalibrationResult")
        .def_readonly("params", &SsviCalibrationResult::params)
        .def_readonly("thetas", &SsviCalibrationResult::thetas)
        .def_readonly("maturities", &SsviCalibrationResult::maturities)
        .def_readonly("rmse", &SsviCalibrationResult::rmse)
        .def_readonly("converged", &SsviCalibrationResult::converged)
        .def_readonly("arbitrage_free", &SsviCalibrationResult::arbitrageFree)
        .def("build_surface", &SsviCalibrationResult::buildSurface,
             py::arg("strike_grid"), py::arg("forwards"),
             py::arg("discount_curve"))
        .def("build_analytical_surface", &SsviCalibrationResult::buildAnalyticalSurface,
             py::arg("forwards"),
             "Build SsviSurface for analytical Dupire local vol (no grid)");

    // SsviCalibrationOptions
    py::class_<SsviCalibrationOptions>(m, "SsviCalibrationOptions")
        .def(py::init<>())
        .def_readwrite("lm_opts", &SsviCalibrationOptions::lmOpts)
        .def_readwrite("enforce_monotonicity", &SsviCalibrationOptions::enforceMonotonicity)
        .def_readwrite("use_grid_search", &SsviCalibrationOptions::useGridSearch);

    // calibrateSsvi
    m.def("calibrate_ssvi", &calibrateSsvi,
          py::arg("strikes_per_maturity"), py::arg("vols_per_maturity"),
          py::arg("forwards"), py::arg("maturities"),
          py::arg("opts") = SsviCalibrationOptions{},
          "Calibrate SSVI surface: per-slice SVI -> extract thetas -> global fit");

    // SsviSurface -- analytical Dupire local vol from SSVI
    py::class_<SsviSurface>(m, "SsviSurface",
        "Analytical Dupire local vol from SSVI parametric form")
        .def(py::init<const SsviCalibrationResult&, const std::vector<double>&>(),
             py::arg("result"), py::arg("forwards"))
        .def("local_volatility", &SsviSurface::localVolatility,
             py::arg("spot"), py::arg("time"))
        .def("implied_volatility", &SsviSurface::impliedVolatility,
             py::arg("strike"), py::arg("time"))
        .def("local_variance", &SsviSurface::localVariance,
             py::arg("k"), py::arg("time"));

}
