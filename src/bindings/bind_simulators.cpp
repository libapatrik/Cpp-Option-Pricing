// bindings for HestonSLVSimulator
#include "bindings_common.h"
#include <cppfm/models/Model.h>
#include <cppfm/market/VolatilitySurface.h>
#include <cppfm/simulators/PathSimulator2D.h>

// wrapper owns both time_steps vector and simulator (C++ stores by reference)
struct HestonSLVSimulatorWrapper
{
    std::vector<double> time_steps;
    std::unique_ptr<HestonSLVPathSimulator2D> simulator;
    std::unique_ptr<VolatilitySurface> owned_vol_surface; // for external surface mode

    // Simulate with control variate (Heston as control)
    py::dict simulate_with_cv(double strike, double analytical_price)
    {
        auto result = simulator->simulateWithControlVariate(strike, analytical_price);
        py::dict d;
        d["mean"] = result.mean;
        d["std_error"] = result.stdError;
        d["variance"] = result.variance;
        d["ci95_lower"] = result.ci95Lower;
        d["ci95_upper"] = result.ci95Upper;
        d["num_paths"] = result.numPaths;
        return d;
    }

    // constructor 1: HestonModel only (uses internal HestonLocalVol)
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

    // constructor 2: HestonModel + external VolatilitySurface
    HestonSLVSimulatorWrapper(const HestonModel &model,
                              const VolatilitySurface &volSurface,
                              std::vector<double> ts,
                              size_t num_paths,
                              size_t num_bins,
                              size_t seed)
        : time_steps(std::move(ts))
    {
        // clone so wrapper owns it (Python object may go out of scope)
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

    // Mixing factor: η=0 pure Heston, η=1 full SLV
    void set_mixing_factor(double eta)
    {
        simulator->setMixingFactor(eta);
    }

    double get_mixing_factor() const
    {
        return simulator->getMixingFactor();
    }
};

void bind_simulators(py::module_ &m)
{
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
        // constructor 1: HestonModel only (internal Dupire)
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
        // constructor 2: HestonModel + external VolatilitySurface
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
             "Reset calibration to single-pass mode")
        .def_property("mixing_factor",
             &HestonSLVSimulatorWrapper::get_mixing_factor,
             &HestonSLVSimulatorWrapper::set_mixing_factor,
             "Mixing factor η ∈ [0,1]: η=0 pure Heston, η=1 full SLV")
        .def("simulate_with_cv", &HestonSLVSimulatorWrapper::simulate_with_cv,
             py::arg("strike"), py::arg("analytical_price"),
             R"pbdoc(
            Simulate with control variate (Heston as control).

            Evolves SLV and pure Heston paths simultaneously using same Brownians.
            High correlation → effective variance reduction (typically 2-5x).

            Args:
                strike: Option strike for payoff computation
                analytical_price: Analytical Heston price (e.g., from COS method)

            Returns:
                dict with keys: mean, std_error, variance, ci95_lower, ci95_upper, num_paths
        )pbdoc");
}
