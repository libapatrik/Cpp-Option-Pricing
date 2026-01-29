// common includes for all binding modules
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

// forward declarations for bind functions
void bind_market(py::module_ &m);
void bind_models(py::module_ &m);
void bind_simulators(py::module_ &m);
void bind_pricers(py::module_ &m);
