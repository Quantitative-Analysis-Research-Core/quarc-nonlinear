// Python bindings for the C++ core. See cpp/include/quarc/quarc.hpp.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include "quarc/quarc.hpp"

namespace nb = nanobind;
using Vec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

NB_MODULE(_core, m) {
    m.doc() = "Compiled kernels shared with the MATLAB and R packages.";
    // std::invalid_argument from the core arrives in Python as ValueError.
    m.def("ent_samp", [](Vec x, int dim, double radius) {
        return quarc::ent_samp(x.data(), x.shape(0), dim, radius);
    }, nb::arg("x"), nb::arg("dim"), nb::arg("radius"),
       "Sample entropy with an absolute radius.");
}
