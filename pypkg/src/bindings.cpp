// Python bindings for the C++ core. See cpp/include/quarc/quarc.hpp.
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include "quarc/quarc.hpp"

#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;
using Vec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Bytes = nb::ndarray<const std::uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

NB_MODULE(_core, m) {
    m.doc() = "Compiled kernels shared with the MATLAB and R packages.";
    // std::invalid_argument from the core arrives in Python as ValueError.
    m.def("ent_samp", [](Vec x, int dim, double radius) {
        return quarc::ent_samp(x.data(), x.shape(0), dim, radius);
    }, nb::arg("x"), nb::arg("dim"), nb::arg("radius"),
       "Sample entropy with an absolute radius.");

    // Matrices cross as flat column-major (Fortran-order) buffers.
    m.def("line_hist", [](Bytes a, std::size_t n, std::size_t theiler, bool cross) {
        if (a.shape(0) != n * n) throw std::invalid_argument("line_hist: a must hold n*n values.");
        auto h = quarc::line_hist(a.data(), n, theiler, cross);
        return std::make_pair(std::move(h.diag), std::move(h.vert));
    }, nb::arg("a"), nb::arg("n"), nb::arg("theiler"), nb::arg("cross"),
       "Diagonal and vertical line lengths of an n-by-n column-major 0/1 matrix.");

    m.def("corr_dim_range", [](Vec x, int dim, std::size_t delay) {
        double lo = 0, hi = 0;
        quarc::corr_dim_range(x.data(), x.shape(0), dim, delay, lo, hi);
        return std::make_pair(lo, hi);
    }, nb::arg("x"), nb::arg("dim"), nb::arg("delay"),
       "Smallest and largest distance over corr_dim's pairs.");

    m.def("corr_dim_counts", [](Vec x, int dim, std::size_t delay, Vec edges) {
        return quarc::corr_dim_counts(x.data(), x.shape(0), dim, delay, edges.data(), edges.shape(0));
    }, nb::arg("x"), nb::arg("dim"), nb::arg("delay"), nb::arg("edges"),
       "histc counts of corr_dim's pair distances against edges.");

    m.def("fnn_false_fraction", [](Vec x, std::size_t delay, std::size_t j, std::size_t n,
                                   double rtol, double atol, double ra) {
        return quarc::fnn_false_fraction(x.data(), x.shape(0), delay, j, n, rtol, atol, ra);
    }, nb::arg("x"), nb::arg("delay"), nb::arg("j"), nb::arg("n"), nb::arg("rtol"),
       nb::arg("atol"), nb::arg("ra"),
       "Fraction of false nearest neighbours when dimension j + 1 is added.");

    m.def("lye_r_divergence", [](Vec Y, std::size_t M, std::size_t D, std::size_t theiler) {
        if (Y.shape(0) != M * D) throw std::invalid_argument("lye_r: Y must hold M*D values.");
        auto r = quarc::lye_r_divergence(Y.data(), M, D, theiler);
        return std::make_tuple(std::move(r.nn), std::move(r.avg), r.count);
    }, nb::arg("Y"), nb::arg("M"), nb::arg("D"), nb::arg("theiler"),
       "Rosenstein nearest neighbours (0-based) and mean log divergence of an M-by-D column-major Y.");
}
