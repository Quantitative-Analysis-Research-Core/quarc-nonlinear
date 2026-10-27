// False nearest neighbours: the k-d tree and search of matlab/fnn.m. See quarc.hpp.
//
// A faithful port, not just an exact nearest-neighbour search: the tree is
// split, and the candidate list merged, exactly as fnn.m does, so that when two
// neighbours are equally near the same one wins (MATLAB's sort is stable).
//
// Copyright (c) 2021-2026 Quantitative Analysis Research Core,
// Center for Human Movement Variability, University of Nebraska at Omaha.
// MIT licence. See LICENSE.txt.

#include "quarc/quarc.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace quarc {

namespace {

constexpr std::size_t kBin = 512;   // leaf size, as kd_part(y, z, 512)
constexpr double kInf = std::numeric_limits<double>::infinity();

struct Node {
    std::size_t lo, hi;      // row range in the model, inclusive
    std::size_t left = 0, right = 0;  // children (0 = terminal)
    std::size_t disc = 0;    // splitting dimension
    double cut = 0;
};

// Rows are stored point-major: row r is y[r*d .. r*d+d-1].
struct Tree {
    std::size_t d;
    std::vector<double> y, z;
    std::vector<Node> nodes;
};

Tree kd_part(std::vector<double> y, std::vector<double> z, std::size_t n, std::size_t d) {
    Tree t{d, std::move(y), std::move(z), {}};
    t.nodes.push_back({0, n - 1});
    std::vector<double> range(d), key, ybuf, zbuf;
    std::vector<std::size_t> order, idx(d);
    for (std::size_t node = 0; node < t.nodes.size(); ++node) {
        const std::size_t lo = t.nodes[node].lo, hi = t.nodes[node].hi;
        const std::size_t len = hi - lo + 1;
        for (std::size_t c = 0; c < d; ++c) {
            double mx = -kInf, mn = kInf;
            for (std::size_t r = lo; r <= hi; ++r) {
                const double v = t.y[r * d + c];
                if (v > mx) mx = v;
                if (v < mn) mn = v;
            }
            range[c] = mx - mn;
        }
        if (!(*std::max_element(range.begin(), range.end()) > 0) || len < kBin) continue;

        // [~, index] = sort(range); split on index(d): the widest dimension,
        // the last of equals under a stable ascending sort.
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(),
                         [&](std::size_t a, std::size_t b) { return range[a] < range[b]; });
        const std::size_t disc = idx[d - 1];

        key.resize(len);
        for (std::size_t r = 0; r < len; ++r) key[r] = t.y[(lo + r) * d + disc];
        order.resize(len);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) { return key[a] < key[b]; });

        double cut = (len % 2) ? key[order[(len - 1) / 2]]
                               : (key[order[len / 2 - 1]] + key[order[len / 2]]) / 2;
        std::size_t nleft = 0;
        for (std::size_t r = 0; r < len; ++r) nleft += key[order[r]] <= cut;
        if (nleft == len) {  // right child would be empty: use y < cut instead
            double mx = -kInf;
            nleft = 0;
            for (std::size_t r = 0; r < len; ++r) {
                const double v = key[order[r]];
                if (v < cut) { ++nleft; if (v > mx) mx = v; }
            }
            cut = (cut + mx) / 2;
        }

        ybuf.resize(len * d);
        zbuf.resize(len);
        for (std::size_t r = 0; r < len; ++r) {
            std::copy_n(&t.y[(lo + order[r]) * d], d, &ybuf[r * d]);
            zbuf[r] = t.z[lo + order[r]];
        }
        std::copy(ybuf.begin(), ybuf.end(), t.y.begin() + lo * d);
        std::copy(zbuf.begin(), zbuf.end(), t.z.begin() + lo);

        const std::size_t L = t.nodes.size();
        t.nodes[node].disc = disc;
        t.nodes[node].cut = cut;
        t.nodes[node].left = L;
        t.nodes[node].right = L + 1;
        t.nodes.push_back({lo, lo + nleft - 1});
        t.nodes.push_back({lo + nleft, hi});
    }
    return t;
}

// The two nearest candidates (pqd, pqz in fnn.m); `have` counts real entries.
struct Best {
    double d[2] = {kInf, kInf};
    double z[2] = {0, 0};
    int have = 0;
};

struct Search {
    const Tree& t;
    const double* yq;
    std::vector<double> lower, upper;
    Best best;
    bool done = false;

    // Stable top-2 insert: a later candidate only displaces a strictly farther one.
    static void insert(Best& b, double d, double z) {
        if (b.have == 0 || d < b.d[0]) {
            b.d[1] = b.d[0]; b.z[1] = b.z[0];
            b.d[0] = d; b.z[0] = z;
        } else if (b.have == 1 || d < b.d[1]) {
            b.d[1] = d; b.z[1] = z;
        } else {
            return;
        }
        if (b.have < 2) ++b.have;
    }

    void leaf(const Node& nd) {
        // [dist; pqd] sorted stably: this leaf's points first, then the old best.
        Best nb;
        const std::size_t d = t.d;
        for (std::size_t r = nd.lo; r <= nd.hi; ++r) {
            double s = 0;
            for (std::size_t c = 0; c < d; ++c) {
                const double e = t.y[r * d + c] - yq[c];
                s += e * e;
            }
            insert(nb, std::sqrt(s), t.z[r]);
        }
        for (int k = 0; k < best.have; ++k) insert(nb, best.d[k], best.z[k]);
        // Placeholders (Inf, no z) only fill what the real candidates leave empty.
        if (nb.have < 2) nb.d[1] = kInf;
        best = nb;

        const double r2 = best.d[1];
        bool inside = true;
        for (std::size_t c = 0; c < d; ++c) {
            if (!(std::fabs(yq[c] - lower[c]) > r2 && std::fabs(yq[c] - upper[c]) > r2)) {
                inside = false;
                break;
            }
        }
        if (inside) done = true;
    }

    bool overlap() const {
        const double dist = best.d[1] * best.d[1];
        double s = 0;
        for (std::size_t c = 0; c < t.d; ++c) {
            if (yq[c] < lower[c]) {
                s += (yq[c] - lower[c]) * (yq[c] - lower[c]);
                if (s > dist) return false;
            } else if (yq[c] > upper[c]) {
                s += (yq[c] - upper[c]) * (yq[c] - upper[c]);
                if (s > dist) return false;
            }
        }
        return true;
    }

    void run(std::size_t node) {
        if (done) return;
        const Node& nd = t.nodes[node];
        if (nd.left == 0) { leaf(nd); return; }
        const std::size_t c = nd.disc;
        const bool goLeft = yq[c] <= nd.cut;
        double& near = goLeft ? upper[c] : lower[c];
        double tmp = near;
        near = nd.cut;
        run(goLeft ? nd.left : nd.right);
        near = tmp;
        if (done) return;
        double& far = goLeft ? lower[c] : upper[c];
        tmp = far;
        far = nd.cut;
        if (overlap()) run(goLeft ? nd.right : nd.left);
        far = tmp;
    }
};

}  // namespace

double fnn_false_fraction(const double* x, std::size_t len, std::size_t delay, std::size_t j,
                          std::size_t n, double rtol, double atol, double ra) {
    if (j < 1 || delay < 1 || n < 2) throw std::invalid_argument("fnn: bad dimension, delay or length.");
    if (n + delay * j > len) throw std::invalid_argument("fnn: series too short.");
    // y = [x(1:n), x(1+delay:n+delay), ...] (j columns); z = x(1+delay*j : n+delay*j).
    std::vector<double> y(n * j), z(n);
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < j; ++c) y[r * j + c] = x[r + delay * c];
        z[r] = x[r + delay * j];
    }
    const std::vector<double> yorig = y, zorig = z;
    const Tree t = kd_part(std::move(y), std::move(z), n, j);

    std::size_t falses = 0;
    for (std::size_t i = 0; i < n; ++i) {
        Search s{t, &yorig[i * j], std::vector<double>(j, -kInf), std::vector<double>(j, kInf)};
        s.run(0);
        const double distance = s.best.z[0] - s.best.z[1];
        bool L = std::fabs(distance) > s.best.d[1] * rtol;
        if (std::sqrt(s.best.d[1] * s.best.d[1] + distance * distance) / ra > atol) L = true;
        falses += L;
    }
    return static_cast<double>(falses) / static_cast<double>(n);
}

}  // namespace quarc
