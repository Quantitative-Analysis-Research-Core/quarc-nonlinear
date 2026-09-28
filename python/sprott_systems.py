"""The 62 systems of Sprott (2003) Appendix A, as vectorized numpy vector fields.

WHY A PORT. The reference corpus recomputes every exponent under one method, and
the dysts vector fields live in Python. Rather than move 126 dysts systems into
MATLAB -- which python/export_dysts.py rejected, because a mistyped coefficient
still produces a plausible chaotic series -- the 62 Sprott systems come here, so
a single engine covers both catalogues.

HOW THE TRANSCRIPTION IS CHECKED. The same objection applies in this direction,
and the answer is that these systems are not judged by eye. Every entry carries
Sprott's published lambda1, and the acceptance gate is reproducing it to the
precision he prints. A mistyped coefficient moves the exponent; four matching
decimals is a far stronger check on the transcription than reading the line
twice. Entries whose lambda is a closed form (logistic, cusp, LCG, Arnold's cat)
are checked against the closed form instead, and the dissipative systems are
checked additionally against the trace identity, sum(lambda) = mean divergence.

Ported from tests/matlab/+quarctest/sprott_catalog.m, which holds the provenance
notes for the transcription decisions -- the Burke-Shaw sign, the Chua bracket,
the Ricker initial condition -- and those notes are repeated here where they
affect the equations.

VERIFIED AGAINST APPENDIX A, 24 September 2026. Every equation, parameter and
initial condition in this file was read against the printed appendix, not
merely inferred from agreeing with the published exponent. That distinction
earned its keep: ACT (A.5.7) carried a real error inherited from the MATLAB
catalogue -- dz/dt takes -delta*alpha*z, and the alpha had been dropped, giving
-1.5z instead of -2.7z. It passed every internal check, because the wrong
system is a perfectly good dynamical system; only the published value caught
it, and only at 31% error.

Three deviations from the printed page are deliberate and stand:

  burke_shaw    A.5.14 prints dz/dt = -Uxy + V. That form diverges: every
                trajectory leaves the attractor at any step size tried. The
                +Uxy form used here gives lambda1 = 2.2442 and a third exponent
                of -13.245, against the appendix's own 2.2499 and -13.2499, so
                the printed sign is a typo in the book.
  gauss_map     A.1.7 gives X0 = 0.1, whose reciprocal is exactly 10 in binary,
                so mod(1/x, 1) is exactly 0 and the orbit dies on iteration 1.
                An irrational start is used instead.
  ricker        The same reciprocal-of-0.1 argument, per the MATLAB notes.

Everything else -- all 62 systems across A.1 to A.6 -- matches the appendix
exactly.

CONVENTIONS. Every f is vectorized over leading axes: y has shape (..., D) and
the return has the same shape, so an ensemble of trajectories advances in one
call. Flows take (t, y) and return dy/dt; maps take (t, y) and return y_next,
ignoring t. All lambdas are nats per unit time, per iteration for maps; Sprott
states that Appendix A is base-e throughout.

Copyright (c) 2021-2026 Quantitative Analysis Research Core,
Center for Human Movement Variability, University of Nebraska at Omaha.
MIT licence. See LICENSE.txt.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf, erfinv

__all__ = ["SYSTEMS", "by_name", "usable"]

TWO_PI = 2.0 * np.pi


def _s(*components):
    """Stack per-component expressions back into a (..., D) state."""
    return np.stack(components, axis=-1)


# ------------------------------------------------------------------ A.1 maps

def f_logistic(t, y):
    x = y[..., 0]
    return _s(4 * x * (1 - x))


def f_sine_map(t, y):
    return _s(np.sin(np.pi * y[..., 0]))


def f_tent(t, y):
    x = y[..., 0]
    return _s(2 * np.minimum(x, 1 - x))


def f_lcg(t, y):
    return _s(np.mod(7141 * y[..., 0] + 54773, 259200))


def f_cubic_map(t, y):
    x = y[..., 0]
    return _s(3 * x * (1 - x ** 2))


def f_ricker(t, y):
    x = y[..., 0]
    return _s(20 * x * np.exp(-x))


def f_gauss_map(t, y):
    return _s(np.mod(1.0 / y[..., 0], 1.0))


def f_cusp(t, y):
    return _s(1 - 2 * np.sqrt(np.abs(y[..., 0])))


def f_gauss_white(t, y):
    x = y[..., 0]
    return _s(np.sqrt(2) * erfinv(1 - 2 * erf(x / np.sqrt(2))))


def f_pinchers(t, y):
    return _s(np.abs(np.tanh(2 * (y[..., 0] - 0.5))))


def f_spence(t, y):
    return _s(np.abs(np.log(y[..., 0])))


def f_sine_circle(t, y):
    x = y[..., 0]
    return _s(np.mod(x + 0.5 - (2 / TWO_PI) * np.sin(TWO_PI * x), 1.0))


# ------------------------------------------------------- A.2 dissipative maps

def f_henon(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(1 - 1.4 * x1 ** 2 + 0.3 * x2, x1)


def f_lozi(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(1 - 1.7 * np.abs(x1) + 0.5 * x2, x1)


def f_delayed_logistic(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(2.27 * x1 * (1 - x2), x1)


def f_tinkerbell(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x1 ** 2 - x2 ** 2 + 0.9 * x1 - 0.6 * x2,
              2 * x1 * x2 + 2 * x1 + 0.5 * x2)


def f_burgers(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(0.75 * x1 - x2 ** 2, 1.75 * x2 + x1 * x2)


def f_holmes_cubic(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -0.2 * x1 + 2.77 * x2 - x2 ** 3)


def f_kaplan_yorke(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(np.mod(2 * x1, 1.0), np.mod(0.2 * x2 + np.cos(4 * np.pi * x1), 1.0))


def f_dissipative_standard(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    yn = np.mod(0.1 * x2 + 8.8 * np.sin(x1), TWO_PI)
    return _s(np.mod(x1 + yn, TWO_PI), yn)


def f_ikeda(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    phi = 0.4 - 6.0 / (1 + x1 ** 2 + x2 ** 2)
    return _s(1 + 0.9 * (x1 * np.cos(phi) - x2 * np.sin(phi)),
              0.9 * (x1 * np.sin(phi) + x2 * np.cos(phi)))


def f_sinai(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(np.mod(x1 + x2 + 0.1 * np.cos(TWO_PI * x2), 1.0),
              np.mod(x1 + 2 * x2, 1.0))


def f_predator_prey(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x1 * np.exp(3 * (1 - x1) - 5 * x2), x1 * (1 - np.exp(-5 * x2)))


# ------------------------------------------------------ A.3 conservative maps

def f_chirikov(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    yn = np.mod(x2 + 1.0 * np.sin(x1), TWO_PI)
    return _s(np.mod(x1 + yn, TWO_PI), yn)


def f_henon_area(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    ca = 0.24
    sa = np.sqrt(1 - ca ** 2)
    return _s(x1 * ca - (x2 - x1 ** 2) * sa, x1 * sa + (x2 - x1 ** 2) * ca)


def f_arnold_cat(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(np.mod(x1 + x2, 1.0), np.mod(x1 + 2 * x2, 1.0))


def f_gingerbreadman(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(1 + np.abs(x1) - x2, x1)


def f_chaotic_web(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    al = np.pi / 2
    w = x2 + 1.0 * np.sin(x1)
    return _s(x1 * np.cos(al) - w * np.sin(al), x1 * np.sin(al) + w * np.cos(al))


def f_lorenz3d_map(t, y):
    x1, x2, x3 = y[..., 0], y[..., 1], y[..., 2]
    return _s(x1 * x2 - x3, x1, x2)


# ---------------------------------------------------------- A.4 driven flows

def f_damped_pendulum(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -np.sin(x1) - 0.05 * x2 + 0.6 * np.sin(0.7 * t))


def f_driven_vdp(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -x1 + 3 * (1 - x1 ** 2) * x2 + 5 * np.sin(1.788 * t))


def f_shaw_vdp(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2 + np.sin(2 * t), -x1 + 1.0 * (1 - x1 ** 2) * x2)


def f_brusselator(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x1 ** 2 * x2 - (1.2 + 1) * x1 + 0.4 + 0.05 * np.sin(0.8 * t),
              -x1 ** 2 * x2 + 1.2 * x1)


def f_ueda(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -x1 ** 3 - 0.05 * x2 + 7.5 * np.sin(t))


def f_duffing_two_well(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -x1 ** 3 + x1 - 0.25 * x2 + 0.4 * np.sin(t))


def f_duffing_vdp(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, 0.2 * (1 - 8 * x1 ** 2) * x2 - x1 ** 3 + 0.35 * np.sin(1.02 * t))


def f_rayleigh_duffing(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, 0.2 * (1 - 4 * x2 ** 2) * x2 - x1 ** 3 + 0.3 * np.sin(1.1 * t))


# ------------------------------------------------------ A.5 autonomous flows

def f_lorenz(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(10 * (yy - x), -x * z + 28 * x - yy, x * yy - (8 / 3) * z)


def f_rossler(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-yy - z, x + 0.2 * yy, 0.2 + z * (x - 5.7))


def f_diffusionless_lorenz(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-yy - x, -x * z, x * yy + 1)


def f_complex_butterfly(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(0.55 * (yy - x), -z * np.sign(x), np.abs(x) - 1)


def f_chen(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(35 * (yy - x), (28 - 35) * x - x * z + 28 * yy, x * yy - 3 * z)


def f_hadley(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-yy ** 2 - z ** 2 - 0.25 * x + 0.25 * 8,
              x * yy - 4 * x * z - yy + 1,
              4 * x * yy + x * z - z)


def f_act(t, y):
    # dz/dt takes -delta*alpha*z, per A.5.7: with alpha = 1.8 and delta = 1.5
    # that is -2.7z. The MATLAB catalogue dropped the alpha and used -1.5z,
    # which this port inherited; the error showed up as a converged exponent
    # 31% below the published one, with every internal check passing, because
    # the wrong system is a perfectly good system.
    alpha, beta, delta, mu = 1.8, -0.07, 1.5, 0.02
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(alpha * (x - yy),
              -4 * alpha * yy + x * z + mu * x ** 3,
              -delta * alpha * z + x * yy + beta * z ** 2)


def f_rabinovich_fabrikant(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy * (z - 1 + x ** 2) + 0.87 * x,
              x * (3 * z + 1 - x ** 2) + 0.87 * yy,
              -2 * z * (1.1 + x * yy))


def f_rigid_body(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-0.4 * x + yy + 10 * yy * z,
              -x - 0.4 * yy + 5 * x * z,
              0.175 * z - 5 * x * yy)


def f_chua(t, y):
    # The bracketed terms are ADDED, per the appendix. Subtracting them
    # instead drives the system to a fixed point (sprott_catalog.m).
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    al, be, a, b = 9.0, 100 / 7, 8 / 7, 5 / 7
    h = b * x + 0.5 * (a - b) * (np.abs(x + 1) - np.abs(x - 1))
    return _s(al * (yy - x + h), x - yy + z, -be * yy)


def f_moore_spiegel(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -z - (6 - 20 + 20 * x ** 2) * yy - 6 * x)


def f_thomas(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-0.18 * x + np.sin(yy), -0.18 * yy + np.sin(z), -0.18 * z + np.sin(x))


def f_halvorsen(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-1.27 * x - 4 * yy - 4 * z - yy ** 2,
              -1.27 * yy - 4 * z - 4 * x - z ** 2,
              -1.27 * z - 4 * x - 4 * yy - x ** 2)


def f_burke_shaw(t, y):
    # dz/dt takes +U*x*y; the appendix line reads "-Uxy + V", which diverges
    # at every step size tried (sprott_catalog.m).
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-10 * (x + yy), -10 * x * z - yy, 10 * x * yy + 13)


def f_rucklidge(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(-2 * x + 6.7 * yy - yy * z, x, -z + yy ** 2)


def f_windmi(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -0.7 * z - yy + 2.5 - np.exp(x))


def f_simplest_quadratic(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -2.017 * z + yy ** 2 - x)


def f_simplest_cubic(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -2.028 * z + x * yy ** 2 - x)


def f_simplest_piecewise(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -0.6 * z - yy + np.abs(x) - 1)


def f_double_scroll(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, z, -0.8 * (z + yy + x - np.sign(x)))


# ---------------------------------------------------- A.6 conservative flows

def f_driven_pendulum(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -np.sin(x1) + 1.0 * np.sin(0.5 * t))


def f_simplest_driven(t, y):
    x1, x2 = y[..., 0], y[..., 1]
    return _s(x2, -x1 ** 3 + np.sin(1.88 * t))


def f_nose_hoover(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(yy, -x + yy * z, 1 - yy ** 2)


def f_labyrinth(t, y):
    x, yy, z = y[..., 0], y[..., 1], y[..., 2]
    return _s(np.sin(yy), np.sin(z), np.sin(x))


def f_henon_heiles(t, y):
    x1, x2, x3, x4 = y[..., 0], y[..., 1], y[..., 2], y[..., 3]
    return _s(x3, x4, -x1 - 2 * x1 * x2, -x2 - x1 ** 2 + x2 ** 2)


# ------------------------------------------------------------------ catalogue

def _e(name, pretty, section, category, kind, f, x0, lam1, tier,
       dt=None, usable=True, note="", wrap=None):
    # wrap gives the period of each state component for systems closed with
    # mod, so separations are measured on the torus rather than across the
    # seam. See lyap_spectrum._separation.
    return dict(name=name, pretty=pretty, section=section, category=category,
                kind=kind, f=f, x0=np.asarray(x0, dtype=np.float64), lam1=lam1,
                tier=tier, dt=dt, usable=usable, note=note,
                wrap=None if wrap is None else np.asarray(wrap, dtype=np.float64))


_LN2 = float(np.log(2))

SYSTEMS = [
    # A.1 noninvertible maps
    _e("logistic", "Logistic map", "A.1.1", "noninvertible_map", "map",
       f_logistic, [0.1], _LN2, "exact"),
    _e("sine_map", "Sine map", "A.1.2", "noninvertible_map", "map",
       f_sine_map, [0.1], 0.689067, "numerical"),
    _e("tent", "Tent map", "A.1.3", "noninvertible_map", "map",
       f_tent, [1 / np.sqrt(2)], _LN2, "exact", usable=False,
       note="Binary shift: exhausts the mantissa and reaches exactly 0 after "
            "~52 iterations. Not fixable in double precision."),
    _e("lcg", "Linear congruential generator", "A.1.4", "noninvertible_map", "map",
       f_lcg, [0.0], float(np.log(7141)), "exact", wrap=[259200.0]),
    _e("cubic_map", "Cubic map", "A.1.5", "noninvertible_map", "map",
       f_cubic_map, [0.1], 1.0986122883, "numerical"),
    # Sprott gives x0 = 0.1, but 1/0.1 is exactly 10 in binary, so mod(.,1) is
    # exactly 0 and the orbit dies on iteration 1 for the Gauss map below.
    _e("ricker", "Ricker's population model", "A.1.6", "noninvertible_map", "map",
       f_ricker, [0.1], 0.384846, "numerical"),
    _e("gauss_map", "Gauss map", "A.1.7", "noninvertible_map", "map",
       f_gauss_map, [0.3141592653589793], 2.373445, "numerical", wrap=[1.0]),
    _e("cusp", "Cusp map", "A.1.8", "noninvertible_map", "map",
       f_cusp, [0.5], 0.5, "exact"),
    _e("gauss_white", "Gaussian white chaotic map", "A.1.9", "noninvertible_map", "map",
       f_gauss_white, [1.0], _LN2, "exact", usable=False,
       note="Leaves the domain on step 2 as transcribed; erfinv returns NaN."),
    _e("pinchers", "Pinchers map", "A.1.10", "noninvertible_map", "map",
       f_pinchers, [0.0], 0.467944, "numerical"),
    _e("spence", "Spence map", "A.1.11", "noninvertible_map", "map",
       f_spence, [0.5], float("inf"), "exact", usable=False,
       note="lambda -> Inf. No finite reference exists."),
    _e("sine_circle", "Sine-circle map", "A.1.12", "noninvertible_map", "map",
       f_sine_circle, [0.1], 0.353863, "numerical", wrap=[1.0]),

    # A.2 dissipative maps
    _e("henon", "Henon map", "A.2.1", "dissipative_map", "map",
       f_henon, [0.0, 0.9], 0.41922, "numerical"),
    _e("lozi", "Lozi map", "A.2.2", "dissipative_map", "map",
       f_lozi, [-0.1, 0.1], 0.47023, "numerical"),
    _e("delayed_logistic", "Delayed logistic map", "A.2.3", "dissipative_map", "map",
       f_delayed_logistic, [0.001, 0.001], 0.18312, "numerical"),
    _e("tinkerbell", "Tinkerbell map", "A.2.4", "dissipative_map", "map",
       f_tinkerbell, [0.0, 0.5], 0.18997, "numerical"),
    _e("burgers", "Burgers' map", "A.2.5", "dissipative_map", "map",
       f_burgers, [-0.1, 0.1], 0.12076, "numerical"),
    _e("holmes_cubic", "Holmes cubic map", "A.2.6", "dissipative_map", "map",
       f_holmes_cubic, [1.6, 0.0], 0.59458, "numerical"),
    _e("kaplan_yorke", "Kaplan-Yorke map", "A.2.7", "dissipative_map", "map",
       f_kaplan_yorke, [1 / np.sqrt(2), -0.4], _LN2, "exact", wrap=[1.0, 1.0], usable=False,
       note="x -> 2x mod 1 is the same binary shift as the tent map."),
    _e("dissipative_standard", "Dissipative standard map", "A.2.8", "dissipative_map", "map",
       f_dissipative_standard, [0.1, 0.1], 1.46995, "numerical", wrap=[TWO_PI, TWO_PI]),
    _e("ikeda", "Ikeda map", "A.2.9", "dissipative_map", "map",
       f_ikeda, [0.0, 0.0], 0.50760, "numerical"),
    _e("sinai", "Sinai map", "A.2.10", "dissipative_map", "map",
       f_sinai, [0.5, 0.5], 0.95946, "numerical", wrap=[1.0, 1.0]),
    _e("predator_prey", "Discrete predator-prey map", "A.2.11", "dissipative_map", "map",
       f_predator_prey, [0.5, 0.5], 0.19664, "numerical"),

    # A.3 conservative maps
    _e("chirikov", "Chirikov (standard) map", "A.3.1", "conservative_map", "map",
       f_chirikov, [0.0, 6.0], 0.10497, "numerical", wrap=[TWO_PI, TWO_PI]),
    _e("henon_area", "Henon area-preserving quadratic map", "A.3.2", "conservative_map", "map",
       f_henon_area, [0.6, 0.13], 0.00643, "numerical"),
    _e("arnold_cat", "Arnold's cat map", "A.3.3", "conservative_map", "map",
       f_arnold_cat, [0.0, 1 / np.sqrt(2)], float(np.log((3 + np.sqrt(5)) / 2)), "exact", wrap=[1.0, 1.0]),
    _e("gingerbreadman", "Gingerbreadman map", "A.3.4", "conservative_map", "map",
       f_gingerbreadman, [0.5, 3.7], 0.07339, "numerical"),
    _e("chaotic_web", "Chaotic web map", "A.3.5", "conservative_map", "map",
       f_chaotic_web, [0.0, 3.0], 0.04847, "numerical"),
    _e("lorenz3d_map", "Lorenz three-dimensional chaotic map", "A.3.6", "conservative_map", "map",
       f_lorenz3d_map, [0.5, 0.5, -1.0], 0.07456, "numerical"),

    # A.4 driven flows
    _e("damped_pendulum", "Damped driven pendulum", "A.4.1", "driven_flow", "flow",
       f_damped_pendulum, [0.0, 0.0], 0.1414, "numerical", dt=0.02),
    _e("driven_vdp", "Driven van der Pol oscillator", "A.4.2", "driven_flow", "flow",
       f_driven_vdp, [-1.9, 0.0], 0.1933, "numerical", dt=0.005),
    _e("shaw_vdp", "Shaw-van der Pol oscillator", "A.4.3", "driven_flow", "flow",
       f_shaw_vdp, [1.3, 0.0], 0.1180, "numerical", dt=0.01),
    _e("brusselator", "Forced Brusselator", "A.4.4", "driven_flow", "flow",
       f_brusselator, [0.3, 2.0], 0.0140, "numerical", dt=0.02),
    _e("ueda", "Ueda oscillator", "A.4.5", "driven_flow", "flow",
       f_ueda, [2.5, 0.0], 0.1034, "numerical", dt=0.01),
    _e("duffing_two_well", "Duffing's two-well oscillator", "A.4.6", "driven_flow", "flow",
       f_duffing_two_well, [0.2, 0.0], 0.1572, "numerical", dt=0.02),
    _e("duffing_vdp", "Duffing-van der Pol oscillator", "A.4.7", "driven_flow", "flow",
       f_duffing_vdp, [0.2, -0.2], 0.0963, "numerical", dt=0.02),
    _e("rayleigh_duffing", "Rayleigh-Duffing oscillator", "A.4.8", "driven_flow", "flow",
       f_rayleigh_duffing, [0.3, 0.0], 0.0912, "numerical", dt=0.02),

    # A.5 autonomous flows
    _e("lorenz", "Lorenz attractor", "A.5.1", "autonomous_flow", "flow",
       f_lorenz, [0.0, -0.01, 9.0], 0.9056, "numerical", dt=0.003),
    _e("rossler", "Rossler attractor", "A.5.2", "autonomous_flow", "flow",
       f_rossler, [-9.0, 0.0, 0.0], 0.0714, "numerical", dt=0.01),
    _e("diffusionless_lorenz", "Diffusionless Lorenz attractor", "A.5.3", "autonomous_flow", "flow",
       f_diffusionless_lorenz, [1.0, -1.0, 0.01], 0.2101, "numerical", dt=0.01),
    _e("complex_butterfly", "Complex butterfly", "A.5.4", "autonomous_flow", "flow",
       f_complex_butterfly, [0.2, 0.0, 0.0], 0.1690, "numerical", dt=0.01),
    _e("chen", "Chen's system", "A.5.5", "autonomous_flow", "flow",
       f_chen, [-10.0, 0.0, 37.0], 2.0272, "numerical", dt=0.001),
    _e("hadley", "Hadley circulation", "A.5.6", "autonomous_flow", "flow",
       f_hadley, [0.0, 0.0, 1.3], 0.1665, "numerical", dt=0.01),
    _e("act", "ACT attractor", "A.5.7", "autonomous_flow", "flow",
       f_act, [0.5, 0.0, 0.0], 0.1634, "numerical", dt=0.01),
    _e("rabinovich_fabrikant", "Rabinovich-Fabrikant attractor", "A.5.8", "autonomous_flow", "flow",
       f_rabinovich_fabrikant, [-1.0, 0.0, 0.5], 0.1981, "numerical", dt=0.005),
    _e("rigid_body", "Linear feedback rigid body motion system", "A.5.9", "autonomous_flow", "flow",
       f_rigid_body, [0.6, 0.0, 0.0], 0.1421, "numerical", dt=0.01),
    _e("chua", "Chua's circuit", "A.5.10", "autonomous_flow", "flow",
       f_chua, [0.0, 0.0, 0.6], 0.3271, "numerical", dt=0.005),
    _e("moore_spiegel", "Moore-Spiegel oscillator", "A.5.11", "autonomous_flow", "flow",
       f_moore_spiegel, [0.1, 0.0, 0.0], 0.1119, "numerical", dt=0.005),
    _e("thomas", "Thomas' cyclically symmetric attractor", "A.5.12", "autonomous_flow", "flow",
       f_thomas, [0.1, 0.0, 0.0], 0.0349, "numerical", dt=0.05),
    _e("halvorsen", "Halvorsen's cyclically symmetric attractor", "A.5.13", "autonomous_flow", "flow",
       f_halvorsen, [-5.0, 0.0, 0.0], 0.7899, "numerical", dt=0.005),
    _e("burke_shaw", "Burke-Shaw attractor", "A.5.14", "autonomous_flow", "flow",
       f_burke_shaw, [0.6, 0.0, 0.0], 2.2499, "numerical", dt=0.001),
    _e("rucklidge", "Rucklidge attractor", "A.5.15", "autonomous_flow", "flow",
       f_rucklidge, [1.0, 0.0, 4.5], 0.0643, "numerical", dt=0.01),
    _e("windmi", "WINDMI attractor", "A.5.16", "autonomous_flow", "flow",
       f_windmi, [0.0, 0.8, 0.0], 0.0755, "numerical", dt=0.01),
    _e("simplest_quadratic", "Simplest quadratic chaotic flow", "A.5.17", "autonomous_flow", "flow",
       f_simplest_quadratic, [-0.9, 0.0, 0.5], 0.0551, "numerical", dt=0.01),
    _e("simplest_cubic", "Simplest cubic chaotic flow", "A.5.18", "autonomous_flow", "flow",
       f_simplest_cubic, [0.0, 0.96, 0.0], 0.0837, "numerical", dt=0.01),
    _e("simplest_piecewise", "Simplest piecewise linear chaotic flow", "A.5.19", "autonomous_flow", "flow",
       f_simplest_piecewise, [0.0, -0.7, 0.0], 0.0362, "numerical", dt=0.02),
    _e("double_scroll", "Double scroll", "A.5.20", "autonomous_flow", "flow",
       f_double_scroll, [0.01, 0.01, 0.0], 0.0497, "numerical", dt=0.02),

    # A.6 conservative flows
    _e("driven_pendulum", "Driven pendulum", "A.6.1", "conservative_flow", "flow",
       f_driven_pendulum, [0.0, 0.0], 0.1633, "numerical", dt=0.02),
    _e("simplest_driven", "Simplest driven chaotic flow", "A.6.2", "conservative_flow", "flow",
       f_simplest_driven, [0.0, 0.0], 0.0971, "numerical", dt=0.02),
    _e("nose_hoover", "Nose-Hoover oscillator", "A.6.3", "conservative_flow", "flow",
       f_nose_hoover, [0.0, 5.0, 0.0], 0.0138, "numerical", dt=0.02),
    _e("labyrinth", "Labyrinth chaos", "A.6.4", "conservative_flow", "flow",
       f_labyrinth, [0.1, 0.0, 0.0], 0.1402, "numerical", dt=0.05),
    _e("henon_heiles", "Henon-Heiles system", "A.6.5", "conservative_flow", "flow",
       f_henon_heiles, [0.499, 0.0, 0.0, 0.03160676], 0.0450, "numerical", dt=0.02),
]


def ic_dependent(entry):
    """Whether this system's exponent is a property of the starting point.

    A dissipative system draws every trajectory onto the same attractor, so its
    exponent belongs to the system and one number describes it. A conservative
    one does not: its state space is a patchwork of chaotic regions and regular
    islands that trajectories never cross between, so the exponent belongs to
    the initial condition. Henon-Heiles gave 0.00415 and 0.04487 here under
    nothing more than a different integration step, because the step nudged the
    trajectory into a different region. Both are correct measurements of
    different things.

    Such an entry still belongs in the catalogue -- excluding it would hide a
    real class of system -- but it carries this flag so a table can mark it
    rather than present a single number as if it were a property of the system.

    Uniformly hyperbolic exceptions are not flagged. Arnold's cat map is
    conservative but has no islands at all: every orbit has the same exponent,
    it is a closed form, and this port reproduces it to 3e-9.
    """
    return entry["category"].startswith("conservative") and entry["tier"] != "exact"


def by_name(name):
    for s in SYSTEMS:
        if s["name"] == name:
            return s
    raise KeyError(f"no Sprott system named {name!r}")


def usable():
    """The systems that can serve as benchmarks, in catalogue order."""
    return [s for s in SYSTEMS if s["usable"]]
