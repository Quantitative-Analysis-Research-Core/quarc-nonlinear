"""QUARC Nonlinear: nonlinear time-series analysis, ported from MATLAB.

Function names, arguments and defaults follow the MATLAB functions of the
same name, which are the specification. Indices in results are 0-based.
"""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

from .ami import ami
from .ami_copula import ami_copula
from .ami_histogram import ami_histogram
from .ami_kde import ami_kde
from .ami_knn import ami_knn
from .ami_partition import ami_partition
from .chaos_library import chaos_library
from .corr_dim import corr_dim
from .crqa import crqa
from .dfa import dfa
from .embed import embed
from .ent import ent
from .ent_ap import ent_ap
from .ent_ms_plus import ent_ms_plus
from .ent_permu import ent_permu
from .ent_samp import ent_samp
from .ent_symbolic import ent_symbolic
from .ent_weighted import ent_weighted
from .ent_xap import ent_xap
from .ent_xsamp import ent_xsamp
from .fgn_sim import fgn_sim
from .fnn import fnn
from .hilbert import hilbert
from .jrqa import jrqa
from .line_hist import line_hist
from .lyapunov import lyapunov
from .lye_r import lye_r
from .lye_w import lye_w
from .mdrqa import mdrqa
from .ode23 import ode23
from .ode45 import ode45
from .psr import psr
from .rel_phase_cont import rel_phase_cont
from .rel_phase_disc import rel_phase_disc
from .rqa import rqa
from .set_radius import set_radius
from .surr_find_rho import surr_find_rho
from .surr_pseudo_periodic import surr_pseudo_periodic
from .surr_theiler import surr_theiler

__all__ = [
    "ami",
    "ami_copula",
    "ami_histogram",
    "ami_kde",
    "ami_knn",
    "ami_partition",
    "chaos_library",
    "corr_dim",
    "crqa",
    "dfa",
    "embed",
    "ent",
    "ent_ap",
    "ent_ms_plus",
    "ent_permu",
    "ent_samp",
    "ent_symbolic",
    "ent_weighted",
    "ent_xap",
    "ent_xsamp",
    "fgn_sim",
    "fnn",
    "hilbert",
    "jrqa",
    "line_hist",
    "lyapunov",
    "lye_r",
    "lye_w",
    "mdrqa",
    "ode23",
    "ode45",
    "psr",
    "rel_phase_cont",
    "rel_phase_disc",
    "rqa",
    "set_radius",
    "surr_find_rho",
    "surr_pseudo_periodic",
    "surr_theiler",
]
