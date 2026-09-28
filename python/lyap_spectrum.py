"""Lyapunov spectra by Sprott's method: nearby trajectories, Gram-Schmidt, renormalization.

WHY THIS EXISTS. The characterization suite compares estimator output against two
reference catalogues whose exponents were produced by different methods: Sprott
(2003) Appendix A used trajectory separation with periodic renormalization, and
the dysts database of Gilpin (2021) integrated the tangent space with QR
orthonormalization. On the 13 systems the catalogues share, the published
exponents disagree by a median of 22%, and by a factor of 21 for Thomas. Some of
that is method and some is control parameters, and as long as both sources are
mixed into one table no ratio can separate estimator error from reference
disagreement. This module recomputes every exponent under one method, from
recorded initial conditions and parameters, so the references stop being a
confound.

THE METHOD, following Sprott (2003) chapter 5. A fiducial trajectory is
integrated alongside k perturbed copies, each displaced from it by d0. After
every renormalization interval the k separation vectors are Gram-Schmidt
orthonormalized; the log of each vector's length before normalization
accumulates into the corresponding exponent, and each copy is placed back at
distance d0 from the fiducial along the orthonormalized direction. k = 1
recovers the two-trajectory largest-exponent method; k = D gives the spectrum.

Nothing here needs a Jacobian. That is the point: the same code runs over
Sprott's 62 hand-written vector fields and over the dysts catalogue, which
supplies no analytic derivatives for many of its systems.

VECTORIZATION. The state carries a leading ensemble axis, so E independent
initial conditions advance in one set of numpy operations. Wall-clock is set by
the number of sequential steps, not by E, so ensemble averaging is nearly free
and is how the ergodic average is driven down to Sprott's reported precision.

Copyright (c) 2021-2026 Quantitative Analysis Research Core,
Center for Human Movement Variability, University of Nebraska at Omaha.
MIT licence. See LICENSE.txt.
"""

from __future__ import annotations

import numpy as np

__all__ = ["spectrum", "SpectrumResult"]


class SpectrumResult:
    """Exponents plus the evidence needed to judge whether they converged.

    Attributes
    ----------
    exponents : (k,) float
        Ensemble-mean exponents, largest first, in nats per unit time (per
        iteration for maps).
    per_member : (E, k) float
        Each ensemble member's own time average. Their spread is the honest
        uncertainty: members are independent trajectories on the same
        attractor, so the standard error over members is a direct estimate of
        how far the reported mean may sit from the true exponent.
    sem : (k,) float
        Standard error of the mean across ensemble members.
    history : (n_checkpoints, k) float or None
        Running exponents at logarithmic checkpoints, for convergence plots.
    checkpoint_time : (n_checkpoints,) float or None
        Elapsed time (or iterations) at each checkpoint.
    """

    __slots__ = ("exponents", "per_member", "sem", "history", "checkpoint_time",
                 "n_steps", "elapsed", "d0", "renorm_every", "divergence",
                 "n_escaped", "n_saturated", "n_collapsed", "n_died")

    def __init__(self, exponents, per_member, sem, history, checkpoint_time,
                 n_steps, elapsed, d0, renorm_every, divergence=np.nan,
                 n_escaped=0, n_saturated=0, n_collapsed=0, n_died=0):
        self.exponents = exponents
        self.per_member = per_member
        self.sem = sem
        self.history = history
        self.checkpoint_time = checkpoint_time
        self.n_steps = n_steps
        self.elapsed = elapsed
        self.d0 = d0
        self.renorm_every = renorm_every
        # Time-averaged div f (flows) or log|det J| (maps) along the fiducial.
        # sum(exponents) must equal this; the gap is trace_error.
        self.divergence = divergence
        # Ensemble members whose state went non-finite and were dropped.
        self.n_escaped = n_escaped
        # Members dropped because their separation left the linear regime,
        # almost always a fiducial trajectory absorbed by a fixed point.
        self.n_saturated = n_saturated
        # Renormalization intervals where a separation underflowed to exactly
        # zero and was re-seeded. A handful over a long run is expected on maps
        # with a critical point; a large count means d0 is too small.
        self.n_collapsed = n_collapsed
        # Members whose fiducial trajectory reached a fixed point of the
        # floating-point map, or left the attractor, and were frozen there.
        self.n_died = n_died

    @property
    def trace_error(self):
        """sum(lambda) minus the time-averaged divergence.

        Only meaningful when the full spectrum was computed (k == D): a partial
        spectrum has no reason to sum to the divergence.
        """
        return float(np.sum(self.exponents) - self.divergence)

    def __repr__(self):
        lam = ", ".join(f"{x:.6f}" for x in self.exponents)
        sem = ", ".join(f"{x:.2g}" for x in self.sem)
        return f"SpectrumResult(exponents=[{lam}], sem=[{sem}], elapsed={self.elapsed:g})"


def _rk4(f, t, y, dt):
    """One fixed-step RK4 step. y has shape (..., D); f returns the same shape.

    Fixed step, not adaptive: Sprott's Appendix A values were computed with a
    fixed-step integrator, and an adaptive one would make the renormalization
    interval depend on the local dynamics, which is precisely what must be held
    constant across systems.
    """
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _gram_schmidt(v):
    """Orthogonalize separation vectors and return them with their lengths.

    v has shape (E, k, D): E ensemble members, k separation vectors of
    dimension D. Returns (unit vectors, lengths), both per member.

    Modified Gram-Schmidt, not classical: the vectors collapse towards the
    leading direction between renormalizations, so they arrive nearly parallel
    and the classical form loses orthogonality exactly when it matters most.
    """
    E, k, D = v.shape
    q = np.empty_like(v)
    norms = np.empty((E, k), dtype=v.dtype)
    for i in range(k):
        w = v[:, i, :]
        for j in range(i):
            w = w - np.sum(w * q[:, j, :], axis=-1, keepdims=True) * q[:, j, :]
        n = np.linalg.norm(w, axis=-1)
        norms[:, i] = n
        # A SEPARATION THAT UNDERFLOWS TO EXACTLY ZERO. Near a critical point
        # of a map, where f' vanishes, the two neighbouring states map to the
        # same double and the separation becomes exactly 0. The direction is
        # then undefined, and it must be re-seeded rather than left at zero:
        # placing the perturbed copy exactly on the fiducial makes the collapse
        # absorbing, the separation stays 0 for the rest of the run, and the
        # member accumulates a spurious log(1/d0) -- about +20 nats -- every
        # step thereafter. For the logistic map at r = 4 this needs an approach
        # within about 1e-8 of x = 0.5: invisible over 200k iterations and
        # near-certain over 2M, which is exactly the regime where longer runs
        # were supposed to buy precision. The event is reported by the caller
        # as n_collapsed, and the interval contributes nothing to the average.
        dead = n <= 0
        if dead.any():
            seed = np.zeros_like(w)
            seed[:, i % w.shape[-1]] = 1.0
            w = np.where(dead[:, None], seed, w)
            n_safe = np.where(dead, 1.0, n)
        else:
            n_safe = n
        q[:, i, :] = w / n_safe[:, None]
    return q, norms


def _separation(pert, fid, wrap):
    """Displacement of each perturbed copy from the fiducial, on the torus.

    Several of Sprott's maps close their state with mod: the standard map and
    the chaotic web on 2*pi, the cat and Sinai maps on 1. When the fiducial and
    a perturbed copy fall on opposite sides of such a seam they are neighbours
    on the torus but a full period apart in coordinates, and a naive difference
    reads that as an enormous divergence -- which would inflate every exponent
    on exactly the systems whose exponents are smallest. Wrapping the
    difference into the half-period measures the distance the dynamics
    actually sees.
    """
    v = pert - fid[:, None, :]
    if wrap is None:
        return v
    w = np.asarray(wrap, dtype=np.float64)
    live = np.isfinite(w) & (w > 0)
    if not live.any():
        return v
    ws = np.where(live, w, 1.0)
    return np.where(live, np.mod(v + 0.5 * ws, ws) - 0.5 * ws, v)


def _divergence(f, t, y, h_rel=1e-6):
    """Trace of the Jacobian at y, by central differences. y is (..., D).

    The independent check on a computed spectrum. For a flow, the sum of the
    Lyapunov exponents equals the time average of div f along the trajectory;
    for a map it equals the time average of log|det J|. That identity does not
    pass through the Gram-Schmidt at all -- it is read straight off the vector
    field -- so agreement between sum(lambda) and this average tests the
    renormalization machinery rather than restating it. Lorenz has constant
    divergence -(sigma + 1 + beta), which is why it lands at 1e-11 there and
    merely close on systems whose divergence varies over the attractor.

    Costs 2D evaluations of f, so it is sampled at renormalization times rather
    than every step.
    """
    D = y.shape[-1]
    out = np.zeros(y.shape[:-1], dtype=np.float64)
    for i in range(D):
        h = h_rel * np.maximum(1.0, np.abs(y[..., i]))
        e = np.zeros_like(y)
        e[..., i] = h
        out += (f(t, y + e)[..., i] - f(t, y - e)[..., i]) / (2 * h)
    return out


def _log_abs_det_jac(f, t, y, h_rel=1e-6):
    """log|det J| at y for a map, by central differences. y is (..., D)."""
    D = y.shape[-1]
    cols = []
    for i in range(D):
        h = h_rel * np.maximum(1.0, np.abs(y[..., i]))
        e = np.zeros_like(y)
        e[..., i] = h
        cols.append((f(t, y + e) - f(t, y - e)) / (2 * h[..., None]))
    J = np.stack(cols, axis=-1)                     # (..., D, D)
    sign, logdet = np.linalg.slogdet(J)
    return np.where(sign == 0, -np.inf, logdet)


def spectrum(f, x0, *, kind, k=None, dt=None, n_steps=100_000,
             transient_steps=10_000, renorm_every=1, d0=1e-8, ensemble=1,
             spread=1e-3, seed=0, checkpoints=0, t0=0.0, wrap=None,
             divergence=False):
    """Lyapunov spectrum of one system by nearby-trajectory Gram-Schmidt.

    Parameters
    ----------
    f : callable
        Flow: ``f(t, y) -> dy/dt``. Map: ``f(t, y) -> y_next``. In both cases y
        has shape (..., D) and f must be vectorized over the leading axes.
    x0 : (D,) array
        Initial condition, as published for the system.
    kind : {"flow", "map"}
        Flows are advanced by fixed-step RK4; maps are iterated directly.
    k : int, optional
        How many exponents. Defaults to D, the full spectrum.
    dt : float
        Integration step. Required for flows, ignored for maps.
    n_steps : int
        Steps over which exponents are accumulated, after the transient.
    transient_steps : int
        Steps discarded first, so the trajectory sits on the attractor and the
        separation vectors have aligned with the local Lyapunov directions.
        Both matter: starting accumulation with arbitrary separation directions
        biases every exponent.
    renorm_every : int
        Steps between Gram-Schmidt renormalizations. For maps this is 1. For
        flows it must be short enough that the separation stays in the linear
        regime and long enough that the vectors separate above roundoff.
    d0 : float
        Separation held between the fiducial and each perturbed trajectory.
    ensemble : int
        Independent initial conditions advanced together. Their mean is the
        reported exponent and their spread is its uncertainty.
    spread : float
        Ensemble members 2..E start from x0 perturbed by this fraction, then
        relax through the transient onto the attractor. Member 1 is x0 verbatim.
    seed : int
        Seed for those perturbations, so the corpus is reproducible.
    checkpoints : int
        If > 0, record running exponents at this many logarithmically spaced
        points, to show convergence rather than assert it.
    t0 : float
        Initial time, for nonautonomous systems.
    wrap : (D,) array, optional
        Period of each state component for systems closed with mod, or 0/inf
        for components that are not periodic. Separations are measured on that
        torus. See _separation.

    Returns
    -------
    SpectrumResult
    """
    if kind not in ("flow", "map"):
        raise ValueError(f"kind must be 'flow' or 'map', got {kind!r}")
    if kind == "flow" and not dt:
        raise ValueError("flows need an integration step dt")

    x0 = np.asarray(x0, dtype=np.float64).ravel()
    D = x0.size
    k = D if k is None else int(k)
    if not 1 <= k <= D:
        raise ValueError(f"k must be in 1..{D}, got {k}")

    E = int(ensemble)
    rng = np.random.default_rng(seed)

    # Member 1 is the published initial condition verbatim; the rest are
    # perturbed around it and land on the attractor during the transient.
    fid = np.repeat(x0[None, :], E, axis=0)
    if E > 1:
        scale = np.where(np.abs(x0) > 0, np.abs(x0), 1.0)
        fid[1:] += rng.normal(size=(E - 1, D)) * spread * scale

    step = (lambda t, y: _rk4(f, t, y, dt)) if kind == "flow" else (lambda t, y: f(t, y))
    h = dt if kind == "flow" else 1.0

    # Relax the fiducial ensemble onto the attractor before any perturbed copy
    # is created: a separation grown during the transient measures the approach
    # to the attractor, not the dynamics on it.
    t = t0
    seen = []
    for _j in range(int(transient_steps)):
        fid = step(t, fid)
        t += h
        if _j % 100 == 0:
            seen.append(fid.copy())

    # d0 IS RELATIVE TO THE ATTRACTOR, NOT ABSOLUTE. The separation has to be
    # small enough to stay linear and large enough to survive the subtraction
    # of two nearby states. Those are both statements about the state's own
    # magnitude, not about 1e-8. The LCG makes the point: its state lives near
    # 1e5, so an absolute d0 of 1e-9 is 1e-14 in relative terms, a few bits
    # above machine epsilon, and the exponent came out 7.4e-4 below ln(7141)
    # with a standard error of 3e-7 -- a bias, not noise. Scaling d0 by the
    # spread of the attractor makes the same argument hold for every system.
    if seen:
        _S = np.stack(seen, axis=0)
        _scale = np.linalg.norm(_S.std(axis=0), axis=-1)
        _scale = np.where(np.isfinite(_scale) & (_scale > 0), _scale, 1.0)
        d0 = d0 * _scale[:, None]
    else:
        d0 = np.full((E, 1), float(d0))

    # Separation directions start orthonormal. They are not yet aligned with
    # the Lyapunov directions, which is why a second, shorter alignment phase
    # runs before accumulation begins.
    q = np.repeat(np.eye(D, dtype=np.float64)[None, :k, :], E, axis=0)
    pert = fid[:, None, :] + d0[:, :, None] * q

    # Running total of div f (flows) or log|det J| (maps) along the fiducial,
    # and how many samples went into it, for the independent check on sum(lambda).
    div_sum = np.zeros(E, dtype=np.float64)
    div_n = 0
    div_stride = max(1, 100 // renorm_every)

    # THE LINEAR REGIME, AND HOW IT FAILS SILENTLY. The method assumes the
    # separation stays infinitesimal, so that its growth measures the
    # linearized dynamics. One event breaks that assumption completely: if the
    # fiducial trajectory lands exactly on an unstable fixed point -- which
    # happens in double precision on maps with a critical point, the logistic
    # map at r = 4 reaching exactly 0 after a few hundred thousand iterations
    # -- then the fiducial stops moving while its perturbed copy escapes to the
    # attractor. The separation then saturates at the size of the attractor and
    # is renormalized from d0 every step, so the member accumulates
    # log(attractor/d0) per iteration, about 20 nats, and reports an exponent
    # in the tens. The mean over the ensemble is not robust to even one such
    # member, while the median is untouched; that asymmetry is how this was
    # found, and it is not something to leave to a robust statistic.
    #
    # A cap on growth alone cannot catch it, because legitimate growth can be
    # large: the LCG stretches by 7141 every iteration. What distinguishes the
    # failure is that the separation approaches the size of the attractor
    # itself. So the state scale is measured during alignment, and any
    # renormalization where the separation exceeds a small fraction of it is
    # counted as saturated.
    # Renormalization intervals that actually contributed to each exponent.
    # Intervals where the separation underflowed are excluded, so the average
    # is over what was measured rather than over wall-clock.
    contributed = np.zeros((E, k), dtype=np.int64)
    sat_count = np.zeros(E, dtype=np.int64)

    # WHEN THE FIDUCIAL TRAJECTORY ITSELF DIES. Distinct from a collapsed
    # separation, and not fixable by re-seeding. In double precision a map's
    # orbit can land exactly on a fixed point -- the logistic map at r = 4
    # reaching exactly 0, the Gauss map reaching exactly 0 and dividing by it
    # -- after which the member is no longer sampling the attractor. Its
    # separation still evolves, perfectly linearly, about the fixed point, so
    # it reports that point's multiplier (ln 4 for the logistic map) instead of
    # the exponent, and nothing about the separation looks wrong.
    #
    # Sprott documents this for the tent map and the Kaplan-Yorke map, where it
    # happens within ~52 iterations and disqualifies the system outright. What
    # the long runs show is that it is not confined to those two: it is a
    # question of horizon, and at a few million iterations it reaches maps that
    # are perfectly well behaved at a few hundred thousand. A member that dies
    # is frozen at its death, keeping the shorter time average it had earned,
    # which is a valid estimate over a shorter record rather than a corrupted
    # one over a long record.
    member_alive = np.ones(E, dtype=bool)
    death_step = np.full(E, -1, dtype=np.int64)
    sat_scale = np.full(E, np.inf, dtype=np.float64)
    # Half the attractor's own size. An earlier, far tighter threshold (1% of
    # it) disqualified every member of the Gauss map, whose derivative -1/x^2
    # is unbounded near zero: a single legitimate iteration there stretches the
    # separation by 1e8, and over two million steps every member meets one.
    # That growth is the map, not a failure, so the threshold marks only a
    # separation comparable to the attractor itself, and a member is rejected
    # only if it is chronically outside the linear regime rather than for
    # isolated excursions.
    SAT_FRAC = 0.5
    SAT_TOLERATED = 0.10

    def advance(n_steps_local, t_local, fid_local, pert_local, accum,
                sample_into=None):
        """Run n steps, renormalizing on schedule. accum may be None."""
        nonlocal div_sum, div_n, sat_count, contributed
        nonlocal member_alive, death_step
        prev_fid = fid_local
        for i in range(int(n_steps_local)):
            prev_fid = fid_local
            fid_local = step(t_local, fid_local)
            # The perturbed copies see the same t: they are trajectories of the
            # same nonautonomous system, not of a shifted one.
            pert_local = step(t_local, pert_local)
            t_local += h
            if (i + 1) % renorm_every == 0:
                v = _separation(pert_local, fid_local, wrap)
                qq, norms = _gram_schmidt(v)
                if sample_into is not None:
                    sample_into.append(fid_local.copy())
                elif accum is not None:
                    sat_count += (norms[:, 0] > SAT_FRAC * sat_scale)
                    # A fiducial that has stopped moving is at a fixed point of
                    # the floating-point map and will never leave it; one that
                    # has gone non-finite has left the attractor entirely.
                    died = member_alive & (
                        ~np.isfinite(fid_local).all(axis=-1)
                        | np.all(fid_local == prev_fid, axis=-1))
                    if died.any():
                        death_step[died] = i
                        member_alive &= ~died
                if accum is not None:
                    # Dead members contribute nothing further; their average is
                    # frozen at the record they earned while alive.
                    live_n = (norms > 0) & member_alive[:, None]
                    accum += np.where(live_n,
                                       np.log(np.where(live_n, norms, 1.0) / d0),
                                       0.0)
                    contributed += live_n
                    # The divergence is a time average, so it is subsampled:
                    # maps renormalize every step and evaluating a determinant
                    # there would cost more than the exponents themselves.
                    if divergence and (i // renorm_every) % div_stride == 0:
                        g = (_divergence(f, t_local, fid_local) if kind == "flow"
                             else _log_abs_det_jac(f, t_local, fid_local))
                        div_sum += np.where(np.isfinite(g), g, 0.0)
                        div_n += 1
                pert_local = fid_local[:, None, :] + d0[:, :, None] * qq
        return t_local, fid_local, pert_local, accum

    # Alignment phase: renormalize without accumulating, so the vectors reach
    # the Lyapunov directions before they start contributing to the average.
    align = max(1, int(transient_steps))
    scale_samples = []
    t, fid, pert, _ = advance(align, t, fid, pert, None,
                              sample_into=scale_samples)
    if scale_samples:
        # Spread of the fiducial over the alignment phase, per member: the
        # size of the set the trajectory actually visits.
        S_ = np.stack(scale_samples, axis=0)                  # (n, E, D)
        sat_scale = np.linalg.norm(S_.std(axis=0), axis=-1)
        sat_scale = np.where(sat_scale > 0, sat_scale, np.inf)

    accum = np.zeros((E, k), dtype=np.float64)
    n_steps = int(n_steps)

    if checkpoints and checkpoints > 1:
        marks = np.unique(np.geomspace(max(renorm_every, n_steps // 1000),
                                       n_steps, checkpoints).astype(np.int64))
        history, cptime = [], []
        done = 0
        for m in marks:
            t, fid, pert, accum = advance(m - done, t, fid, pert, accum)
            done = int(m)
            history.append((accum / (done * h)).mean(axis=0).copy())
            cptime.append(done * h)
        history = np.asarray(history)
        cptime = np.asarray(cptime)
        n_steps = done
    else:
        t, fid, pert, accum = advance(n_steps, t, fid, pert, accum)
        history = cptime = None

    n_renorm_total = max(1, n_steps // renorm_every)
    n_died = int((death_step >= 0).sum())
    # A member that died early earned too short a record to average. Keeping it
    # at face value is how a frozen member with no record at all entered the
    # mean as a zero; requiring a minimum share of the run makes the cut
    # explicit instead.
    MIN_RECORD = 0.01
    enough = contributed[:, 0] >= MIN_RECORD * n_renorm_total
    # Intervals lost to a collapsed separation, not counting those a dead
    # member never reached.
    reachable = np.where(death_step >= 0, death_step // max(1, renorm_every),
                         n_renorm_total)
    n_collapsed = int(reachable.sum() * k - contributed.sum())
    # Divide by the time actually measured, not the time elapsed.
    measured = np.maximum(contributed, 1) * renorm_every * h
    per_member = accum / measured

    # A member whose trajectory left the attractor and overflowed carries no
    # exponent. Dropping it silently would report an average over whichever
    # members happened to survive, so the count is returned and the caller
    # decides whether the run is admissible.
    n_renorm = max(1, n_steps // renorm_every)
    saturated = sat_count > max(1, int(SAT_TOLERATED * n_renorm))
    alive = np.isfinite(per_member).all(axis=1) & ~saturated & enough
    n_escaped = int((~alive).sum())
    n_saturated = int(saturated.sum())
    if not alive.any():
        raise FloatingPointError("every ensemble member left the attractor")
    per_member = per_member[alive]

    exponents = per_member.mean(axis=0)
    n_alive = per_member.shape[0]
    sem = (per_member.std(axis=0, ddof=1) / np.sqrt(n_alive)
           if n_alive > 1 else np.zeros(k))

    div = np.nan
    if divergence and div_n:
        dm = (div_sum / div_n)[alive]
        div = float(dm.mean())

    return SpectrumResult(exponents, per_member, sem, history, cptime,
                          n_steps, n_steps * h, d0, renorm_every,
                          divergence=div, n_escaped=n_escaped,
                          n_saturated=n_saturated, n_collapsed=n_collapsed,
                          n_died=n_died)


def spectrum_adaptive(f, x0, *, max_halvings=6, min_alive_frac=0.5, **kw):
    """spectrum(), backing off the horizon until the orbits survive it.

    Maps have a maximum usable iteration count in double precision, and it is
    a property of the system, not a constant. Sprott disqualifies the tent and
    Kaplan-Yorke maps because their orbits reach exactly 0 within ~52 steps;
    the same death reaches the logistic map near a million iterations and the
    Gauss map, which divides by its state, sooner. Asking for a horizon longer
    than a system can support does not produce a less precise answer, it
    produces a wrong one or no answer at all, so the horizon is discovered
    rather than assumed: halve it until most members survive, and report what
    was actually used.

    Returns (result, n_steps_used). The caller records n_steps_used, because a
    system that needed a shorter horizon has a correspondingly larger sem and
    that is a fact about the system worth keeping.
    """
    n = int(kw.pop("n_steps"))
    tr = int(kw.pop("transient_steps", 10_000))
    ensemble = int(kw.get("ensemble", 1))
    last = None
    for _ in range(max_halvings + 1):
        try:
            r = spectrum(f, x0, n_steps=n,
                         transient_steps=min(tr, max(1000, n // 10)), **kw)
        except FloatingPointError as exc:
            last = exc
            n //= 2
            if n < 1000:
                break
            continue
        if r.n_died <= (1 - min_alive_frac) * ensemble:
            return r, n
        last = r
        n //= 2
        if n < 1000:
            break
    if isinstance(last, BaseException):
        raise last
    raise FloatingPointError(
        "no horizon short enough kept the orbits alive")
