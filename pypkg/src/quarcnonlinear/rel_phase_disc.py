"""Discrete relative phase, frequency ratios and phase coupling."""

# Copyright (c) 2021-2026 Quantitative Analysis Research Core,
# Center for Human Movement Variability, University of Nebraska at Omaha.
# MIT licence. See LICENSE.txt.

import numpy as np

from ._mround import mround


def _put(v, i, val):
    """v[i] = val with 1-based i, growing v with zeros past its end as
    MATLAB's assignment does."""
    if i > len(v):
        v.extend([0.0] * (int(i) - len(v)))
    v[int(i) - 1] = val


def rel_phase_disc(x, y, S1_events=None, S2_events=None):
    """Event-based relative phase between two oscillators, the frequency
    ratios it implies and the phase coupling at the dominant ratio, as
    MATLAB's rel_phase_disc without its figures.

    S1_events and S2_events are 0-based event indices for x and y; by default
    their peaks, found as MATLAB's findpeaks does.

    Returns (FR, PC, DRP), each a dict of 1-D arrays: FR with Ratio,
    Occurances and Percentage_of_Total; PC with n and PC; DRP with t, n, T
    and DRP (one row per y event in a complete x cycle, DRP in degrees).
    """
    S1 = np.asarray(_find_peaks(x) if S1_events is None else S1_events, dtype=np.float64).ravel()
    S2 = np.asarray(_find_peaks(y) if S2_events is None else S2_events, dtype=np.float64).ravel()
    if S1.size < 2:
        raise ValueError("x needs at least two events.")

    # One row per y event j that lies in a complete x cycle and has a next y
    # event: t is the time from it to the next x event, T the y period
    # starting at it, and n counts the y events still to come in the same x
    # cycle (0 for the cycle's last).
    head = S2[:-1]
    use = np.flatnonzero((head >= S1[0]) & (head < S1[-1]))
    cyc = np.array([np.flatnonzero(S1 <= S2[j]).max() for j in use], dtype=np.int64)
    tcol = S1[cyc + 1] - S2[use]
    Tcol = S2[use + 1] - S2[use]
    n = np.array([np.sum(cyc[k + 1 :] == cyc[k]) for k in range(use.size)], dtype=np.float64)
    DRP = tcol / Tcol * 360
    Raw = {"t": tcol, "n": n, "T": Tcol, "DRP": DRP}
    maxDRP = int(mround(DRP.max() / 360)) + 1

    DRPfirst = [DRP[0]]
    nrow = DRP.size
    for k in range(3, nrow):                  # MATLAB's 3:(nrow-1), 1-based
        if n[k - 1] == 0:
            DRPfirst.append(DRP[k])
    freq = []
    for k, dk in enumerate(DRPfirst, start=1):
        for i1 in range(1, maxDRP + 1):
            if (i1 - 1) * 360 < dk <= i1 * 360:
                _put(freq, k, float(i1))
    if len(freq) < 4:
        raise ValueError("Too few cycles to classify frequency ratios.")
    half = []
    if freq[0] == freq[2] and freq[1] == freq[3] and abs(freq[0] - freq[1]) == 1:
        half.append((freq[0] + freq[1]) / 2)
    else:
        half.append(freq[0])
    last = len(freq)
    counter = 0
    for k in range(2, last):                  # 1-based k, as MATLAB
        if freq[k - 2] == freq[k] and abs(freq[k - 1] - freq[k]) == 1:
            _put(half, k, (freq[k - 1] + freq[k]) / 2)
            counter += 1
        else:
            if counter % 2 == 1 and counter < k and counter > 1:
                _put(half, k - counter - 1, half[k - counter - 1])
            _put(half, k, freq[k - 1])
            counter = 0
    _put(half, last, freq[last - 1])
    half = np.array(half)
    nh = half.size
    withNC = np.zeros(nh)
    withNC[0] = half[0] if half[0] == half[1] else 0
    for k in range(1, nh - 1):
        withNC[k] = half[k] if (half[k] == half[k + 1] or half[k] == half[k - 1]) else 0
    withNC[-1] = half[-1] if half[-1] == half[-2] else 0

    maxn = withNC.max()
    Ratio = [0.0]
    Occ = [float(np.sum(withNC == 0))]
    Pct = [Occ[0] / nh * 100]
    if maxn >= 1:
        for r in np.arange(1, maxn + 0.25, 0.5):
            i = int(round(r * 2))
            _put(Ratio, i, float(r))
            _put(Occ, i, float(np.sum(withNC == r)))
            _put(Pct, i, Occ[i - 1] / nh * 100)
    Ratio, Occ, Pct = np.array(Ratio), np.array(Occ), np.array(Pct)
    lag = int(np.argmax(Occ)) + 1
    if lag % 2 == 0:
        lag = lag // 2
    original = DRP[: DRP.size - lag]
    lagged = DRP[lag:]
    kmax = int(mround(Ratio[Pct == Pct.max()][0]))
    wdn = []
    nPC = np.zeros(max(kmax, 0))
    PC = np.zeros(max(kmax, 0))
    for k in range(1, kmax + 1):
        sel = (original <= 360 * k) & (lagged <= 360 * k)
        dn = (original[sel] - lagged[sel]) * np.cos(np.pi / 4)
        for q, dq in enumerate(dn, start=1):
            _put(wdn, q, 1 - abs(dq) / 40 if abs(dq) <= 40 else 0.0)
        nPC[k - 1] = kmax - k
        PC[k - 1] = np.sum(wdn) / len(wdn) * 100 if wdn else np.nan
    FR = {"Ratio": Ratio, "Occurances": Occ, "Percentage_of_Total": Pct}
    return FR, {"n": nPC, "PC": PC}, Raw


def _find_peaks(x):
    """0-based indices of local maxima, as MATLAB's findpeaks: strictly above
    the preceding sample, and for a flat top the first sample of the plateau,
    provided it is followed by a fall. Endpoints are never peaks."""
    x = np.asarray(x, dtype=np.float64).ravel()
    n = x.size
    peaks = []
    i = 1
    while i < n - 1:
        if x[i] > x[i - 1]:
            j = i
            while j < n - 1 and x[j + 1] == x[i]:
                j += 1
            if j < n - 1 and x[j + 1] < x[i]:
                peaks.append(i)
            i = j + 1
        else:
            i += 1
    return np.array(peaks, dtype=np.int64)
