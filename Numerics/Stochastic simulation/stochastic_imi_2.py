#!/usr/bin/env python3
"""
Stochastic interference mutual information for the eight feed-forward loop types.

For each motif this program runs two exact stochastic simulations of the
chemical master equation:

  (i)  the feed-forward loop,        X -> Y -> Z together with X -> Z,
  (ii) its open-loop equivalent, in which the indirect route is driven by an
       independent but statistically equivalent input Xt, so that the two
       routes no longer share an upstream regulator.

Second moments are accumulated by time averaging and converted to mutual
information through the Gaussian-channel expression.  Three quantities are
reported for every motif: the mutual information of the feed-forward loop, the
mutual information of its open-loop equivalent, and their difference

    I_int(X;Z) = I(X;Z)_FFL - I(X,Xt;Z)_open-loop

Each is written with its standard error across replicates.

Both simulations run inside this program, so it needs no input from any other
code and writes the final data files directly.

Run length is set in units of the input correlation time 1/beta_X, which gives
comparable statistical precision to parameter sets whose degradation rates
differ by more than an order of magnitude.  Because the interference term is
the difference of two mutual informations of similar size, motifs with a small
|I_int| need long runs; the standard-error files say whether a value is
resolved.

Usage
-----
    python3 stochastic_imi.py                 defaults below
    python3 stochastic_imi.py 50000 6         run length, replicate count

Only the standard library is required.  If numba happens to be installed the
two simulation loops are compiled automatically, which speeds them up by
roughly two orders of magnitude; nothing else changes.

Results are written under results/stochastic/ as two-column files,
"motif <TAB> value".
"""

import os
import sys
import random
from math import log, log2, sqrt

# ---------------------------------------------------------------------------
# Defaults, overridable on the command line.
# ---------------------------------------------------------------------------
DEFAULT_NTAU = 10000.0   # run length in units of 1/beta_X
DEFAULT_NREP = 10         # independent replicates per motif
BURN_IN_TAU = 20.0        # discarded transient, in units of 1/beta_X

# ---------------------------------------------------------------------------
# Optional just-in-time compilation.  Without numba the identical pure-Python
# functions are used.
# ---------------------------------------------------------------------------
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]

        def decorate(function):
            return function
        return decorate


# ---------------------------------------------------------------------------
# Regulatory topology of the eight FFL types.
# Each entry is (sign of X->Y, sign of X->Z, sign of Y->Z);
# +1 denotes activation and -1 denotes repression.
# ---------------------------------------------------------------------------
TOPOLOGY = {
    "C1": (+1, +1, +1),
    "C2": (-1, -1, +1),
    "C3": (+1, -1, -1),
    "C4": (-1, +1, -1),
    "I1": (+1, +1, -1),
    "I2": (-1, -1, -1),
    "I3": (+1, -1, +1),
    "I4": (-1, +1, +1),
}

COHERENT = ["C1", "C2", "C3", "C4"]
INCOHERENT = ["I1", "I2", "I3", "I4"]

# ---------------------------------------------------------------------------
# Optimized parameters.
# Order: beta_X, beta_Y, beta_Z [1/min],
#        K_XY, K_XZ, K_YZ, <x>, <y>, <z> [molecules/cell].
# The four motifs of a class share one parameter set.
# ---------------------------------------------------------------------------
OPTIMIZED = {
    ("Ecoli", "AND", "coherent"):
        (5.0007e-3, 9.9884e-2, 4.9986e-1, 63.52, 69.95, 18.40, 10.64, 13.17, 99.94),
    ("Ecoli", "AND", "incoherent"):
        (5.0023e-3, 9.8712e-2, 7.2570e-2, 54.32, 99.55, 37.29, 15.93, 48.82, 56.96),
    ("Ecoli", "OR", "coherent"):
        (1.3426e-2, 9.3002e-2, 2.7098e-1, 64.61, 41.42, 26.90, 10.00, 17.25, 99.95),
    ("Ecoli", "OR", "incoherent"):
        (3.4645e-2, 6.2942e-2, 2.6985e-1, 34.95, 44.20, 21.86, 10.12, 39.82, 99.83),
    ("Yeast", "AND", "coherent"):
        (1.2336e-2, 9.0400e-3, 2.8630e-2, 467.23, 27.69, 499.80, 50.00, 63.42, 499.96),
    ("Yeast", "AND", "incoherent"):
        (1.2494e-3, 6.2408e-2, 2.5360e-2, 217.30, 157.60, 25.69, 50.01, 446.66, 306.79),
    ("Yeast", "OR", "coherent"):
        (1.2230e-3, 7.6122e-2, 4.2529e-2, 134.93, 294.17, 401.08, 51.13, 53.04, 357.06),
    ("Yeast", "OR", "incoherent"):
        (1.0933e-3, 4.4021e-2, 8.8137e-2, 261.24, 272.19, 22.59, 112.80, 431.18, 206.83),
}

# ---------------------------------------------------------------------------
# Fixed representative parameter set, shared by all eight motifs and both gates.
# ---------------------------------------------------------------------------
REPRESENTATIVE = (0.1, 1.0, 10.0, 100.0, 100.0, 100.0, 40.0, 100.0, 100.0)


def hill_coefficients(sign, k):
    """
    Write the Hill regulatory factor as (a*v + b) / (k + v), so that activation
    and repression differ only in the constants a and b.  This keeps the
    simulation loops free of branching on the regulatory sign.
    """
    return (1.0, 0.0) if sign > 0 else (0.0, k)


# ---------------------------------------------------------------------------
# Feed-forward loop.  Returns I(X;Z) in bits.
# ---------------------------------------------------------------------------
@njit(cache=True)
def simulate_ffl(ya, yb, zxa, zxb, zya, zyb, or_gate,
                 bx, by, bz, kxy, kxz, kyz, ax, ay, az,
                 x0, y0, z0, total_time, burn_in, seed):
    random.seed(seed)

    x = x0
    y = y0
    z = z0
    t = 0.0
    sampled = 0.0
    sx = 0.0
    sz = 0.0
    sxx = 0.0
    szz = 0.0
    sxz = 0.0

    while t < total_time:
        xf = float(x)
        yf = float(y)
        hy = (ya * xf + yb) / (kxy + xf)
        hzx = (zxa * xf + zxb) / (kxz + xf)
        hzy = (zya * yf + zyb) / (kyz + yf)
        if or_gate:
            gz = hzx + hzy
        else:
            gz = hzx * hzy

        p0 = ax
        p1 = bx * xf
        p2 = ay * hy
        p3 = by * yf
        p4 = az * gz
        p5 = bz * float(z)
        total = p0 + p1 + p2 + p3 + p4 + p5
        if total <= 0.0:
            break

        dt = -log(random.random()) / total
        step = dt
        if step > total_time - t:
            step = total_time - t

        if t >= burn_in:
            weight = step
        elif t + step > burn_in:
            weight = t + step - burn_in
        else:
            weight = 0.0

        if weight > 0.0:
            zf = float(z)
            sx += xf * weight
            sz += zf * weight
            sxx += xf * xf * weight
            szz += zf * zf * weight
            sxz += xf * zf * weight
            sampled += weight

        t += dt
        if t >= total_time:
            break

        r = random.random() * total
        if r < p0:
            x += 1
        elif r < p0 + p1:
            x -= 1
        elif r < p0 + p1 + p2:
            y += 1
        elif r < p0 + p1 + p2 + p3:
            y -= 1
        elif r < p0 + p1 + p2 + p3 + p4:
            z += 1
        else:
            z -= 1

    mx = sx / sampled
    mz = sz / sampled
    ex = (sxx / sampled - mx * mx) / (mx * mx)
    ez = (szz / sampled - mz * mz) / (mz * mz)
    cxz = (sxz / sampled - mx * mz) / (mx * mz)

    return 0.5 * log2(ex * ez / (ex * ez - cxz * cxz))


# ---------------------------------------------------------------------------
# Open-loop equivalent.  X drives the direct route, Xt the indirect route.
# Returns I(X,Xt;Z) in bits.
# ---------------------------------------------------------------------------
@njit(cache=True)
def simulate_open_loop(ya, yb, zxa, zxb, zya, zyb, or_gate,
                       bx, by, bz, kxy, kxz, kyz, ax, ay, az,
                       x0, y0, z0, total_time, burn_in, seed):
    random.seed(seed)

    x = x0
    xt = x0
    y = y0
    z = z0
    t = 0.0
    sampled = 0.0
    sx = 0.0
    st = 0.0
    sz = 0.0
    sxx = 0.0
    stt = 0.0
    szz = 0.0
    sxz = 0.0
    stz = 0.0

    while t < total_time:
        xf = float(x)
        tf = float(xt)
        yf = float(y)
        hy = (ya * tf + yb) / (kxy + tf)
        hzx = (zxa * xf + zxb) / (kxz + xf)
        hzy = (zya * yf + zyb) / (kyz + yf)
        if or_gate:
            gz = hzx + hzy
        else:
            gz = hzx * hzy

        p0 = ax
        p1 = bx * xf
        p2 = ax
        p3 = bx * tf
        p4 = ay * hy
        p5 = by * yf
        p6 = az * gz
        p7 = bz * float(z)
        total = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7
        if total <= 0.0:
            break

        dt = -log(random.random()) / total
        step = dt
        if step > total_time - t:
            step = total_time - t

        if t >= burn_in:
            weight = step
        elif t + step > burn_in:
            weight = t + step - burn_in
        else:
            weight = 0.0

        if weight > 0.0:
            zf = float(z)
            sx += xf * weight
            st += tf * weight
            sz += zf * weight
            sxx += xf * xf * weight
            stt += tf * tf * weight
            szz += zf * zf * weight
            sxz += xf * zf * weight
            stz += tf * zf * weight
            sampled += weight

        t += dt
        if t >= total_time:
            break

        r = random.random() * total
        c = p0
        if r < c:
            x += 1
        else:
            c += p1
            if r < c:
                x -= 1
            else:
                c += p2
                if r < c:
                    xt += 1
                else:
                    c += p3
                    if r < c:
                        xt -= 1
                    else:
                        c += p4
                        if r < c:
                            y += 1
                        else:
                            c += p5
                            if r < c:
                                y -= 1
                            else:
                                c += p6
                                if r < c:
                                    z += 1
                                else:
                                    z -= 1

    mx = sx / sampled
    mt = st / sampled
    mz = sz / sampled
    ex = (sxx / sampled - mx * mx) / (mx * mx)
    et = (stt / sampled - mt * mt) / (mt * mt)
    ez = (szz / sampled - mz * mz) / (mz * mz)
    cxz = (sxz / sampled - mx * mz) / (mx * mz)
    ctz = (stz / sampled - mt * mz) / (mt * mz)

    conditional = ez - (cxz * cxz / ex + ctz * ctz / et)
    return 0.5 * log2(ez / conditional)


# ---------------------------------------------------------------------------
def motif_arguments(motif, gate, bx, by, bz, kxy, kxz, kyz, xav, yav, zav):
    """Assemble the constant arguments both simulation loops expect."""
    sign_y, sign_zx, sign_zy = TOPOLOGY[motif]
    ya, yb = hill_coefficients(sign_y, kxy)
    zxa, zxb = hill_coefficients(sign_zx, kxz)
    zya, zyb = hill_coefficients(sign_zy, kyz)
    or_gate = (gate == "OR")

    # Synthesis rate constants fixed by the mean-field steady state.
    hy = (ya * xav + yb) / (kxy + xav)
    hzx = (zxa * xav + zxb) / (kxz + xav)
    hzy = (zya * yav + zyb) / (kyz + yav)
    gz = hzx + hzy if or_gate else hzx * hzy

    ax = bx * xav
    ay = by * yav / hy
    az = bz * zav / gz

    return (ya, yb, zxa, zxb, zya, zyb, or_gate,
            bx, by, bz, kxy, kxz, kyz, ax, ay, az,
            int(round(xav)), int(round(yav)), int(round(zav)))


def write_column(directory, quantity, suffix, organism, gate, cls,
                 motifs, values):
    """Write a two-column 'motif <TAB> value' data file."""
    if organism:
        name = "%s-%s-%s-%s%s.dat" % (quantity, cls, gate, organism, suffix)
    else:
        name = "%s-%s-%s%s.dat" % (quantity, cls, gate, suffix)
    with open(os.path.join(directory, name), "w") as handle:
        for motif, value in zip(motifs, values):
            handle.write("%s\t%.6f\n" % (motif, value))


def run_set(params, organism, gate, cls, directory, n_tau, n_rep, salt):
    motifs = COHERENT if cls == "coherent" else INCOHERENT
    bx = params[0]
    total_time = n_tau / bx
    burn_in = BURN_IN_TAU / bx

    imi, sem = [], []
    ffl_mi, ffl_sem = [], []
    open_mi, open_sem = [], []

    def mean_and_error(values):
        """Sample mean and standard error of the mean."""
        average = sum(values) / len(values)
        if len(values) > 1:
            variance = sum((v - average) ** 2 for v in values) / (len(values) - 1)
            return average, sqrt(variance / len(values))
        return average, 0.0

    for index, motif in enumerate(motifs):
        constants = motif_arguments(motif, gate, *params)

        differences, ffl_values, open_values = [], [], []
        for replicate in range(n_rep):
            seed = salt + 1000 * index + 2 * replicate
            a = simulate_ffl(*constants, total_time, burn_in, seed)
            b = simulate_open_loop(*constants, total_time, burn_in, seed + 1)
            differences.append(a - b)
            ffl_values.append(a)
            open_values.append(b)

        mean, error = mean_and_error(differences)
        ffl_mean, ffl_error = mean_and_error(ffl_values)
        open_mean, open_error = mean_and_error(open_values)

        imi.append(mean)
        sem.append(error)
        ffl_mi.append(ffl_mean)
        ffl_sem.append(ffl_error)
        open_mi.append(open_mean)
        open_sem.append(open_error)

        print("  %-3s  IMI = %+.6f +/- %.6f   I_FFL = %.6f +/- %.6f   "
              "I_open = %.6f +/- %.6f"
              % (motif, mean, error, ffl_mean, ffl_error, open_mean, open_error))
        sys.stdout.flush()

    write_column(directory, "IMI", "", organism, gate, cls, motifs, imi)
    write_column(directory, "IMI", "-sem", organism, gate, cls, motifs, sem)
    write_column(directory, "Itotal", "", organism, gate, cls, motifs, ffl_mi)
    write_column(directory, "Itotal", "-sem", organism, gate, cls, motifs, ffl_sem)
    write_column(directory, "Ipath", "", organism, gate, cls, motifs, open_mi)
    write_column(directory, "Ipath", "-sem", organism, gate, cls, motifs, open_sem)


def main():
    n_tau = float(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_NTAU
    n_rep = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_NREP
    n_rep = max(1, n_rep)

    optimized_dir = os.path.join("results", "stochastic", "optimized")
    representative_dir = os.path.join("results", "stochastic", "representative")
    os.makedirs(optimized_dir, exist_ok=True)
    os.makedirs(representative_dir, exist_ok=True)

    print("run length = %.0f input correlation times, %d replicates"
          % (n_tau, n_rep))
    print("just-in-time compilation: %s\n"
          % ("enabled" if HAVE_NUMBA else "not available, using pure Python"))

    for index, ((organism, gate, cls), params) in enumerate(sorted(OPTIMIZED.items())):
        print("optimized  %s  %s  %s" % (organism, gate, cls))
        run_set(params, organism, gate, cls, optimized_dir,
                n_tau, n_rep, 7919 * (index + 1))

    index = 0
    for gate in ("AND", "OR"):
        for cls in ("coherent", "incoherent"):
            print("representative  %s  %s" % (gate, cls))
            run_set(REPRESENTATIVE, None, gate, cls, representative_dir,
                    n_tau, n_rep, 104729 * (index + 1))
            index += 1

    print("\nWritten under results/stochastic/")


if __name__ == "__main__":
    main()
