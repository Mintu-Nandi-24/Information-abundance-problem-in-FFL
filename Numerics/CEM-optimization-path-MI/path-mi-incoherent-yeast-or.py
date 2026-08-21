# =============================================================================
# FILE: Numerics/CEM-optimization-path-MI/path-mi-incoherent-yeast-or.py
# PAPER: An information-theoretic perspective on feed-forward loop abundances in transcriptional networks
# AUTHORS: Mintu Nandi, Sudip Chattopadhyay, and Suman K Banik
# CONTACT: mintunandi@ubi.s.u-tokyo.ac.jp; sudip@chem.iiests.ac.in;
#          skbanik@jcbose.ac.in
#
# PURPOSE
#   Performs abundance-guided cross-entropy optimization of pathway mutual information
#   ($I_{\mathrm{path}}(X;Z)$) for incoherent FFLs in *S. cerevisiae* / yeast with an OR
#   output gate. The four motifs (I1, I2, I3, I4) share one nine-parameter
#   operating-point vector in each candidate evaluation.
#
# OPTIMIZATION SETUP
#   Target windows: I1/I2 in [5, 7], I1/I3 in [18, 22], I1/I4 in [50, 100].
#   Default budget: 1000 pre-scan points; 100 runs;
#   1000 CEM iterations/run; population 1000;
#   elite fraction 0.1; smoothing 0.3; seed 42.
#   The objective first minimizes ratio-constraint violations and then
#   maximizes the mean motif metric within the feasible region.
#
# INPUTS AND EXECUTION
#   No external input files or command-line arguments are used.
#   Bounds, targets, seeds, and plotting settings are defined below.
#   Recommended working directory: Numerics/CEM-optimization-path-MI/
#   Run: python path-mi-incoherent-yeast-or.py
#
# OUTPUTS
#   Relative output directory: incoherent_yeast_OR_data/
#   Files: incoherent-yeast-OR-diagnostics.dat, incoherent-yeast-OR-best_penalty_convergence.dat, incoherent-yeast-OR-imi_trajectory.dat, incoherent-yeast-OR-final_imi_summary.dat, incoherent-yeast-OR-ratio_satisfaction_scatter.dat, incoherent-yeast-OR-avg_imi_by_run.dat, incoherent-yeast-OR-pairwise_distances.dat, incoherent-yeast-OR-parameter_correlations.dat, incoherent-yeast-OR-normalized_sensitivity.dat, incoherent-yeast-OR-parameter_cv.dat, incoherent-yeast-OR-violin_data.dat, incoherent-yeast-OR-final_parameter_summary.dat
#   Matplotlib figures are displayed interactively and are not saved.
#
# CODE-TO-MANUSCRIPT NOTATION
#   theta = (beta_x,beta_y,beta_z,Kxy,Kxz,Kyz,xbar,ybar,zbar)
#         = (beta_X,beta_Y,beta_Z,K_XY,K_XZ,K_YZ,<x>,<y>,<z>).
#   (s_xy,s_xz,s_yz) gives the signs of X->Y, X->Z, and Y->Z.
#   fyxp,fzxp,fzyp = f'_YX,f'_ZX,f'_ZY.
#   eta_xz1,eta_xz2 = zeta_XZ,d and zeta_XZ,ind.
#   eta_zi,eta_zd,eta_zind,eta_zsyn,eta_p = eta_Z,0^2, eta_Z,d^2,
#       eta_Z,ind^2, eta_Z,int^2, and eta_Z,path^2.
#   denom_inner = eta_Z|X,path^2; numer_inner = eta_Z|X,int^2.
# =============================================================================

from __future__ import annotations
import math
import os
import warnings
from dataclasses import dataclass
from typing import Dict, Tuple, List, Any, Optional

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings('ignore')

# ============================================================
# 0) GLOBAL PARAMETERS - TUNE ALL SETTINGS HERE
# ============================================================

# Optimization Parameters
N_RUNS = 100  # 1000
N_SAMPLES_FEASIBILITY = 1000  # 200000
N_ITER_CEM = 1000  # 1000
POPULATION_SIZE = 1000  # 10000
ELITE_FRAC = 0.1
ALPHA = 0.3
COV_FLOOR = 1e-6

# Objective Parameters
INFEASIBLE_OFFSET = 1e6
FEASIBILITY_TOL = 1e-12

# Sensitivity Analysis Parameters
N_PERTURBATIONS = 200
PERTURBATION_FRAC = 0.1

# Visualization Parameters
FIG_SIZE_DIAGNOSTICS = (20, 22)
FIG_SIZE_STANDARD = (10, 6)
FIG_SIZE_IMI_TRAJECTORY = (14, 10)
FIG_SIZE_VIOLIN = (16, 6)
FIG_SIZE_SENSITIVITY = (14, 7)
FIG_SIZE_CV = (12, 6)

# Random Seed
RANDOM_SEED = 42

# Data Export Directory
EXPORT_DIR = "incoherent_yeast_OR_data"

# ============================================================
# 0b) MOTIFS AND CONSTANTS (DO NOT MODIFY)
# ============================================================
INCOHERENT_MOTIFS = {
    "I1": (+1, +1, -1),
    "I2": (-1, -1, -1),
    "I3": (+1, -1, +1),
    "I4": (-1, +1, +1),
}
MOTIF_NAMES = list(INCOHERENT_MOTIFS.keys())

PARAM_NAMES = ["beta_x", "beta_y", "beta_z", "Kxy", "Kxz", "Kyz", "xbar", "ybar", "zbar"]

# Desired ordering windows for incoherent FFLs
IFFL_TARGETS = {
    ("I1", "I2"): (5.0, 7.0),
    ("I1", "I3"): (18.0, 22.0),
    ("I1", "I4"): (50.0, 100.0),
}

# Parameter bounds
LO_THETA = np.array([1e-3, 1e-3, 1e-3, 20.0, 20.0, 20.0, 50.0, 50.0, 50.0], float)
HI_THETA = np.array([1e-1, 1e-1, 1e-1, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0], float)
LO_PHI, HI_PHI = np.log(LO_THETA), np.log(HI_THETA)

if not os.path.exists(EXPORT_DIR):
    os.makedirs(EXPORT_DIR)

# ============================================================
# 1) DATA EXPORT
# ============================================================
def export_to_dat(filename: str, data: np.ndarray, header: List[str], comments: str = ""):
    """Export data to .dat file compatible with OriginPro."""
    filepath = os.path.join(EXPORT_DIR, filename)
    with open(filepath, 'w') as f:
        if comments:
            f.write(f"# {comments}\n")
        f.write("# " + "\t".join(header) + "\n")
        for row in data:
            row_strs = []
            for val in row:
                if isinstance(val, str):
                    row_strs.append(val)
                elif isinstance(val, (int, np.integer)):
                    row_strs.append(str(int(val)))
                elif isinstance(val, (float, np.floating)):
                    row_strs.append(f"{float(val):.8e}")
                else:
                    row_strs.append(str(val))
            f.write("\t".join(row_strs) + "\n")
    print(f"  Exported: {filename}")

# ============================================================
# 2) REGULATION FUNCTIONS
# ============================================================
def hill(u: float, K: float, sign: int) -> float:
    if sign == +1:
        return u / (K + u)
    else:
        return K / (K + u)


def hill_prime(u: float, K: float, sign: int) -> float:
    base = K / ((K + u) ** 2)
    return base if sign == +1 else -base

# ============================================================
# 3) ANALYTICAL IMI CALCULATION FOR INCOHERENT FFL
# ============================================================
def imi_analytical_iffl(theta: np.ndarray, motif: str) -> Tuple[float, bool]:
    beta_x, beta_y, beta_z, Kxy, Kxz, Kyz, xbar, ybar, zbar = map(float, theta)

    if min(beta_x, beta_y, beta_z, Kxy, Kxz, Kyz, xbar, ybar, zbar) <= 0:
        return float("nan"), False

    s_xy, s_xz, s_yz = INCOHERENT_MOTIFS[motif]

    fyx = hill(xbar, Kxy, s_xy)
    fzx = hill(xbar, Kxz, s_xz)
    fzy = hill(ybar, Kyz, s_yz)

    eps = 1e-15
    if fyx <= eps or fzx <= eps or fzy <= eps:
        return float("nan"), False

    alpha_y = (beta_y * ybar) / fyx
    alpha_z = (beta_z * zbar) / (fzx + fzy)

    fyxp = alpha_y * hill_prime(xbar, Kxy, s_xy)
    fzxp = alpha_z * hill_prime(xbar, Kxz, s_xz)
    fzyp = alpha_z * hill_prime(ybar, Kyz, s_yz)

    Gxy = 1.0 / (beta_x + beta_y)
    Gxz = 1.0 / (beta_x + beta_z)
    Gyz = 1.0 / (beta_y + beta_z)
    Gxyz1 = 1.0 / ((beta_x + beta_y) * (beta_x + beta_z))
    Gxyz2 = (beta_x + beta_y + beta_z) / (beta_y * (beta_x + beta_y) * (beta_x + beta_z) * (beta_y + beta_z))
    Gxyz3 = (2.0 * beta_x + beta_y + beta_z) / ((beta_x + beta_y) * (beta_x + beta_z) * (beta_y + beta_z))

    eta_xy = (Gxy / ybar) * fyxp
    eta_xz1 = (Gxz / zbar) * fzxp
    eta_xz2 = (Gxyz1 / zbar) * fyxp * fzyp
    eta_xz = eta_xz1 + eta_xz2
    
    eta_yz1 = (Gyz / zbar) * fzyp
    eta_yz2 = (Gxyz2 * xbar) / (ybar * zbar) * (fyxp ** 2) * fzyp
    eta_yz3 = (Gxyz3 * xbar) / (ybar * zbar) * fyxp * fzxp

    eta_x = 1.0 / xbar

    eta_zi = 1.0 / zbar
    eta_zd = (fzxp * xbar) / (beta_z * zbar) * eta_xz1
    eta_zind = (fzyp * ybar) / (beta_z * zbar) * (eta_yz1 + eta_yz2)

    pref1 = (fzxp * xbar) / (beta_z * zbar)
    eta_zsyn1 = pref1 * eta_xz2
    pref2 = (fzyp * ybar) / (beta_z * zbar)
    eta_zsyn2 = pref2 * eta_yz3
    eta_zsyn = eta_zsyn1 + eta_zsyn2

    eta_p = eta_zi + eta_zd + eta_zind

    eta_z = eta_p + eta_zsyn

    min_denom = 1e-14
    min_arg = 1e-14
    

    denom_inner = eta_p - (eta_xz1 ** 2 + eta_xz2 ** 2) / max(eta_x, eps)
    if denom_inner <= min_denom:
        return float("nan"), False

    numer_inner = eta_zsyn - (2.0 * eta_xz1 * eta_xz2) / max(eta_x, eps)
    term2 = 1.0 + numer_inner / denom_inner
    if term2 <= min_arg:
        return float("nan"), False

    term1 = 1.0 + eta_zsyn / max(eta_p, eps)
    if term1 <= min_arg:
        return float("nan"), False

    ratio = term1 / term2
    if ratio <= min_arg:
        return float("nan"), False

    ixzs = 0.5 * np.log2(ratio)
    
    
    ratio_p = eta_p / denom_inner
    ixzp = 0.5 * np.log2(ratio_p)
    
    
    denom_total = eta_x * eta_z - eta_xz ** 2
    ratio_total = eta_x * eta_z / denom_total
    ixz = 0.5 * np.log2(ratio_total)
    
    return float(ixzp), True



def calculate_all_iffl_imis(theta: np.ndarray) -> Tuple[Dict[str, float], bool]:
    imis: Dict[str, float] = {}
    ok = True
    for motif in INCOHERENT_MOTIFS:
        value, valid = imi_analytical_iffl(theta, motif)
        imis[motif] = value
        ok = ok and valid and np.isfinite(value)
    return imis, ok


def compute_iffl_ratios(imis: Dict[str, float]) -> Dict[str, float]:
    ratios: Dict[str, float] = {}
    for (a, b), _ in IFFL_TARGETS.items():
        ratios[f"|{a}|/|{b}|"] = abs(imis[a]) / max(abs(imis[b]), 1e-9)
    return ratios


def ratio_penalty_from_imis(imis: Dict[str, float]) -> float:
    penalty = 0.0
    for (a, b), (lo, hi) in IFFL_TARGETS.items():
        r = abs(imis[a]) / max(abs(imis[b]), 1e-9)
        if r <= 0:
            return float("inf")
        lr = math.log(r)
        if lr < math.log(lo):
            penalty += (math.log(lo) - lr) ** 2
        elif lr > math.log(hi):
            penalty += (lr - math.log(hi)) ** 2
    return penalty


def evaluate_theta(theta: np.ndarray) -> Dict[str, Any]:
    """
    Objective for INCOHERENT FFL:
      1) IMI values must be negative (incoherent constraint)
      2) Preserve required ratio ordering (feasibility first)
      3) Among feasible solutions, maximize the average absolute IMI magnitude
    """
    imis, ok = calculate_all_iffl_imis(theta)

    if not ok:
        return {
            'loss': INFEASIBLE_OFFSET + 1e3,
            'penalty': float('inf'),
            'avg_imi': -float('inf'),
            'feasible': False,
            'imis': imis,
            'ratios': {},
        }

    ratios = compute_iffl_ratios(imis)
    penalty = ratio_penalty_from_imis(imis)
    # For incoherent, we want to maximize the average absolute IMI magnitude
    avg_imi = float(np.mean([abs(imis[m]) for m in INCOHERENT_MOTIFS]))
    feasible = bool(penalty <= FEASIBILITY_TOL)

    if feasible:
        loss = -avg_imi  # maximize average absolute IMI
    else:
        loss = INFEASIBLE_OFFSET + penalty

    return {
        'loss': float(loss),
        'penalty': float(penalty),
        'avg_imi': float(avg_imi),
        'feasible': feasible,
        'imis': imis,
        'ratios': ratios,
    }

# ============================================================
# 4) NORMALIZED SENSITIVITY ANALYSIS
# ============================================================
def normalized_sensitivity_analysis(
    theta: np.ndarray,
    param_names: List[str],
    n_perturbations: int = N_PERTURBATIONS,
    perturbation_frac: float = PERTURBATION_FRAC,
) -> Dict[str, Any]:
    """
    Normalized sensitivity (elasticity) for each motif:
        (ΔIMI/IMI) / (Δp/p)
    evaluated around a representative parameter vector.
    """
    original_imis, ok = calculate_all_iffl_imis(theta)
    if not ok:
        raise ValueError("Sensitivity analysis requires a valid representative parameter set.")

    results = {
        'param_name': param_names,
        'normalized_sensitivity': {p: 0.0 for p in param_names},
        'motif_sensitivity': {motif: {p: 0.0 for p in param_names} for motif in INCOHERENT_MOTIFS},
    }

    for i, param_name in enumerate(param_names):
        param_value = theta[i]
        for motif in INCOHERENT_MOTIFS:
            sensitivities = []
            for perturbation in np.linspace(-perturbation_frac, perturbation_frac, n_perturbations):
                if abs(perturbation) < 1e-12:
                    continue
                perturbed_theta = theta.copy()
                perturbed_theta[i] = max(LO_THETA[i], min(HI_THETA[i], param_value * (1.0 + perturbation)))
                perturbed_imis, ok = calculate_all_iffl_imis(perturbed_theta)
                if not ok:
                    continue
                if original_imis[motif] == 0:
                    continue
                delta_imi = perturbed_imis[motif] - original_imis[motif]
                norm_sens = (delta_imi / original_imis[motif]) / perturbation
                sensitivities.append(abs(norm_sens))
            avg_sens = float(np.mean(sensitivities)) if sensitivities else np.nan
            results['motif_sensitivity'][motif][param_name] = avg_sens
        across = [results['motif_sensitivity'][m][param_name] for m in INCOHERENT_MOTIFS]
        results['normalized_sensitivity'][param_name] = float(np.nanmean(across))

    return results

# ============================================================
# 5) CEM OPTIMIZER
# ============================================================
@dataclass
class CEMConfig:
    n_iter: int = N_ITER_CEM
    population_size: int = POPULATION_SIZE
    elite_frac: float = ELITE_FRAC
    alpha: float = ALPHA
    cov_floor: float = COV_FLOOR
    seed: int = RANDOM_SEED
    verbose: bool = False


@dataclass
class SingleRunResult:
    best_theta: np.ndarray
    best_loss: float
    best_penalty: float
    best_avg_imi: float
    best_feasible: bool
    best_ratios: Dict[str, float]
    best_imis: Dict[str, float]
    trace_means: np.ndarray
    trace_stds: np.ndarray
    trace_best_losses: np.ndarray
    trace_best_penalties: np.ndarray
    trace_best_avg_imi: np.ndarray
    trace_best_imis: np.ndarray
    trace_feasible_fraction: np.ndarray


class CEMOptimizer:
    def __init__(self, config: CEMConfig):
        self.config = config
        self.rng = np.random.default_rng(config.seed)

    def optimize(self, init_center: Optional[np.ndarray] = None) -> SingleRunResult:
        if init_center is None:
            mean_phi = (LO_PHI + HI_PHI) / 2
        else:
            mean_phi = np.log(init_center)

        init_std = (HI_PHI - LO_PHI) / 4
        cov_matrix = np.diag(init_std ** 2)

        trace_means = []
        trace_stds = []
        trace_best_losses = []
        trace_best_penalties = []
        trace_best_avg_imi = []
        trace_best_imis = []
        trace_feasible_fraction = []

        best_loss = float('inf')
        best_penalty = float('inf')
        best_avg_imi = -float('inf')
        best_theta = None
        best_ratios: Dict[str, float] = {}
        best_imis: Dict[str, float] = {}
        best_feasible = False

        for _ in range(self.config.n_iter):
            try:
                phi_samples = self.rng.multivariate_normal(mean_phi, cov_matrix, size=self.config.population_size)
            except np.linalg.LinAlgError:
                cov_matrix = np.diag(init_std ** 2)
                phi_samples = self.rng.multivariate_normal(mean_phi, cov_matrix, size=self.config.population_size)

            phi_samples = np.clip(phi_samples, LO_PHI, HI_PHI)
            theta_samples = np.exp(phi_samples)

            evals = [evaluate_theta(theta) for theta in theta_samples]
            losses = np.array([ev['loss'] for ev in evals], float)
            feasible_flags = np.array([ev['feasible'] for ev in evals], bool)

            sorted_indices = np.argsort(losses)
            elite_size = max(2, int(self.config.elite_frac * self.config.population_size))
            elite_indices = sorted_indices[:elite_size]
            elite_phis = phi_samples[elite_indices]

            current_best_idx = elite_indices[0]
            current_best_eval = evals[current_best_idx]
            current_best_loss = current_best_eval['loss']
            current_best_avg_imi = current_best_eval['avg_imi']

            if (current_best_loss < best_loss - 1e-15) or (
                abs(current_best_loss - best_loss) <= 1e-15 and current_best_avg_imi > best_avg_imi
            ):
                best_loss = current_best_loss
                best_penalty = current_best_eval['penalty']
                best_avg_imi = current_best_avg_imi
                best_theta = theta_samples[current_best_idx].copy()
                best_ratios = current_best_eval['ratios'].copy()
                best_imis = current_best_eval['imis'].copy()
                best_feasible = bool(current_best_eval['feasible'])

            new_mean = np.mean(elite_phis, axis=0)
            if elite_size > 1:
                X = elite_phis - new_mean
                new_cov = (X.T @ X) / (elite_size - 1)
            else:
                new_cov = np.diag(init_std ** 2)

            new_cov = np.atleast_2d(new_cov)
            if new_cov.shape != (9, 9):
                new_cov = np.diag(np.diag(new_cov) if new_cov.ndim == 2 else np.full(9, init_std[0] ** 2))

            mean_phi = (1 - self.config.alpha) * mean_phi + self.config.alpha * new_mean
            cov_matrix = (1 - self.config.alpha) * cov_matrix + self.config.alpha * new_cov
            cov_matrix += self.config.cov_floor * np.eye(9)
            cov_matrix = (cov_matrix + cov_matrix.T) / 2

            trace_means.append(mean_phi.copy())
            trace_stds.append(np.sqrt(np.diag(cov_matrix)).copy())
            trace_best_losses.append(best_loss)
            trace_best_penalties.append(best_penalty if np.isfinite(best_penalty) else np.nan)
            trace_best_avg_imi.append(best_avg_imi if np.isfinite(best_avg_imi) else np.nan)
            trace_best_imis.append([best_imis.get(m, np.nan) for m in MOTIF_NAMES])
            trace_feasible_fraction.append(np.mean(feasible_flags))

        if best_theta is None:
            best_theta = np.exp(mean_phi)
            fallback_eval = evaluate_theta(best_theta)
            best_loss = fallback_eval['loss']
            best_penalty = fallback_eval['penalty']
            best_avg_imi = fallback_eval['avg_imi']
            best_feasible = bool(fallback_eval['feasible'])
            best_ratios = fallback_eval['ratios'].copy()
            best_imis = fallback_eval['imis'].copy()

        return SingleRunResult(
            best_theta=best_theta,
            best_loss=best_loss,
            best_penalty=best_penalty,
            best_avg_imi=best_avg_imi,
            best_feasible=best_feasible,
            best_ratios=best_ratios,
            best_imis=best_imis,
            trace_means=np.array(trace_means),
            trace_stds=np.array(trace_stds),
            trace_best_losses=np.array(trace_best_losses),
            trace_best_penalties=np.array(trace_best_penalties),
            trace_best_avg_imi=np.array(trace_best_avg_imi),
            trace_best_imis=np.array(trace_best_imis),
            trace_feasible_fraction=np.array(trace_feasible_fraction),
        )

# ============================================================
# 6) MULTI-RUN ANALYSIS
# ============================================================
def run_multiple_optimizations(n_runs: int = N_RUNS, n_samples: int = N_SAMPLES_FEASIBILITY) -> Dict[str, Any]:
    print("=" * 80)
    print(f"MULTI-RUN IFFL OPTIMIZATION ({n_runs} runs)")
    print("Objective: IMI < 0, preserve ratio ordering, then maximize average |IMI|")
    print("=" * 80)

    print("\n[Phase 1] Feasibility pre-scan...")
    rng = np.random.default_rng(RANDOM_SEED)
    phi_samples = LO_PHI + (HI_PHI - LO_PHI) * rng.random((n_samples, 9))
    theta_samples = np.exp(phi_samples)

    best_scan_loss = float('inf')
    best_scan_theta = None
    best_scan_avg_imi = -float('inf')

    for theta in theta_samples:
        ev = evaluate_theta(theta)
        if (ev['loss'] < best_scan_loss - 1e-15) or (
            abs(ev['loss'] - best_scan_loss) <= 1e-15 and ev['avg_imi'] > best_scan_avg_imi
        ):
            best_scan_loss = ev['loss']
            best_scan_avg_imi = ev['avg_imi']
            best_scan_theta = theta.copy()

    print(f"  Best pre-scan loss: {best_scan_loss:.6g}")
    print(f"  Best pre-scan average |IMI|: {best_scan_avg_imi:.6g}")

    print(f"\n[Phase 2] Running {n_runs} CEM optimizations...")
    all_results: List[SingleRunResult] = []
    all_trace_means = []
    all_trace_best_losses = []
    all_trace_best_penalties = []
    all_trace_best_avg_imi = []
    all_trace_best_imis = []
    all_trace_feasible_fraction = []

    for run_id in range(n_runs):
        seed = 100 + run_id * 10
        config = CEMConfig(
            n_iter=N_ITER_CEM,
            population_size=POPULATION_SIZE,
            elite_frac=ELITE_FRAC,
            alpha=ALPHA,
            cov_floor=COV_FLOOR,
            seed=seed,
            verbose=False,
        )
        optimizer = CEMOptimizer(config)
        result = optimizer.optimize(init_center=best_scan_theta)
        all_results.append(result)
        all_trace_means.append(result.trace_means)
        all_trace_best_losses.append(result.trace_best_losses)
        all_trace_best_penalties.append(result.trace_best_penalties)
        all_trace_best_avg_imi.append(result.trace_best_avg_imi)
        all_trace_best_imis.append(result.trace_best_imis)
        all_trace_feasible_fraction.append(result.trace_feasible_fraction)

        print(
            f"  Run {run_id + 1:3d}/{n_runs}: "
            f"feasible = {result.best_feasible}, "
            f"best penalty = {result.best_penalty:.6g}, "
            f"best avg |IMI| = {result.best_avg_imi:.6g}"
        )

    successful_results = [r for r in all_results if r.best_feasible]
    if len(successful_results) == 0:
        print("\nWARNING: No feasible run found. Summary will use all runs.")
        successful_results = all_results[:]

    all_thetas = np.array([r.best_theta for r in successful_results])
    trace_means_array = np.array(all_trace_means)
    trace_best_losses_array = np.array(all_trace_best_losses)
    trace_best_penalties_array = np.array(all_trace_best_penalties)
    trace_best_avg_imi_array = np.array(all_trace_best_avg_imi)
    trace_best_imis_array = np.array(all_trace_best_imis)
    trace_feasible_fraction_array = np.array(all_trace_feasible_fraction)

    mean_across_runs = np.mean(trace_means_array, axis=0)
    std_across_runs = np.std(trace_means_array, axis=0, ddof=1)
    best_loss_mean = np.mean(trace_best_losses_array, axis=0)
    best_loss_std = np.std(trace_best_losses_array, axis=0, ddof=1)
    best_penalty_mean = np.nanmean(trace_best_penalties_array, axis=0)
    best_penalty_std = np.nanstd(trace_best_penalties_array, axis=0, ddof=1)
    best_avg_imi_mean = np.nanmean(trace_best_avg_imi_array, axis=0)
    best_avg_imi_std = np.nanstd(trace_best_avg_imi_array, axis=0, ddof=1)
    best_imi_mean_by_motif = np.nanmean(trace_best_imis_array, axis=0)  # (iter, motif)
    best_imi_std_by_motif = np.nanstd(trace_best_imis_array, axis=0, ddof=1)
    feasible_fraction_mean = np.mean(trace_feasible_fraction_array, axis=0)
    feasible_fraction_std = np.std(trace_feasible_fraction_array, axis=0, ddof=1)

    theta_means = np.mean(all_thetas, axis=0)
    theta_stds = np.std(all_thetas, axis=0, ddof=1)
    theta_cv = np.where(theta_means != 0, theta_stds / theta_means, np.nan)

    best_avg_imi_by_run = np.array([r.best_avg_imi for r in all_results], float)
    best_penalty_by_run = np.array([r.best_penalty for r in all_results], float)
    feasible_by_run = np.array([int(r.best_feasible) for r in all_results], int)

    all_imis_by_run = {
        motif: np.array([r.best_imis.get(motif, np.nan) for r in all_results], float)
        for motif in INCOHERENT_MOTIFS
    }
    all_ratios_by_run = {
        f"|{a}|/|{b}|": np.array([r.best_ratios.get(f"|{a}|/|{b}|", np.nan) for r in all_results], float)
        for (a, b) in IFFL_TARGETS
    }

    imi_means = {m: float(np.nanmean([r.best_imis.get(m, np.nan) for r in successful_results])) for m in INCOHERENT_MOTIFS}
    imi_stds = {m: float(np.nanstd([r.best_imis.get(m, np.nan) for r in successful_results], ddof=1)) for m in INCOHERENT_MOTIFS}
    ratio_means = {k: float(np.nanmean([r.best_ratios.get(k, np.nan) for r in successful_results])) for k in all_ratios_by_run}
    ratio_stds = {k: float(np.nanstd([r.best_ratios.get(k, np.nan) for r in successful_results], ddof=1)) for k in all_ratios_by_run}

    sensitivity_results = normalized_sensitivity_analysis(theta_means, PARAM_NAMES)

    return {
        'all_results': all_results,
        'successful_results': successful_results,
        'all_thetas': all_thetas,
        'theta_means': theta_means,
        'theta_stds': theta_stds,
        'theta_cv': theta_cv,
        'imi_means': imi_means,
        'imi_stds': imi_stds,
        'ratio_means': ratio_means,
        'ratio_stds': ratio_stds,
        'best_avg_imi_by_run': best_avg_imi_by_run,
        'best_penalty_by_run': best_penalty_by_run,
        'feasible_by_run': feasible_by_run,
        'all_imis_by_run': all_imis_by_run,
        'all_ratios_by_run': all_ratios_by_run,
        'mean_across_runs': mean_across_runs,
        'std_across_runs': std_across_runs,
        'best_loss_mean': best_loss_mean,
        'best_loss_std': best_loss_std,
        'best_penalty_mean': best_penalty_mean,
        'best_penalty_std': best_penalty_std,
        'best_avg_imi_mean': best_avg_imi_mean,
        'best_avg_imi_std': best_avg_imi_std,
        'best_imi_mean_by_motif': best_imi_mean_by_motif,
        'best_imi_std_by_motif': best_imi_std_by_motif,
        'feasible_fraction_mean': feasible_fraction_mean,
        'feasible_fraction_std': feasible_fraction_std,
        'sensitivity_results': sensitivity_results,
        'n_runs': n_runs,
        'n_success': len(successful_results),
        'n_iter': N_ITER_CEM,
        'best_scan_theta': best_scan_theta,
    }

# ============================================================
# 7) EXPLORATION ANALYSIS
# ============================================================
def compute_pairwise_distances(multi_results: Dict[str, Any]) -> np.ndarray:
    """Pairwise distances among final feasible best solutions (run-by-run)."""
    all_thetas = np.array(multi_results['all_thetas'])
    n_runs = len(all_thetas)

    if n_runs == 0:
        return np.zeros((0, 0))

    log_thetas = np.log(all_thetas)
    log_lo = np.log(LO_THETA)
    log_hi = np.log(HI_THETA)
    log_range = log_hi - log_lo

    normalized = (log_thetas - log_lo) / log_range
    normalized = np.clip(normalized, 0, 1)

    dist_matrix = np.zeros((n_runs, n_runs))
    for i in range(n_runs):
        for j in range(n_runs):
            raw_dist = np.sqrt(np.mean((normalized[i] - normalized[j]) ** 2))
            dist_matrix[i, j] = raw_dist

    print("\n" + "=" * 80)
    print("EXPLORATION ANALYSIS (Run-by-run distances among final feasible solutions)")
    print("=" * 80)
    print("Distance = 0: identical solutions")
    print("Distance = 1: solutions at opposite ends of parameter range")
    print(f"Number of successful runs: {n_runs}")

    if n_runs > 1:
        off_diag = np.array([dist_matrix[i, j] for i in range(n_runs) for j in range(i + 1, n_runs)])
        print(f"Mean pairwise distance: {np.mean(off_diag):.4f}")
        print(f"Std of pairwise distances: {np.std(off_diag, ddof=1):.4f}")
        print(f"Min pairwise distance: {np.min(off_diag):.4f}")
        print(f"Max pairwise distance: {np.max(off_diag):.4f}")
    else:
        print("Only one successful run, so pairwise summary is not defined.")

    return dist_matrix


def compute_parameter_correlation_matrix(multi_results: Dict[str, Any]) -> np.ndarray:
    all_thetas = np.array(multi_results['all_thetas'])
    if len(all_thetas) < 2:
        return np.eye(len(PARAM_NAMES))

    corr_matrix = np.corrcoef(all_thetas.T)

    print("\nHigh parameter correlations (|r| > 0.8):")
    high_corr_found = False
    for i in range(len(PARAM_NAMES)):
        for j in range(i + 1, len(PARAM_NAMES)):
            if abs(corr_matrix[i, j]) > 0.8:
                print(f"  {PARAM_NAMES[i]} ↔ {PARAM_NAMES[j]}: r = {corr_matrix[i, j]:.3f}")
                high_corr_found = True
    if not high_corr_found:
        print("  No strong correlations found (|r| <= 0.8)")

    return corr_matrix

# ============================================================
# 8) DATA EXPORT
# ============================================================
def export_all_data(multi_results: Dict[str, Any]):
    print("\n" + "=" * 80)
    print("EXPORTING DATA TO .DAT FILES")
    print("=" * 80)

    iterations = np.arange(1, multi_results['n_iter'] + 1)

    # File 1: diagnostics figure data
    log_means = multi_results['mean_across_runs']
    log_stds = multi_results['std_across_runs']
    orig_means = np.exp(log_means)
    orig_stds = orig_means * log_stds
    header1 = [
        "Iteration",
        "Best_Avg_IMI_Mean", "Best_Avg_IMI_SD",
        "Feasible_Fraction_Mean", "Feasible_Fraction_SD",
    ]
    for p in PARAM_NAMES:
        header1.extend([p, f"{p}_SD"])

    data1 = []
    for i in range(multi_results['n_iter']):
        row = [
            iterations[i],
            multi_results['best_avg_imi_mean'][i],
            multi_results['best_avg_imi_std'][i],
            multi_results['feasible_fraction_mean'][i],
            multi_results['feasible_fraction_std'][i],
        ]
        for j in range(len(PARAM_NAMES)):
            row.extend([orig_means[i, j], orig_stds[i, j]])
        data1.append(row)
    export_to_dat(
        "incoherent-yeast-OR-diagnostics.dat",
        np.array(data1, dtype=object),
        header1,
        "Diagnostics figure data: best average |IMI|, feasible fraction, and parameter-center evolution",
    )

    # File 2: best penalty convergence plot data
    data2 = [[iterations[i], multi_results['best_penalty_mean'][i], multi_results['best_penalty_std'][i]]
             for i in range(multi_results['n_iter'])]
    export_to_dat(
        "incoherent-yeast-OR-best_penalty_convergence.dat",
        np.array(data2, dtype=object),
        ["Iteration", "Best_Penalty_Mean", "Best_Penalty_SD"],
        "Best penalty versus iteration (mean ± SD across runs)",
    )

    # File 3: IMI trajectories by motif
    header3 = ["Iteration"]
    for motif in MOTIF_NAMES:
        header3.extend([f"{motif}_Mean", f"{motif}_SD"])
    data3 = []
    for i in range(multi_results['n_iter']):
        row = [iterations[i]]
        for j, motif in enumerate(MOTIF_NAMES):
            row.extend([
                multi_results['best_imi_mean_by_motif'][i, j],
                multi_results['best_imi_std_by_motif'][i, j],
            ])
        data3.append(row)
    export_to_dat(
        "incoherent-yeast-OR-imi_trajectory.dat",
        np.array(data3, dtype=object),
        header3,
        "Best-so-far IMI values for each motif versus iteration (mean ± SD across runs)",
    )

    # File 4: final IMI bar plot data
    data4 = [[motif, multi_results['imi_means'][motif], multi_results['imi_stds'][motif]] for motif in MOTIF_NAMES]
    export_to_dat(
        "incoherent-yeast-OR-final_imi_summary.dat",
        np.array(data4, dtype=object),
        ["Motif", "IMI_Mean", "IMI_SD"],
        "Final IMI values for each motif (mean ± SD across successful runs)",
    )

    # File 5: ratio satisfaction scatter-ready data
    ratio_keys = [f"|{a}|/|{b}|" for (a, b) in IFFL_TARGETS]
    data5 = []
    for i, key in enumerate(ratio_keys):
        y = multi_results['all_ratios_by_run'][key]
        valid = np.isfinite(y)
        rng = np.random.default_rng(123 + i)
        jitter = rng.normal(0, 0.04, size=np.sum(valid))
        a, b = key.replace('|', '').split('/')
        lo, hi = IFFL_TARGETS[(a, b)]
        mean_val = np.nanmean(y)
        std_val = np.nanstd(y, ddof=1)

        valid_indices = np.where(valid)[0]
        for k, run_idx in enumerate(valid_indices):
            x_base = float(i)
            x_jittered = float(i + jitter[k])
            data5.append([
                key,
                run_idx + 1,
                x_base,
                x_jittered,
                y[run_idx],
                lo,
                hi,
                mean_val,
                std_val,
            ])

    export_to_dat(
        "incoherent-yeast-OR-ratio_satisfaction_scatter.dat",
        np.array(data5, dtype=object),
        ["Ratio", "Run", "X_Base", "X_Jittered", "Ratio_Value", "Target_Low", "Target_High", "Mean", "SD"],
        "Scatter-ready ratio data with target intervals and summary statistics",
    )

    # File 6: average |IMI| distribution plot data
    data6 = [[i + 1, multi_results['best_avg_imi_by_run'][i]] for i in range(multi_results['n_runs'])]
    export_to_dat(
        "incoherent-yeast-OR-avg_imi_by_run.dat",
        np.array(data6, dtype=object),
        ["Run", "Best_Avg_IMI"],
        "Final best average |IMI| in each run",
    )

    # File 7: pairwise distances heatmap data
    dist_matrix = compute_pairwise_distances(multi_results)
    if dist_matrix.size > 0:
        n_success = dist_matrix.shape[0]
        header7 = [f"Run_{i + 1}" for i in range(n_success)]
        export_to_dat(
            "incoherent-yeast-OR-pairwise_distances.dat",
            dist_matrix,
            header7,
            "Normalized pairwise distances among final feasible best solutions",
        )

    # File 8: parameter correlations heatmap data
    corr_matrix = compute_parameter_correlation_matrix(multi_results)
    export_to_dat(
        "incoherent-yeast-OR-parameter_correlations.dat",
        corr_matrix,
        PARAM_NAMES,
        "Parameter correlation matrix across final feasible solutions",
    )

    # File 9: normalized sensitivity bar plot data
    sens = multi_results['sensitivity_results']
    header9 = ["Parameter"] + MOTIF_NAMES + ["Avg_Across_Motifs"]
    data9 = []
    for p in PARAM_NAMES:
        row = [p] + [sens['motif_sensitivity'][m][p] for m in MOTIF_NAMES] + [sens['normalized_sensitivity'][p]]
        data9.append(row)
    export_to_dat(
        "incoherent-yeast-OR-normalized_sensitivity.dat",
        np.array(data9, dtype=object),
        header9,
        "Average normalized sensitivity of each parameter for each motif",
    )

    # File 10: parameter CV bar plot data
    data10 = [[PARAM_NAMES[i], multi_results['theta_cv'][i]] for i in range(len(PARAM_NAMES))]
    export_to_dat(
        "incoherent-yeast-OR-parameter_cv.dat",
        np.array(data10, dtype=object),
        ["Parameter", "CV"],
        "Parameter coefficient of variation across successful runs",
    )

    # File 11: violin plot data
    data11 = []
    for i, theta in enumerate(multi_results['all_thetas']):
        data11.append([i + 1] + list(theta))
    export_to_dat(
        "incoherent-yeast-OR-violin_data.dat",
        np.array(data11, dtype=object),
        ["Successful_Run"] + PARAM_NAMES,
        "Raw final feasible parameter values",
    )

    # File 12: final parameter summary table data
    data12 = [[PARAM_NAMES[i], multi_results['theta_means'][i], multi_results['theta_stds'][i]]
              for i in range(len(PARAM_NAMES))]
    export_to_dat(
        "incoherent-yeast-OR-final_parameter_summary.dat",
        np.array(data12, dtype=object),
        ["Parameter", "Mean", "Std_Dev"],
        "Mean and standard deviation of final feasible parameter values",
    )

    print("\n" + "=" * 80)
    print(f"All data exported to '{EXPORT_DIR}/' directory")
    print("=" * 80)

# ============================================================
# 9) VISUALIZATION
# ============================================================
def plot_diagnostics_figure(multi_results: Dict[str, Any]):
    """
    Figure: best average |IMI| convergence + feasible fraction + parameter-center evolution
            + pairwise distances + parameter correlations
    """
    n_iter = multi_results['n_iter']
    iterations = np.arange(1, n_iter + 1)

    fig = plt.figure(figsize=FIG_SIZE_DIAGNOSTICS)
    gs = fig.add_gridspec(6, 3, height_ratios=[1, 1, 1.4, 1.4, 1, 1], hspace=0.6, wspace=0.4)

    ax_avg = fig.add_subplot(gs[0, :2])
    mean_avg = multi_results['best_avg_imi_mean']
    std_avg = multi_results['best_avg_imi_std']
    ax_avg.plot(iterations, mean_avg, linewidth=2, label='Mean across runs')
    ax_avg.fill_between(iterations, mean_avg - std_avg, mean_avg + std_avg, alpha=0.3, label='±1 SD')
    ax_avg.set_xlabel('Iteration')
    ax_avg.set_ylabel('Best average |IMI|')
    ax_avg.set_title('Objective Convergence: Best Average |IMI|', fontweight='bold')
    ax_avg.grid(True, alpha=0.3)
    ax_avg.legend(loc='best')

    ax_feas = fig.add_subplot(gs[0, 2])
    feas_mean = multi_results['feasible_fraction_mean']
    feas_std = multi_results['feasible_fraction_std']
    ax_feas.plot(iterations, feas_mean, linewidth=2, label='Mean across runs')
    ax_feas.fill_between(iterations, np.maximum(0, feas_mean - feas_std), np.minimum(1, feas_mean + feas_std), alpha=0.3)
    ax_feas.set_xlabel('Iteration')
    ax_feas.set_ylabel('Feasible fraction')
    ax_feas.set_title('Population Feasibility', fontweight='bold')
    ax_feas.set_ylim(0, 1.05)
    ax_feas.grid(True, alpha=0.3)

    mean_across_runs = multi_results['mean_across_runs']
    std_across_runs = multi_results['std_across_runs']
    for i, name in enumerate(PARAM_NAMES):
        row = 1 + (i // 3)
        col = i % 3
        ax = fig.add_subplot(gs[row, col])
        param_mean = mean_across_runs[:, i]
        param_std = std_across_runs[:, i]
        ax.plot(iterations, param_mean, linewidth=2)
        ax.fill_between(iterations, param_mean - param_std, param_mean + param_std, alpha=0.3)
        ax.set_xlabel('Iteration', fontsize=9)
        ax.set_ylabel(f'log({name})', fontsize=9)
        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)

    ax_dist = fig.add_subplot(gs[4:6, 0])
    dist_matrix = compute_pairwise_distances(multi_results)
    if dist_matrix.size > 0:
        im_dist = ax_dist.imshow(dist_matrix, cmap='hot', interpolation='nearest')
        ax_dist.set_xlabel('Successful run index')
        ax_dist.set_ylabel('Successful run index')
        ax_dist.set_title('Pairwise Distances Among Final Feasible Solutions', fontsize=11, fontweight='bold')
        plt.colorbar(im_dist, ax=ax_dist)
    else:
        ax_dist.text(0.5, 0.5, 'No feasible run', ha='center', va='center')
        ax_dist.set_axis_off()

    ax_corr = fig.add_subplot(gs[4:6, 1:3])
    corr_matrix = compute_parameter_correlation_matrix(multi_results)
    im_corr = ax_corr.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, interpolation='nearest')
    ax_corr.set_xticks(np.arange(len(PARAM_NAMES)))
    ax_corr.set_yticks(np.arange(len(PARAM_NAMES)))
    ax_corr.set_xticklabels(PARAM_NAMES, rotation=45, ha='right', fontsize=9)
    ax_corr.set_yticklabels(PARAM_NAMES, fontsize=9)
    ax_corr.set_title('Parameter Correlations Across Final Feasible Solutions', fontsize=11, fontweight='bold')
    plt.colorbar(im_corr, ax=ax_corr, label='Correlation')

    plt.tight_layout()
    plt.show()
    return fig


def plot_best_penalty_convergence(multi_results: Dict[str, Any]):
    fig, ax = plt.subplots(figsize=FIG_SIZE_STANDARD)
    iterations = np.arange(1, multi_results['n_iter'] + 1)
    mean_penalty = multi_results['best_penalty_mean']
    std_penalty = multi_results['best_penalty_std']
    ax.plot(iterations, mean_penalty, linewidth=2, label='Mean across runs')
    ax.fill_between(iterations, np.maximum(0, mean_penalty - std_penalty), mean_penalty + std_penalty, alpha=0.3, label='±1 SD')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Best penalty so far')
    ax.set_title('Best Penalty Convergence', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    plt.tight_layout()
    plt.show()
    return fig


def plot_imi_trajectory(multi_results: Dict[str, Any]):
    fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE_IMI_TRAJECTORY)
    iterations = np.arange(1, multi_results['n_iter'] + 1)
    axes = axes.ravel()
    for j, motif in enumerate(MOTIF_NAMES):
        ax = axes[j]
        mean_vals = multi_results['best_imi_mean_by_motif'][:, j]
        std_vals = multi_results['best_imi_std_by_motif'][:, j]
        ax.plot(iterations, mean_vals, linewidth=2, label='Mean across runs')
        ax.fill_between(iterations, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3, label='±1 SD')
        ax.set_xlabel('Iteration')
        ax.set_ylabel(f'{motif} IMI')
        ax.set_title(f'{motif}: IMI trajectory', fontweight='bold')
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.legend(loc='best')
    plt.tight_layout()
    plt.show()
    return fig


def plot_final_imi_bar(multi_results: Dict[str, Any]):
    fig, ax = plt.subplots(figsize=FIG_SIZE_STANDARD)
    means = [multi_results['imi_means'][m] for m in MOTIF_NAMES]
    stds = [multi_results['imi_stds'][m] for m in MOTIF_NAMES]
    bars = ax.bar(MOTIF_NAMES, means, yerr=stds, capsize=5, edgecolor='black', alpha=0.75)
    ax.set_xlabel('Motif')
    ax.set_ylabel('Final IMI')
    ax.set_title('Final IMI of Each Incoherent FFL Type', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    for bar, mean in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{mean:.4f}', ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.show()
    return fig


def plot_ratio_satisfaction(multi_results: Dict[str, Any]):
    """Scatter of final ratio values across runs, with target bands."""
    ratio_keys = [f"|{a}|/|{b}|" for (a, b) in IFFL_TARGETS]
    fig, ax = plt.subplots(figsize=FIG_SIZE_STANDARD)

    x_positions = np.arange(len(ratio_keys))

    for i, key in enumerate(ratio_keys):
        y = multi_results['all_ratios_by_run'][key]
        valid = np.isfinite(y)
        jitter = np.random.default_rng(123 + i).normal(0, 0.04, size=np.sum(valid))
        ax.scatter(np.full(np.sum(valid), x_positions[i]) + jitter, y[valid], alpha=0.4, s=20)

        a, b = key.replace('|', '').split('/')
        lo, hi = IFFL_TARGETS[(a, b)]
        ax.fill_between([x_positions[i] - 0.25, x_positions[i] + 0.25], lo, hi, alpha=0.25)
        ax.hlines([lo, hi], x_positions[i] - 0.25, x_positions[i] + 0.25, linewidth=2)

        mean_val = np.nanmean(y)
        std_val = np.nanstd(y, ddof=1)
        ax.errorbar(x_positions[i], mean_val, yerr=std_val, fmt='o', capsize=5, linewidth=2)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(ratio_keys)
    ax.set_ylabel('Final ratio value')
    ax.set_title('Ratio Satisfaction Across Runs\n(Shaded bands = target intervals)', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.show()
    return fig


def plot_avg_imi_distribution(multi_results: Dict[str, Any]):
    fig, ax = plt.subplots(figsize=FIG_SIZE_STANDARD)
    values = multi_results['best_avg_imi_by_run']
    finite = np.isfinite(values)
    ax.hist(values[finite], bins=30, alpha=0.75, edgecolor='black')
    ax.set_xlabel('Final best average |IMI|')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Final Best Average |IMI| Across Runs', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.show()
    return fig


def plot_normalized_sensitivity_bar(multi_results: Dict[str, Any]):
    sens = multi_results['sensitivity_results']
    fig, ax = plt.subplots(figsize=FIG_SIZE_SENSITIVITY)
    x = np.arange(len(PARAM_NAMES))
    width = 0.18
    for j, motif in enumerate(MOTIF_NAMES):
        vals = [sens['motif_sensitivity'][motif][p] for p in PARAM_NAMES]
        ax.bar(x + (j - 1.5) * width, vals, width=width, edgecolor='black', alpha=0.75, label=motif)
    ax.set_xticks(x)
    ax.set_xticklabels(PARAM_NAMES, rotation=45, ha='right')
    ax.set_ylabel('Average normalized sensitivity')
    ax.set_title('Average Normalized Sensitivity by Motif', fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.show()
    return fig


def plot_parameter_cv_bar(multi_results: Dict[str, Any]):
    fig, ax = plt.subplots(figsize=FIG_SIZE_CV)
    vals = multi_results['theta_cv']
    bars = ax.bar(PARAM_NAMES, vals, edgecolor='black', alpha=0.75)
    ax.set_ylabel('CV across successful runs')
    ax.set_title('Parameter CV Across Runs', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(axis='x', rotation=45)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{val:.3f}', ha='center', va='bottom', fontsize=8)
    plt.tight_layout()
    plt.show()
    return fig


def plot_parameter_violin(multi_results: Dict[str, Any]):
    all_thetas = np.array(multi_results['all_thetas'])
    if len(all_thetas) == 0:
        print("No feasible solutions available for violin plot.")
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_SIZE_VIOLIN)

    beta_indices = [0, 1, 2]
    other_indices = [3, 4, 5, 6, 7, 8]
    beta_names = [PARAM_NAMES[i] for i in beta_indices]
    other_names = [PARAM_NAMES[i] for i in other_indices]

    beta_data = [all_thetas[:, i] for i in beta_indices]
    other_data = [all_thetas[:, i] for i in other_indices]

    pos1 = np.arange(len(beta_indices))
    parts1 = ax1.violinplot(beta_data, positions=pos1, showmeans=True, showmedians=True, widths=0.7)
    for pc in parts1['bodies']:
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
    parts1['cmeans'].set_linewidth(2)
    parts1['cmedians'].set_linewidth(2)
    ax1.set_xticks(pos1)
    ax1.set_xticklabels(beta_names, rotation=45, ha='right')
    ax1.set_ylabel('Parameter value')
    ax1.set_title('Degradation-rate parameters', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    pos2 = np.arange(len(other_indices))
    parts2 = ax2.violinplot(other_data, positions=pos2, showmeans=True, showmedians=True, widths=0.7)
    for pc in parts2['bodies']:
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
    parts2['cmeans'].set_linewidth(2)
    parts2['cmedians'].set_linewidth(2)
    ax2.set_yscale('log')
    ax2.set_xticks(pos2)
    ax2.set_xticklabels(other_names, rotation=45, ha='right')
    ax2.set_ylabel('Parameter value (log scale)')
    ax2.set_title('Threshold and mean-expression parameters', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.suptitle('Distributions of Final Feasible Parameter Values', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
    return fig

# ============================================================
# 10) MAIN EXECUTION
# ============================================================
def run_iffl_multi_optimization(n_runs: int = N_RUNS, n_samples: int = N_SAMPLES_FEASIBILITY):
    multi_results = run_multiple_optimizations(n_runs=n_runs, n_samples=n_samples)

    export_all_data(multi_results)

    print("\n[Phase] Plotting diagnostics figure...")
    plot_diagnostics_figure(multi_results)

    print("\n[Phase] Plotting best penalty convergence...")
    plot_best_penalty_convergence(multi_results)

    print("\n[Phase] Plotting IMI trajectory for each motif...")
    plot_imi_trajectory(multi_results)

    print("\n[Phase] Plotting final IMI bar plot...")
    plot_final_imi_bar(multi_results)

    print("\n[Phase] Plotting ratio satisfaction...")
    plot_ratio_satisfaction(multi_results)

    print("\n[Phase] Plotting average |IMI| distribution...")
    plot_avg_imi_distribution(multi_results)

    print("\n[Phase] Plotting normalized sensitivity bar plot...")
    plot_normalized_sensitivity_bar(multi_results)

    print("\n[Phase] Plotting parameter CV bar plot...")
    plot_parameter_cv_bar(multi_results)

    print("\n[Phase] Plotting parameter violin plots...")
    plot_parameter_violin(multi_results)

    return multi_results


if __name__ == "__main__":
    results = run_iffl_multi_optimization()