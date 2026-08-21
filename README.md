# An information-theoretic perspective on feed-forward loop abundances in transcriptional networks

This repository contains the symbolic derivations and numerical codes associated with the manuscript **“An information-theoretic framework for understanding feed-forward loop abundances in transcriptional networks.”** The codes implement the linear-noise calculations, mutual-information decomposition, cross-entropy optimization, open-loop reference calculations, and Gillespie simulations used in the study.

## Authors and affiliations

1. **Mintu Nandi**  
   Universal Biology Institute, The University of Tokyo, 7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033, Japan  
   Email: `mintunandi@ubi.s.u-tokyo.ac.jp`
2. **Sudip Chattopadhyay**  
   Department of Chemistry, Indian Institute of Engineering Science and Technology, Shibpur, Howrah 711103, India  
   Email: `sudip@chem.iiests.ac.in`
3. **Suman K Banik**  
   Department of Chemical Sciences, Bose Institute, EN 80, Sector V, Bidhan Nagar, Kolkata 700091, India  
   Email: `skbanik@jcbose.ac.in`

## Repository contents

The supplied repository contains **95 code files**:

- **5 Wolfram Mathematica notebooks** for symbolic and analytical calculations;
- **24 Python scripts** for cross-entropy-method (CEM) optimization;
- **64 Python scripts** for Gillespie simulations of FFL and open-loop circuits; and
- **2 Python scripts** that combine the stochastic FFL and open-loop results to obtain IMI.

Code filenames identify their metric, organism, circuit, motif, and gate where applicable. The directory hierarchy and generated-data filenames follow the workflow described below.

Each of the 95 code files begins with a file-specific documentation block describing its purpose, inputs, outputs, execution path, model condition, and correspondence with the notation used in the manuscript and Supplementary Material.

## Scientific workflow

1. **Symbolic moment calculation**  
   `Theory/Moments-calculations-FFL-OpenLoop.nb` solves the steady-state Lyapunov equations for the FFL and its open-loop reference.

2. **Abundance-guided CEM optimization**  
   The scripts in `Numerics/CEM-optimization-IMI/` optimize IMI. The two control directories repeat the same search with pathway MI or total MI.

3. **Stochastic validation**  
   The scripts in `Numerics/Gillespie-after-optimization-IMI/` simulate each optimized FFL and its open-loop reference. The organism-level `Calculations-*.py` scripts calculate

   $$
   I_{\mathrm{int}}(X;Z)=I_{\mathrm{FFL}}(X;Z)-I_{\mathrm{OL}}(X,\widetilde X;Z).
   $$

4. **Analytical FFL calculations**  
   `Theory/FFL/AND.nb` and `Theory/FFL/OR.nb` evaluate the analytical fluctuation terms, total MI, pathway MI, IMI, pathway-interference strength, and local pathway sensitivities for AND- and OR-gated FFLs using the optimized parameter set given in Supplementary Material.

5. **Analytical open-loop calculations**  
   `Theory/OL/Open-loop-AND.nb` and `Theory/OL/Open-loop-OR.nb` calculate the joint information between the two independent open-loop inputs and the output using the optimized parameter set for FFL, see Supplementary Material.

## Directory structure

```text
.
├── README.md
├── Theory/
│   ├── Moments-calculations-FFL-OpenLoop.nb
│   ├── FFL/
│   │   ├── AND.nb
│   │   └── OR.nb
│   └── OL/
│       ├── Open-loop-AND.nb
│       └── Open-loop-OR.nb
└── Numerics/
    ├── CEM-optimization-IMI/                 # 8 scripts
    ├── CEM-optimization-path-MI/             # 8 control scripts
    ├── CEM-optimization-total-MI/            # 8 control scripts
    └── Gillespie-after-optimization-IMI/
        ├── Ecoli/
        │   ├── Calculations-ecoli.py
        │   ├── FFL/AND/                      # 8 motif scripts
        │   ├── FFL/OR/                       # 8 motif scripts
        │   ├── OL/AND/                       # 8 open-loop scripts
        │   └── OL/OR/                        # 8 open-loop scripts
        └── Yeast/
            ├── Calculations-yeast.py
            ├── FFL/AND/                      # 8 motif scripts
            ├── FFL/OR/                       # 8 motif scripts
            ├── OL/AND/                       # 8 open-loop scripts
            └── OL/OR/                        # 8 open-loop scripts
```

All paths are case-sensitive on Linux and macOS. In particular, the repository uses `Ecoli`, `Yeast`, `FFL`, `OL`, `AND`, and `OR` exactly as shown above.

## Software requirements

### Python codes

- Python **3.9 or newer** is recommended.
- NumPy
- Matplotlib, required by the CEM scripts only

The exact Python version used to generate the supplied results is not recorded in the archive. The source files use standard-library modules plus NumPy and Matplotlib.

Install the Python dependencies with

```bash
python -m pip install numpy matplotlib
```

### Mathematica notebooks

- `Theory/FFL/AND.nb` and `Theory/FFL/OR.nb` were created with Wolfram 14.3.
- `Theory/OL/Open-loop-AND.nb` and `Theory/OL/Open-loop-OR.nb` were created with Wolfram 14.3.
- `Theory/Moments-calculations-FFL-OpenLoop.nb` was created with Mathematica 14.3.

Wolfram Mathematica 14.3 or a compatible later release is recommended for opening all notebooks.

## Running the Mathematica notebooks

Open a notebook in Mathematica and evaluate its input cells. The notebooks use hard-coded analytical expressions and parameter sets; they do not read external files.

| Notebook | Purpose |
|---|---|
| `Theory/Moments-calculations-FFL-OpenLoop.nb` | Symbolically solves the FFL and open-loop Lyapunov equations for variances and covariances. The notebook contains `Quit[]`; evaluate the FFL and open-loop calculation groups separately if the kernel closes. |
| `Theory/FFL/AND.nb` | Analytical FFL calculations for all eight motifs, both organisms, and an AND output gate. |
| `Theory/FFL/OR.nb` | Analytical FFL calculations for all eight motifs, both organisms, and an OR output gate. |
| `Theory/OL/Open-loop-AND.nb` | Analytical open-loop reference calculations for the AND gate. |
| `Theory/OL/Open-loop-OR.nb` | Analytical open-loop reference calculations for the OR gate. |

The FFL notebooks use the optimized parameter values directly and generate evaluated notebook outputs and plots. The open-loop notebooks use the corresponding matched parameter sets.

## Running the CEM optimization codes

### Important working-directory rule

Each CEM script creates a relative output directory such as `coherent_ecoli_AND_data/`. The same output-directory names are reused in the IMI, pathway-MI, and total-MI folders. To keep the three analyses separate, **run each script from its own containing directory**.

Example:

```bash
cd Numerics/CEM-optimization-IMI
python imi-coherent-ecoli-and.py
```

For the pathway-MI control:

```bash
cd Numerics/CEM-optimization-path-MI
python path-mi-coherent-ecoli-and.py
```

For the total-MI control:

```bash
cd Numerics/CEM-optimization-total-MI
python total-mi-coherent-ecoli-and.py
```

The scripts have no command-line interface. Search bounds, target windows, CEM settings, random seeds, sensitivity settings, and plotting sizes are defined in the global-parameter section near the beginning of each file.

### Common CEM settings in the supplied scripts

| Setting | Value |
|---|---:|
| Independent runs | 100 |
| Feasibility pre-scan samples | 1,000 |
| CEM iterations per run | 1,000 |
| Population size per iteration | 1,000 |
| Elite fraction | 0.1 |
| Smoothing parameter (`ALPHA`) | 0.3 |
| Covariance floor | $10^{-6}$ |
| Feasibility tolerance | $10^{-12}$ |
| Base random seed | 42 |
| Sensitivity perturbations | 200 over a relative range of $\pm 0.1$ |

### Statistical conventions used by the CEM scripts

- Across-run variability is reported using the sample standard deviation with denominator $N_{\mathrm{run}}-1$.
- Run–run distances are calculated as the root-mean-square distance in normalized log-parameter space, $D_{rr'}=[p^{-1}\sum_i(q_{r,i}-q_{r',i})^2]^{1/2}$, with $0\le D_{rr'}\le1$.

The default budget is computationally intensive: one script evaluates approximately $100\times1000\times1000=10^8$ candidate parameter vectors, in addition to the pre-scan and post-optimization analyses.

### Search bounds

The optimized vector is

$$
\boldsymbol\theta=(\beta_X,\beta_Y,\beta_Z,K_{XY},K_{XZ},K_{YZ},\langle x\rangle,\langle y\rangle,\langle z\rangle).
$$

| Organism | Lower bounds | Upper bounds |
|---|---|---|
| *E. coli* | `(0.005, 0.005, 0.005, 5, 5, 5, 10, 10, 10)` | `(0.1, 0.1, 0.5, 100, 100, 100, 100, 100, 100)` |
| Yeast | `(0.001, 0.001, 0.001, 20, 20, 20, 50, 50, 50)` | `(0.1, 0.1, 0.1, 500, 500, 500, 500, 500, 500)` |

### CEM scripts and target windows

| Script | Optimized quantity | Organism | Class | Gate | Target ratio windows |
|---|---:|---|---|---|---|
| `Numerics/CEM-optimization-IMI/imi-coherent-ecoli-and.py` | $I_{\mathrm{int}}(X;Z)$ | ecoli | coherent | AND | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-IMI/imi-coherent-ecoli-or.py` | $I_{\mathrm{int}}(X;Z)$ | ecoli | coherent | OR | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-IMI/imi-incoherent-ecoli-and.py` | $I_{\mathrm{int}}(X;Z)$ | ecoli | incoherent | AND | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-IMI/imi-incoherent-ecoli-or.py` | $I_{\mathrm{int}}(X;Z)$ | ecoli | incoherent | OR | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-IMI/imi-coherent-yeast-and.py` | $I_{\mathrm{int}}(X;Z)$ | yeast | coherent | AND | `C1/C2` = 3–8; `C1/C3` = 20–150; `C1/C4` = 20–150 |
| `Numerics/CEM-optimization-IMI/imi-coherent-yeast-or.py` | $I_{\mathrm{int}}(X;Z)$ | yeast | coherent | OR | `C1/C2` = 2–9; `C1/C3` = 10–150; `C1/C4` = 10–150 |
| `Numerics/CEM-optimization-IMI/imi-incoherent-yeast-and.py` | $I_{\mathrm{int}}(X;Z)$ | yeast | incoherent | AND | `I1/I2` = 5–7; `I1/I3` = 18–22; `I1/I4` = 50–100 |
| `Numerics/CEM-optimization-IMI/imi-incoherent-yeast-or.py` | $I_{\mathrm{int}}(X;Z)$ | yeast | incoherent | OR | `I1/I2` = 4–8; `I1/I3` = 15–25; `I1/I4` = 20–150 |
| `Numerics/CEM-optimization-path-MI/path-mi-coherent-ecoli-and.py` | $I_{\mathrm{path}}(X;Z)$ | ecoli | coherent | AND | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-path-MI/path-mi-coherent-ecoli-or.py` | $I_{\mathrm{path}}(X;Z)$ | ecoli | coherent | OR | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-path-MI/path-mi-incoherent-ecoli-and.py` | $I_{\mathrm{path}}(X;Z)$ | ecoli | incoherent | AND | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-path-MI/path-mi-incoherent-ecoli-or.py` | $I_{\mathrm{path}}(X;Z)$ | ecoli | incoherent | OR | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-path-MI/path-mi-coherent-yeast-and.py` | $I_{\mathrm{path}}(X;Z)$ | yeast | coherent | AND | `C1/C2` = 3–8; `C1/C3` = 20–150; `C1/C4` = 20–150 |
| `Numerics/CEM-optimization-path-MI/path-mi-coherent-yeast-or.py` | $I_{\mathrm{path}}(X;Z)$ | yeast | coherent | OR | `C1/C2` = 3–8; `C1/C3` = 20–150; `C1/C4` = 20–150 |
| `Numerics/CEM-optimization-path-MI/path-mi-incoherent-yeast-and.py` | $I_{\mathrm{path}}(X;Z)$ | yeast | incoherent | AND | `I1/I2` = 5–7; `I1/I3` = 18–22; `I1/I4` = 50–100 |
| `Numerics/CEM-optimization-path-MI/path-mi-incoherent-yeast-or.py` | $I_{\mathrm{path}}(X;Z)$ | yeast | incoherent | OR | `I1/I2` = 5–7; `I1/I3` = 18–22; `I1/I4` = 50–100 |
| `Numerics/CEM-optimization-total-MI/total-mi-coherent-ecoli-and.py` | $I(X;Z)$ | ecoli | coherent | AND | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-total-MI/total-mi-coherent-ecoli-or.py` | $I(X;Z)$ | ecoli | coherent | OR | `C1/C2` = 8–10; `C1/C3` = 2–3; `C1/C4` = 4–5 |
| `Numerics/CEM-optimization-total-MI/total-mi-incoherent-ecoli-and.py` | $I(X;Z)$ | ecoli | incoherent | AND | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-total-MI/total-mi-incoherent-ecoli-or.py` | $I(X;Z)$ | ecoli | incoherent | OR | `I1/I2` = 8–10; `I1/I3` = 3.5–4.5; `I1/I4` = 4.5–5.5 |
| `Numerics/CEM-optimization-total-MI/total-mi-coherent-yeast-and.py` | $I(X;Z)$ | yeast | coherent | AND | `C1/C2` = 3–8; `C1/C3` = 20–150; `C1/C4` = 20–150 |
| `Numerics/CEM-optimization-total-MI/total-mi-coherent-yeast-or.py` | $I(X;Z)$ | yeast | coherent | OR | `C1/C2` = 3–8; `C1/C3` = 20–150; `C1/C4` = 20–150 |
| `Numerics/CEM-optimization-total-MI/total-mi-incoherent-yeast-and.py` | $I(X;Z)$ | yeast | incoherent | AND | `I1/I2` = 5–7; `I1/I3` = 18–22; `I1/I4` = 50–100 |
| `Numerics/CEM-optimization-total-MI/total-mi-incoherent-yeast-or.py` | $I(X;Z)$ | yeast | incoherent | OR | `I1/I2` = 5–7; `I1/I3` = 18–22; `I1/I4` = 50–100 |

### CEM outputs

Each script writes 12 tab-delimited, OriginPro-compatible `.dat` files to its condition-specific export directory. The prefix is the script condition, for example `coherent-ecoli-AND-`.

| Output pattern | Contents |
|---|---|
| `<condition>-diagnostics.dat` | Best objective trajectory, feasible fraction, and CEM parameter-center trajectories. |
| `<condition>-best_penalty_convergence.dat` | Mean and standard deviation of the best-so-far penalty versus iteration. |
| `<condition>-imi_trajectory.dat` | Motif-wise best-so-far metric trajectories. The name remains `imi_trajectory` in all three metric folders. |
| `<condition>-final_imi_summary.dat` | Final motif-wise metric mean and standard deviation. The filename retains `imi` in the control folders. |
| `<condition>-ratio_satisfaction_scatter.dat` | Run-level metric ratios with target intervals and summary statistics. |
| `<condition>-avg_imi_by_run.dat` | Final best average optimized metric for every independent run. |
| `<condition>-pairwise_distances.dat` | Run–run distances between final feasible parameter vectors. |
| `<condition>-parameter_correlations.dat` | Correlation matrix of final feasible parameter vectors. |
| `<condition>-normalized_sensitivity.dat` | Finite-perturbation normalized sensitivity for each motif and parameter. |
| `<condition>-parameter_cv.dat` | Coefficient of variation of optimized parameters across successful runs. |
| `<condition>-violin_data.dat` | Run-level final feasible parameter vectors. |
| `<condition>-final_parameter_summary.dat` | Mean and standard deviation of final feasible parameters. |

The scripts also display Matplotlib figures interactively. They do not save the figures automatically. On a headless system, the numerical exports can be generated with a noninteractive Matplotlib backend, for example:

```bash
MPLBACKEND=Agg python imi-coherent-ecoli-and.py
```

### Metric interpretation in the control folders

In the pathway-MI and total-MI control directories, internal identifiers such as `imi_analytical_*`, `imis`, `avg_imi`, and output filenames containing `imi` refer to the metric evaluated in that directory:

- `CEM-optimization-IMI/`: the stored quantity is $I_{\mathrm{int}}(X;Z)$;
- `CEM-optimization-path-MI/`: the stored quantity is $I_{\mathrm{path}}(X;Z)$;
- `CEM-optimization-total-MI/`: the stored quantity is $I(X;Z)$.

## Running the Gillespie simulations

The Gillespie scripts are standalone programs with hard-coded optimized parameter sets. They do not accept command-line arguments. Each script writes one two-column `.dat` file containing the motif label and a Gaussian mutual-information estimate calculated from simulated first and second moments.

Because all output filenames are relative to the current working directory, run the 32 scripts for one organism **from that organism's root directory**. This places all FFL and open-loop output files together, which is required by `Calculations-*.py`.

### *E. coli*

```bash
cd Numerics/Gillespie-after-optimization-IMI/Ecoli
for f in FFL/AND/*.py FFL/OR/*.py OL/AND/*.py OL/OR/*.py; do
    python "$f"
done
python Calculations-ecoli.py --dir .
```

### Yeast

```bash
cd Numerics/Gillespie-after-optimization-IMI/Yeast
for f in FFL/AND/*.py FFL/OR/*.py OL/AND/*.py OL/OR/*.py; do
    python "$f"
done
python Calculations-yeast.py --dir .
```

PowerShell equivalent for one organism directory:

```powershell
$files = Get-ChildItem FFL\AND\*.py, FFL\OR\*.py, OL\AND\*.py, OL\OR\*.py
foreach ($file in $files) { python $file.FullName }
python Calculations-ecoli.py --dir .   # use Calculations-yeast.py in the Yeast directory
```

The FFL scripts write files named

```text
total-MI-abund-<motif>ffl-<gate>-num.dat
```

and the open-loop scripts write

```text
total-MI-abund-<motif>ol-<gate>-num.dat
```

where `<motif>` is `c1`–`c4` or `i1`–`i4`, and `<gate>` is `and` or `or`.

The aggregation scripts subtract the open-loop value from the FFL value and write

```text
sMI-coherent-and-<Organism>-num.dat
sMI-incoherent-and-<Organism>-num.dat
sMI-coherent-or-<Organism>-num.dat
sMI-incoherent-or-<Organism>-num.dat
```

Here `sMI` denotes the signed difference $I_{\mathrm{FFL}}-I_{\mathrm{OL}}$, which corresponds to IMI in the manuscript.

### Stochastic-code behavior

- FFL simulations use 6 reaction channels; open-loop simulations use 8 reaction channels.
- Production-rate constants are calculated from the hard-coded operating point and regulatory factors.
- The reported MI is the Gaussian second-moment expression, not a histogram-based or nonparametric MI estimator.
- Every Gillespie script sets a fixed, context-specific NumPy seed. The seed is listed in the header and assigned to `RANDOM_SEED` near the top of the file.
- No separate burn-in interval is implemented.
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-c1-and.py` computes continuous-time residence-time-weighted moments. The remaining stochastic scripts compute moments from the sequence of simulated post-reaction states.
- The simulation time is `1_000_000` for FFL scripts and `5_000_000` for open-loop scripts.

The documentation header at the top of every stochastic script states its model condition, parameter set, output file, execution path, and notation correspondence.

### Included yeast data

The repository includes precomputed yeast `.dat` outputs. The 32 motif-level total-MI files in the `Yeast/` root are read directly by `Calculations-yeast.py`; corresponding copies are also stored beside the motif-level simulation scripts. The four `sMI-...-Yeast-num.dat` files are located in the `Yeast/` root. The E. coli output files are generated by running the E. coli simulation workflow.

## Correspondence between code and manuscript notation

### Parameters and regulatory functions

| Manuscript | CEM Python | Gillespie FFL | Gillespie open loop | Mathematica | Meaning |
|---|---|---|---|---|---|
| $\beta_X,\beta_Y,\beta_Z$ | `beta_x`, `beta_y`, `beta_z` | `b_x`, `b_y`, `b_z` | `b_x1`, `b_x2`, `b_y`, `b_z` | `βx`, `βy`, `βz` | Degradation rates |
| $\alpha_X,\alpha_Y,\alpha_Z$ | `alpha_y`, `alpha_z` are calculated locally; $\alpha_X=\beta_X\langle x\rangle$ | `a_x`, `a_y`, `a_z` | `a_x1`, `a_x2`, `a_y`, `a_z` | `αx`, `αy`, `αz` | Production-rate scales |
| $K_{XY},K_{XZ},K_{YZ}$ | `Kxy`, `Kxz`, `Kyz` | `Kxy`, `Kxz`, `Kyz` | `Kx2y`, `Kx1z`, `Kyz` | `Kxy`, `Kxz`, `Kyz` | Regulatory thresholds |
| $\langle x\rangle,\langle y\rangle,\langle z\rangle$ | `xbar`, `ybar`, `zbar` | `xav_param`, `yav_param`, `zav_param` | `x1av_param`, `x2av_param`, `yav_param`, `zav_param` | `xav`, `yav`, `zav` | Prescribed operating-point copy numbers |
| $f'_{YX},f'_{ZX},f'_{ZY}$ | `fyxp`, `fzxp`, `fzyp` | derivatives are not stored explicitly | `fyx2p`, `fzx1p`, `fzyp` in notebooks | `fyxp`, `fzxp`, `fzyp` | Local production sensitivities |
| signs of $X\to Y$, $X\to Z$, $Y\to Z$ | `(s_xy, s_xz, s_yz)` | encoded by `cal_func_f` | encoded by `cal_func_f` | encoded motif by motif | Edge signs in that order |

In the open-loop scripts, `x1` is the direct-path input $X$, whereas `x2` is the independent indirect-path input $\widetilde X$. Therefore, `Kx1z` corresponds to $K_{XZ}$ and `Kx2y` corresponds to the matched $K_{XY}$.

### Fluctuation and information quantities

| Manuscript | Main code name(s) | Meaning |
|---|---|---|
| $\eta_X^2$ | `eta_x`, `xcv2`, `x1cv2`, `x2cv2` | Normalized input variance |
| $\zeta_{XZ,d}$ | `eta_xz1` | Direct-path normalized covariance; the historical code name uses `eta` |
| $\zeta_{XZ,ind}$ | `eta_xz2` | Indirect-path normalized covariance; the historical code name uses `eta` |
| $\eta_{Z,0}^2$ | `eta_zi` | Intrinsic output contribution |
| $\eta_{Z,d}^2$ | `eta_zd` | Direct-path output contribution |
| $\eta_{Z,ind}^2$ | `eta_zind` | Indirect-path output contribution |
| $\eta_{Z,int}^2$ | `eta_zsyn` | Direct–indirect interference contribution; `syn` is a historical code label |
| $\eta_{Z,path}^2$ | `eta_p` | Pathway contribution to output fluctuation |
| $\eta_Z^2$ | `eta_z`, `zcv2` | Total normalized output variance |
| $\eta_{Z|X,path}^2$ | `denom_inner` | Conditional pathway fluctuation |
| $\eta_{Z|X,int}^2$ | `numer_inner` | Conditional interference fluctuation |
| $I(X;Z)$ | `ixz`, `Ixz` | Total input-output MI |
| $I_{\mathrm{path}}(X;Z)$ | `ixzp`, `Ixz0`, `Ix1x2z` | Pathway MI; in the open loop, $I(X,\widetilde X;Z)$ |
| $I_{\mathrm{int}}(X;Z)$ | `ixzs`, `Ixzs`, `sMI` | Interference mutual information |
| $S_{XZ}=2\zeta_{XZ,d}\zeta_{XZ,ind}/\eta_X^2$ | `Df` in the FFL notebooks | Pathway-interference strength |

## Exact code inventory

<details>
<summary><strong>Theory notebooks (5 files)</strong></summary>

- `Theory/FFL/AND.nb`
- `Theory/FFL/OR.nb`
- `Theory/Moments-calculations-FFL-OpenLoop.nb`
- `Theory/OL/Open-loop-AND.nb`
- `Theory/OL/Open-loop-OR.nb`

</details>

<details>
<summary><strong>CEM optimization: IMI (8 files)</strong></summary>

- `Numerics/CEM-optimization-IMI/imi-coherent-ecoli-and.py`
- `Numerics/CEM-optimization-IMI/imi-coherent-ecoli-or.py`
- `Numerics/CEM-optimization-IMI/imi-coherent-yeast-and.py`
- `Numerics/CEM-optimization-IMI/imi-coherent-yeast-or.py`
- `Numerics/CEM-optimization-IMI/imi-incoherent-ecoli-and.py`
- `Numerics/CEM-optimization-IMI/imi-incoherent-ecoli-or.py`
- `Numerics/CEM-optimization-IMI/imi-incoherent-yeast-and.py`
- `Numerics/CEM-optimization-IMI/imi-incoherent-yeast-or.py`

</details>

<details>
<summary><strong>CEM optimization: pathway MI (8 files)</strong></summary>

- `Numerics/CEM-optimization-path-MI/path-mi-coherent-ecoli-and.py`
- `Numerics/CEM-optimization-path-MI/path-mi-coherent-ecoli-or.py`
- `Numerics/CEM-optimization-path-MI/path-mi-coherent-yeast-and.py`
- `Numerics/CEM-optimization-path-MI/path-mi-coherent-yeast-or.py`
- `Numerics/CEM-optimization-path-MI/path-mi-incoherent-ecoli-and.py`
- `Numerics/CEM-optimization-path-MI/path-mi-incoherent-ecoli-or.py`
- `Numerics/CEM-optimization-path-MI/path-mi-incoherent-yeast-and.py`
- `Numerics/CEM-optimization-path-MI/path-mi-incoherent-yeast-or.py`

</details>

<details>
<summary><strong>CEM optimization: total MI (8 files)</strong></summary>

- `Numerics/CEM-optimization-total-MI/total-mi-coherent-ecoli-and.py`
- `Numerics/CEM-optimization-total-MI/total-mi-coherent-ecoli-or.py`
- `Numerics/CEM-optimization-total-MI/total-mi-coherent-yeast-and.py`
- `Numerics/CEM-optimization-total-MI/total-mi-coherent-yeast-or.py`
- `Numerics/CEM-optimization-total-MI/total-mi-incoherent-ecoli-and.py`
- `Numerics/CEM-optimization-total-MI/total-mi-incoherent-ecoli-or.py`
- `Numerics/CEM-optimization-total-MI/total-mi-incoherent-yeast-and.py`
- `Numerics/CEM-optimization-total-MI/total-mi-incoherent-yeast-or.py`

</details>

<details>
<summary><strong>Gillespie aggregation scripts (2 files)</strong></summary>

- `Numerics/Gillespie-after-optimization-IMI/Ecoli/Calculations-ecoli.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/Calculations-yeast.py`

</details>

<details>
<summary><strong>Gillespie FFL scripts (32 files)</strong></summary>

- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-c1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-c2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-c3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-c4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-i1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-i2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-i3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/AND/ecoli-ffl-i4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-c1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-c2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-c3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-c4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-i1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-i2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-i3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/FFL/OR/ecoli-ffl-i4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-c1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-c2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-c3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-c4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-i1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-i2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-i3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/AND/yeast-ffl-i4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-c1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-c2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-c3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-c4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-i1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-i2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-i3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/FFL/OR/yeast-ffl-i4-or.py`

</details>

<details>
<summary><strong>Gillespie open-loop scripts (32 files)</strong></summary>

- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-c1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-c2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-c3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-c4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-i1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-i2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-i3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/AND/ecoli-ol-i4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-c1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-c2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-c3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-c4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-i1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-i2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-i3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Ecoli/OL/OR/ecoli-ol-i4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-c1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-c2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-c3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-c4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-i1-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-i2-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-i3-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/AND/yeast-ol-i4-and.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-c1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-c2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-c3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-c4-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-i1-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-i2-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-i3-or.py`
- `Numerics/Gillespie-after-optimization-IMI/Yeast/OL/OR/yeast-ol-i4-or.py`

</details>


## Reproducibility and output conventions

- The CEM scripts use deterministic NumPy generators with fixed seeds. Individual runs use seeds `100 + 10 × run_id` after a pre-scan initialized with seed 42.
- The Gillespie scripts initialize NumPy's global random generator with a fixed, context-specific seed. The 64 simulation seeds are unique across organism, circuit, gate, and motif.
- CEM `.dat` files are tab-delimited and begin with comment/header lines prefixed by `#`.
- Gillespie `.dat` files contain a motif label and one numerical value.
- All generated paths are relative to the current working directory. Follow the working-directory instructions above to avoid mixing or overwriting outputs.
- The complete CEM runs and long Gillespie simulations can require substantial computation time. No parallel execution is implemented in the supplied scripts.

## Citation

If you use this repository, please cite the associated manuscript:

> M. Nandi, S. Chattopadhyay, and S. K. Banik, “An information-theoretic framework for understanding feed-forward loop abundances in transcriptional networks.”

Update this section with the journal, year, volume, pages/article number, and DOI after publication.

## License

No software license was included in the supplied archive. Add a `LICENSE` file before public release if reuse terms are to be specified.

## Contact

Questions about the codes can be directed to:

- Mintu Nandi: `mintunandi@ubi.s.u-tokyo.ac.jp`
- Sudip Chattopadhyay: `sudip@chem.iiests.ac.in`
- Suman K Banik: `skbanik@jcbose.ac.in`
