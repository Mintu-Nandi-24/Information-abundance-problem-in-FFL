# An information-theoretic perspective on feed-forward loop abundances in transcriptional networks

This repository contains the analytical notebooks and numerical codes accompanying the manuscript **"An information-theoretic perspective on feed-forward loop abundances in transcriptional networks"** by Mintu Nandi, Sudip Chattopadhyay, and Suman K. Banik.

The implementation covers the steady-state linear noise approximation (LNA), decomposition of input-output mutual information (MI) into pathway MI and interference MI (IMI), representative-parameter calculations, abundance-guided cross-entropy-method (CEM) optimization, matched open-loop reference calculations, and Gillespie simulations at both representative and optimized operating points.

## Authors and contact

**Mintu Nandi - corresponding author**  
Universal Biology Institute, The University of Tokyo, 7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033, Japan  
Email: `mintunandi@ubi.s.u-tokyo.ac.jp`

**Sudip Chattopadhyay**  
Department of Chemistry, Indian Institute of Engineering Science and Technology, Shibpur, Howrah 711103, India  
Email: `sudip@chem.iiests.ac.in`

**Suman K. Banik**  
Department of Chemical Sciences, Bose Institute, EN 80, Sector V, Bidhan Nagar, Kolkata 700091, India  
Email: `skbanik@jcbose.ac.in`

## Contents and organization

The source distribution contains **32 code files: 7 Mathematica notebooks and 25 Python scripts**. The Python programs comprise 24 CEM scripts and **one consolidated stochastic-simulation script**. This simulation script runs both FFLs and their open-loop equivalents (OLEs), computes their information difference, and writes the final mean and standard-error files directly. No separate aggregation script is required.

| Location | Contents |
|---|---|
| [`Theory/Moments-calculations-FFL-OpenLoop.nb`](Theory/Moments-calculations-FFL-OpenLoop.nb) | Symbolic steady-state covariance calculations for the FFL and its open-loop reference. |
| [`Theory/Pre-optimization/`](Theory/Pre-optimization/) | Two notebooks evaluating total MI, pathway MI, and IMI at representative parameters, before CEM optimization. |
| [`Theory/FFL/`](Theory/FFL/) | Two notebooks evaluating optimized FFLs, interference strength, pathway sensitivities, and parameter dependence. |
| [`Theory/OL/`](Theory/OL/) | Two notebooks evaluating the joint input-output information of matched optimized open-loop references. |
| [`Numerics/CEM-optimization-IMI/`](Numerics/CEM-optimization-IMI/) | Eight IMI optimization scripts. |
| [`Numerics/CEM-optimization-path-MI/`](Numerics/CEM-optimization-path-MI/) | Eight pathway-MI comparison scripts. |
| [`Numerics/CEM-optimization-total-MI/`](Numerics/CEM-optimization-total-MI/) | Eight total-MI comparison scripts. |
| [`Numerics/Stochastic simulation/stochastic_imi_2.py`](Numerics/Stochastic%20simulation/stochastic_imi_2.py) | A standalone Gillespie program for all motifs, both gates, and representative and optimized parameters. |

Place this `README.md` at the repository root, alongside `Theory/` and `Numerics/`. The source archive contains no generated `.dat` result files or standalone figure exports; these are created by running the programs. Mathematica notebooks include saved evaluation outputs.

### Complete source directory tree

The paths below reproduce the supplied source structure. `README.md` is the documentation added at the root.

```text
.
|-- README.md
|-- Numerics/
|   |-- CEM-optimization-IMI/
|   |   |-- imi-coherent-ecoli-and.py
|   |   |-- imi-coherent-ecoli-or.py
|   |   |-- imi-coherent-yeast-and.py
|   |   |-- imi-coherent-yeast-or.py
|   |   |-- imi-incoherent-ecoli-and.py
|   |   |-- imi-incoherent-ecoli-or.py
|   |   |-- imi-incoherent-yeast-and.py
|   |   `-- imi-incoherent-yeast-or.py
|   |-- CEM-optimization-path-MI/
|   |   |-- path-mi-coherent-ecoli-and.py
|   |   |-- path-mi-coherent-ecoli-or.py
|   |   |-- path-mi-coherent-yeast-and.py
|   |   |-- path-mi-coherent-yeast-or.py
|   |   |-- path-mi-incoherent-ecoli-and.py
|   |   |-- path-mi-incoherent-ecoli-or.py
|   |   |-- path-mi-incoherent-yeast-and.py
|   |   `-- path-mi-incoherent-yeast-or.py
|   |-- CEM-optimization-total-MI/
|   |   |-- total-mi-coherent-ecoli-and.py
|   |   |-- total-mi-coherent-ecoli-or.py
|   |   |-- total-mi-coherent-yeast-and.py
|   |   |-- total-mi-coherent-yeast-or.py
|   |   |-- total-mi-incoherent-ecoli-and.py
|   |   |-- total-mi-incoherent-ecoli-or.py
|   |   |-- total-mi-incoherent-yeast-and.py
|   |   `-- total-mi-incoherent-yeast-or.py
|   `-- Stochastic simulation/
|       `-- stochastic_imi_2.py
`-- Theory/
    |-- FFL/
    |   |-- AND.nb
    |   `-- OR.nb
    |-- OL/
    |   |-- Open-loop-AND.nb
    |   `-- Open-loop-OR.nb
    |-- Pre-optimization/
    |   |-- AND.nb
    |   `-- OR.nb
    `-- Moments-calculations-FFL-OpenLoop.nb
```

Directory and filename capitalization should be preserved. In particular, **`Stochastic simulation` contains a space**. Quote that path in shell commands. Use the actual filename `stochastic_imi_2.py`; the shorter filename in the script's introductory usage examples is not the distributed filename.

## Software requirements

| Component | Requirements |
|---|---|
| CEM calculations | Python 3, NumPy, and Matplotlib. The study's author-reported Python version is **3.12.3**. |
| Stochastic calculations | Python standard library only; **Numba is optional** and is used automatically when importable. |
| Symbolic, optimized-FFL, and OLE notebooks | Notebook metadata records **Wolfram Mathematica 14.3**. |
| Representative-parameter notebooks | `Theory/Pre-optimization/AND.nb` and `OR.nb` record **Wolfram 15.0**. |

Install the Python packages needed by the CEM scripts:

```bash
python -m pip install numpy matplotlib
```

To enable optional compilation of the stochastic loops:

```bash
python -m pip install numba
```

The archive does not contain a pinned dependency file or record NumPy, Matplotlib, or Numba versions. Record the versions and whether Numba was enabled when reproducing results. Notebook version metadata identifies the saved front-end version; it is not a separate record of every historical calculation environment.

## Model and information quantities

Each FFL contains an input regulator $X$, an intermediate regulator $Y$, and an output $Z$. The direct route is $X\to Z$; the indirect route is $X\to Y\to Z$. The eight motifs differ in their activating or repressing regulatory signs. The OLE replaces the shared input of the indirect route by an independent, statistically equivalent input $\widetilde X$.

The numerical programs use the noncooperative Hill factors

$$
h_A(u;K)=\frac{u}{K+u},\qquad h_R(u;K)=\frac{K}{K+u}.
$$

The AND output multiplies the two regulatory factors. The OR output **adds** them, as in the source model; it is not a probabilistic Boolean-OR formula.

The prescribed operating-point vector is

$$
\boldsymbol\theta=(\beta_X,\beta_Y,\beta_Z,K_{XY},K_{XZ},K_{YZ},\langle x\rangle,\langle y\rangle,\langle z\rangle).
$$

Removal rates are in inverse minutes and thresholds and mean copy numbers are in molecules per cell. The production scales are determined from mean-field steady-state balance, rather than sampled independently. They remain fixed during each stochastic trajectory.

The information decomposition is

$$
I_{\mathrm{FFL}}(X;Z)=I_{\mathrm{path}}(X;Z)+I_{\mathrm{int}}(X;Z).
$$

Within the manuscript's matched LNA/Gaussian construction,

$$
I_{\mathrm{path}}(X;Z)=I_{\mathrm{OLE}}(X,\widetilde X;Z).
$$

The simulation program therefore reports the FFL information, the OLE joint-input information, and their signed difference. All information values are in **bits**. The simulated information quantities use the Gaussian formula evaluated with stochastic first and second moments; they are not histogram estimates of the full discrete MI.

## Workflow 1: analytical calculations in Mathematica

Open each notebook in Mathematica and evaluate its input cells in order. The calculations use parameter values embedded in the notebooks; no external data file or preceding CEM run is required to evaluate the supplied parameter sets. Save the notebook before evaluating cells that use `NotebookDirectory[]`.

The notebooks contain `Quit[]` commands. Evaluate the calculation groups in order and re-evaluate the required definitions after a kernel restart; an export cell needs the variables created by its preceding calculation cells.

### 1.1 Symbolic moments

Open:

```text
Theory/Moments-calculations-FFL-OpenLoop.nb
```

This notebook derives steady-state variances and covariances from the FFL and open-loop Lyapunov equations. It displays the symbolic results in the notebook and has no standalone file-export commands.

### 1.2 Representative-parameter calculations

Open:

```text
Theory/Pre-optimization/AND.nb
Theory/Pre-optimization/OR.nb
```

Both notebooks use the representative vector

```text
(0.1, 1.0, 10.0, 100.0, 100.0, 100.0, 40.0, 100.0, 100.0)
```

in the parameter order defined above, shared across all eight motifs. These are the representative-parameter calculations associated with the main-text AND-gate comparison and its supplementary OR-gate counterpart. No CEM search is performed in these notebooks.

Each notebook displays information plots and exports six motif/value tables:

```text
IMI-coherent-<GATE>.dat
IMI-incoherent-<GATE>.dat
MI-coherent-<GATE>.dat
MI-incoherent-<GATE>.dat
PMI-coherent-<GATE>.dat
PMI-incoherent-<GATE>.dat
```

Here `<GATE>` is `AND` or `OR`; `MI` denotes total MI, `PMI` denotes pathway MI, and `IMI` denotes signed interference MI.

**Export-path setup is required for these two notebooks.** Their six `Export` calls currently point to the author's local Windows directory:

```text
F:\Colab\FFL-Info-Decomp\Theory\non-CEM-IMI\
```

That path is not a repository directory. Replace the destination in each of those calls with a writable path on your computer. A portable replacement for the six export calls, after the tables `t1` through `t6` have been calculated, is:

```wolfram
gate = "AND"; (* Set to "OR" in Theory/Pre-optimization/OR.nb. *)
exportDir = NotebookDirectory[];
Export[FileNameJoin[{exportDir, "IMI-coherent-" <> gate <> ".dat"}], t1, "Table"];
Export[FileNameJoin[{exportDir, "IMI-incoherent-" <> gate <> ".dat"}], t2, "Table"];
Export[FileNameJoin[{exportDir, "MI-coherent-" <> gate <> ".dat"}], t3, "Table"];
Export[FileNameJoin[{exportDir, "MI-incoherent-" <> gate <> ".dat"}], t4, "Table"];
Export[FileNameJoin[{exportDir, "PMI-coherent-" <> gate <> ".dat"}], t5, "Table"];
Export[FileNameJoin[{exportDir, "PMI-incoherent-" <> gate <> ".dat"}], t6, "Table"];
```

This explicitly changes only the export destination: it places the generated tables beside the representative-parameter notebooks. It is a documented portability adjustment, not a change already applied to the distributed notebooks.

### 1.3 Optimized FFL calculations and biophysical quantities

Open:

```text
Theory/FFL/AND.nb
Theory/FFL/OR.nb
```

These notebooks evaluate the hard-coded optimized FFL operating points for both organisms. They calculate the information decomposition and quantities labeled `Df`, `Ed`, `Ei`, and robustness in the code. Their export blocks create:

```text
Theory/FFL/Data_exports_AND/
Theory/FFL/Data_exports_OR/
```

Each notebook exports **20 `.dat` tables**: five tables for each of four result families. Four tables separate class and organism; the fifth ends in `_ALL.dat` and combines them.

| Prefix/family | Contents |
|---|---|
| `00_Ipath_...` | Analytical pathway information. |
| `01_Df_vs_IMI_...` | Interference strength `Df` and IMI. |
| `02_Ed_Ei_...` | Direct- and indirect-path sensitivity quantities defined in the notebook. |
| `03_Robustness_vs_IMI_...` | The notebook's parameter-robustness score and IMI. |

For example, the AND notebook writes `01_Df_vs_IMI_CFFL_AND_Ecoli.dat`, the corresponding `IFFL` and `Yeast` tables, and `01_Df_vs_IMI_AND_ALL.dat`. The same naming convention applies to the other families and to `OR`.

It also saves the following three figures in **both PDF and PNG** format:

```text
Fig_Df_vs_IMI_<GATE>
Fig_Ed_Ei_<GATE>
Fig_Robustness_vs_IMI_<GATE>
```

These export directories are calculated from `NotebookDirectory[]`, not from the terminal's working directory. Saved notebook output can show a former Windows location; the actual export code uses the notebook's current location.

### 1.4 Matched open-loop information

Open:

```text
Theory/OL/Open-loop-AND.nb
Theory/OL/Open-loop-OR.nb
```

The two inputs are represented separately, with the direct input corresponding to $X$ and the indirect input to $\widetilde X$. These notebooks calculate joint input-output information using the supplied optimized FFL-matched parameters.

The export directories are:

```text
Theory/OL/Data_exports_OL_AND/
Theory/OL/Data_exports_OL_OR/
```

Each contains five `.dat` tables:

```text
00_Ic_CFFL_OL_<GATE>_Ecoli.dat
00_Ic_CFFL_OL_<GATE>_Yeast.dat
00_Ii_IFFL_OL_<GATE>_Ecoli.dat
00_Ii_IFFL_OL_<GATE>_Yeast.dat
00_Ic_Ii_OL_<GATE>_ALL.dat
```

The tables provide the analytical OLE information for comparison with the FFL pathway information. No additional optimization is performed by these notebooks.

## Workflow 2: abundance-guided CEM optimization

### 2.1 Running one condition

Each of the three CEM directories contains eight standalone scripts covering two organisms, two motif classes, and two gates. Their exact filenames are listed in the source tree.

**Run each script from its own CEM directory.** Output paths are relative to the working directory, and different metric directories reuse the same condition-specific output names.

Starting from the repository root, run the IMI example with:

```bash
cd Numerics/CEM-optimization-IMI
python imi-coherent-ecoli-and.py
```

Starting separately from the repository root, run the pathway-MI or total-MI counterpart with:

```bash
cd Numerics/CEM-optimization-path-MI
python path-mi-coherent-ecoli-and.py
```

```bash
cd Numerics/CEM-optimization-total-MI
python total-mi-coherent-ecoli-and.py
```

These programs have no command-line option parser. Edit the global configuration in the selected script to change the optimization budget, bounds, target windows, or sensitivity settings.

Matplotlib plots are displayed using `plt.show()` and are **not saved automatically by the Python CEM programs**. On a noninteractive system, select the `Agg` backend before execution. From the relevant CEM directory:

```bash
MPLBACKEND=Agg python imi-coherent-ecoli-and.py
```

The corresponding PowerShell commands are:

```powershell
$env:MPLBACKEND = "Agg"
python .\imi-coherent-ecoli-and.py
```

### 2.2 Common numerical settings

| Setting | Source value |
|---|---:|
| Independent CEM runs, `N_RUNS` | 100 |
| Pre-scan samples, `N_SAMPLES_FEASIBILITY` | 1,000 |
| Iterations per run, `N_ITER_CEM` | 1,000 |
| Population size, `POPULATION_SIZE` | 1,000 |
| Elite fraction, `ELITE_FRAC` | 0.1 |
| Smoothing coefficient, `ALPHA` | 0.3 |
| Covariance floor, `COV_FLOOR` | $10^{-6}$ |
| Feasibility tolerance, `FEASIBILITY_TOL` | $10^{-12}$ |
| Infeasible-loss offset, `INFEASIBLE_OFFSET` | $10^6$ |
| Pre-scan random seed, `RANDOM_SEED` | 42 |
| CEM seed for zero-based run `r` | `100 + 10*r` |
| Sensitivity perturbations | 200 evenly spaced relative perturbations over $[-0.1,0.1]$ |

The pre-scan selects an initial parameter center that is reused across the independent-seed CEM runs. One script's full CEM budget is $100\times1000\times1000=10^8$ candidate vectors, each involving four motif evaluations, before additional analyses. Running all 24 scripts is therefore a substantial calculation; these are not short demonstration defaults.

Within a class, all four motifs share the same candidate parameter vector. The search prioritizes satisfaction of the prescribed ratio windows, then maximizes the average metric in the feasible region. IMI scripts require all coherent values to be positive or all incoherent values to be negative, respectively; incoherent IMI ratios use magnitudes. These are implemented search restrictions, not universal information-sign claims. The control scripts evaluate their own information quantity even where their inherited comments or variable names still say `IMI`.

### 2.3 Parameter windows

The order is $(\beta_X,\beta_Y,\beta_Z,K_{XY},K_{XZ},K_{YZ},\langle x\rangle,\langle y\rangle,\langle z\rangle)$.

| Window | Lower bounds | Upper bounds |
|---|---|---|
| *E. coli* | `(0.005, 0.005, 0.005, 5, 5, 5, 10, 10, 10)` | `(0.1, 0.1, 0.5, 100, 100, 100, 100, 100, 100)` |
| Yeast | `(0.001, 0.001, 0.001, 20, 20, 20, 50, 50, 50)` | `(0.1, 0.1, 0.1, 500, 500, 500, 500, 500, 500)` |

### 2.4 Target windows as implemented

The three entries in each cell below give the intervals for motif 1 divided by motifs 2, 3, and 4, respectively. Thus, a coherent row lists `C1/C2; C1/C3; C1/C4`, and an incoherent row lists `I1/I2; I1/I3; I1/I4`.

| Organism | Class | Gate | IMI targets | Pathway-MI targets | Total-MI targets |
|---|---|---|---|---|---|
| E. coli | coherent | AND | [8, 10]; [2, 3]; [4, 5] | [8, 10]; [2, 3]; [4, 5] | [8, 10]; [2, 3]; [4, 5] |
| E. coli | coherent | OR | [8, 10]; [2, 3]; [4, 5] | [8, 10]; [2, 3]; [4, 5] | [8, 10]; [2, 3]; [4, 5] |
| E. coli | incoherent | AND | [8, 10]; [3.5, 4.5]; [4.5, 5.5] | [8, 10]; [3.5, 4.5]; [4.5, 5.5] | [8, 10]; [3.5, 4.5]; [4.5, 5.5] |
| E. coli | incoherent | OR | [8, 10]; [3.5, 4.5]; [4.5, 5.5] | [8, 10]; [3.5, 4.5]; [4.5, 5.5] | [8, 10]; [3.5, 4.5]; [4.5, 5.5] |
| Yeast | coherent | AND | [3, 8]; [20, 150]; [20, 150] | [3, 8]; [20, 150]; [20, 150] | [3, 8]; [20, 150]; [20, 150] |
| Yeast | coherent | OR | [2, 9]; [10, 150]; [10, 150] | [3, 8]; [20, 150]; [20, 150] | [3, 8]; [20, 150]; [20, 150] |
| Yeast | incoherent | AND | [5, 7]; [18, 22]; [50, 100] | [5, 7]; [18, 22]; [50, 100] | [5, 7]; [18, 22]; [50, 100] |
| Yeast | incoherent | OR | [4, 8]; [15, 25]; [20, 150] | [5, 7]; [18, 22]; [50, 100] | [5, 7]; [18, 22]; [50, 100] |

The two yeast OR-gate rows use different target intervals between the IMI scripts and the control scripts. The table records the supplied configuration as written; it should not be interpreted as an identical-target comparison in those cases.

### 2.5 CEM output paths and files

A run creates `<class>_<organism>_<GATE>_data/` in its working directory. For example, when the command is run from the recommended IMI directory, the output location is:

```text
Numerics/CEM-optimization-IMI/coherent_ecoli_AND_data/
```

The filename prefix is `<class>-<organism>-<GATE>-`, for example `coherent-ecoli-AND-`. The export routine writes up to twelve tables:

| Filename suffix | Contents |
|---|---|
| `diagnostics.dat` | CEM center trajectories, objective diagnostics, and feasible fractions. |
| `best_penalty_convergence.dat` | Mean and sample SD of the best penalty versus iteration. |
| `imi_trajectory.dat` | Motif-wise best-so-far metric trajectories. |
| `final_imi_summary.dat` | Final motif-wise metric means and sample SDs. |
| `ratio_satisfaction_scatter.dat` | Run-level ratios, target bounds, and summaries. |
| `avg_imi_by_run.dat` | Best average metric and associated run diagnostics. |
| `pairwise_distances.dat` | Pairwise normalized log-parameter distances, when the distance array is nonempty. |
| `parameter_correlations.dat` | Parameter correlation matrix across the selected run results. |
| `normalized_sensitivity.dat` | Finite-perturbation normalized sensitivities. |
| `parameter_cv.dat` | Parameter coefficients of variation across selected results. |
| `violin_data.dat` | Selected run-level parameter vectors. |
| `final_parameter_summary.dat` | Mean and sample SD of the selected parameter vectors. |

The files are tab-delimited, use `#` comment/header lines, and can be imported into OriginPro. Names such as `imi_trajectory` and `final_imi_summary` are retained in all three directories. **The directory identifies the actual information measure**: IMI, pathway MI, or total MI.

Across-run SDs use `ddof=1`, not the standard error. Final summaries normally use feasible runs. If none are feasible, the programs print a warning and fall back to all runs for these summaries; a file labeled `final_imi_summary` or `violin_data` does not by itself establish feasibility.

The run-distance calculation uses

$$
q_{ri}=\frac{\ln\theta_{ri}-\ln L_i}{\ln H_i-\ln L_i},\qquad
D_{rr'}=\sqrt{\frac{1}{9}\sum_{i=1}^9(q_{ri}-q_{r'i})^2}.
$$

Sensitivity calculations operate around the selected mean parameter vector. They average absolute relative metric changes divided by the requested relative parameter perturbation; parameter values are clipped to the bounds. They are distinct from the parameter-robustness calculations implemented separately in the Mathematica notebooks.

CEM output is **not read automatically** by the analytical notebooks or stochastic script. To study a newly optimized point, update the embedded parameter values consistently in the downstream files. Retain the original values separately when reproducing the supplied operating points.

## Workflow 3: Gillespie simulations for representative and optimized parameters

### 3.1 Main command and working directory

From the **repository root**, run:

```bash
python "Numerics/Stochastic simulation/stochastic_imi_2.py" 100000 10
```

The same quoted command works in PowerShell. It uses the manuscript protocol of **$10^5$ input correlation times and 10 independent replicates** for every simulated system.

The two positional arguments are:

```text
python "Numerics/Stochastic simulation/stochastic_imi_2.py" NTAU NREP
```

`NTAU` is the total simulation duration in units of $1/\beta_X$, including the discarded initial transient. `NREP` is the number of replicates per motif and system. No other command-line selector is implemented: each run evaluates **all optimized and representative cases**.

| Setting | Program with no arguments | Manuscript-protocol command above |
|---|---:|---:|
| Input-correlation-time units, `NTAU` | **10,000** | **100,000** |
| Independent replicates, `NREP` | 10 | 10 |
| Discarded transient, `BURN_IN_TAU` | 20 | 20 |
| Total physical duration | `NTAU / beta_X` | $10^5/\beta_X$ |
| Accumulation duration | `(NTAU - 20) / beta_X` | $(10^5-20)/\beta_X$ |

**The distributed script's `DEFAULT_NTAU` is `10000.0`, not `100000.0`.** Pass `100000 10` explicitly for the protocol described in the manuscript; no source change is needed. Use `NTAU > 20` and at least two replicates for an estimated SEM. With one replicate, the code writes a zero SEM by convention, which is not evidence of negligible uncertainty.

The complete run is computationally intensive. A smaller execution test can be run with `1000 2`, but it must not be treated as reproducing the manuscript's numerical precision. All runs reuse the same output filenames, so move or copy existing results before changing the settings.

### 3.2 Conditions evaluated

The program first evaluates the eight entries of the `OPTIMIZED` dictionary, sorted by organism, gate, and class. It then evaluates the representative parameter vector for AND/coherent, AND/incoherent, OR/coherent, and OR/incoherent conditions.

For each condition it simulates four motifs, and each motif has both an FFL and an OLE trajectory per replicate. The complete run therefore contains **48 motif/parameter/gate comparisons** and, for 10 replicates, **960 trajectories**: 640 optimized and 320 representative trajectories.

The representative vector is shared by **all eight motifs and both gates**:

```text
REPRESENTATIVE = (0.1, 1.0, 10.0, 100.0, 100.0, 100.0, 40.0, 100.0, 100.0)
```

The optimized vectors below are embedded in `stochastic_imi_2.py`; the program does not require a CEM export file.

| Organism | Gate | Class | `(beta_X, beta_Y, beta_Z)` | `(K_XY, K_XZ, K_YZ)` | `(<x>, <y>, <z>)` |
|---|---|---|---|---|---|
| Ecoli | AND | coherent | `(0.0050007, 0.099884, 0.49986)` | `(63.52, 69.95, 18.4)` | `(10.64, 13.17, 99.94)` |
| Ecoli | AND | incoherent | `(0.0050023, 0.098712, 0.07257)` | `(54.32, 99.55, 37.29)` | `(15.93, 48.82, 56.96)` |
| Ecoli | OR | coherent | `(0.013426, 0.093002, 0.27098)` | `(64.61, 41.42, 26.9)` | `(10, 17.25, 99.95)` |
| Ecoli | OR | incoherent | `(0.034645, 0.062942, 0.26985)` | `(34.95, 44.2, 21.86)` | `(10.12, 39.82, 99.83)` |
| Yeast | AND | coherent | `(0.012336, 0.00904, 0.02863)` | `(467.23, 27.69, 499.8)` | `(50, 63.42, 499.96)` |
| Yeast | AND | incoherent | `(0.0012494, 0.062408, 0.02536)` | `(217.3, 157.6, 25.69)` | `(50.01, 446.66, 306.79)` |
| Yeast | OR | coherent | `(0.001223, 0.076122, 0.042529)` | `(134.93, 294.17, 401.08)` | `(51.13, 53.04, 357.06)` |
| Yeast | OR | incoherent | `(0.0010933, 0.044021, 0.088137)` | `(261.24, 272.19, 22.59)` | `(112.8, 431.18, 206.83)` |

### 3.3 Trajectories, initialization, and moment estimation

Each FFL has six production/removal channels for $X,Y,Z$. Its OLE has eight channels for $X,\widetilde X,Y,Z$, with $X$ driving the direct route and $\widetilde X$ the indirect route. The OLE input means, production rates, and removal rates are matched to those of the FFL input, while the two input reaction processes are independent.

Production scales are calculated once from the prescribed mean-field balance. Initial molecule counts are `int(round(xav))`, `int(round(yav))`, and `int(round(zav))`; both OLE inputs start at the rounded input mean. Each FFL/OLE pair runs for the same physical duration. Scaling by $1/\beta_X$ fixes the number of input correlation times, rather than guaranteeing identical precision for all parameter sets.

The code accumulates **residence-time-weighted** first and second moments of the state occupied during each inter-reaction interval. An interval crossing the burn-in boundary contributes only its post-burn-in portion, and the final interval is clipped to the requested stop time. No averaging over an equally weighted list of reaction events is used.

For each FFL replicate, the program converts the measured input/output moments to

$$
I_{\mathrm{FFL}}^{(r)}=
\frac12\log_2\!\left[
\frac{\eta_Z^2}{\eta_Z^2-\zeta_{XZ}^2/\eta_X^2}
\right].
$$

For each OLE replicate, it uses

$$
I_{\mathrm{OLE}}^{(r)}=
\frac12\log_2\!\left[
\frac{\eta_{Z,\mathrm{OLE}}^2}
{\eta_{Z,\mathrm{OLE}}^2-
 (\zeta_{XZ}^{\mathrm{OLE}})^2/\eta_{X,\mathrm{OLE}}^2-
 (\zeta_{\widetilde X Z}^{\mathrm{OLE}})^2/\eta_{\widetilde X,\mathrm{OLE}}^2}
\right].
$$

The OLE estimator sets the input cross-covariance to its model value of zero; it does not fit a finite-sample cross-covariance between the independent inputs. All normalized quantities use the simulated means and moments of the relevant replicate.

The replicate-wise interference estimate is

$$
I_{\mathrm{int}}^{(r)}=I_{\mathrm{FFL}}^{(r)}-I_{\mathrm{OLE}}^{(r)}.
$$

For each of the three quantities, the program reports the replicate mean and

$$
\mathrm{SEM}=\sqrt{\frac{1}{N_{\mathrm{rep}}(N_{\mathrm{rep}}-1)}
\sum_{r=1}^{N_{\mathrm{rep}}}(I^{(r)}-\overline I)^2}.
$$

For IMI, this formula is applied to the **replicate-wise differences**, not to a subtraction of separately reported error bars. The manuscript's sign-resolution check, $|\overline I_{\mathrm{int}}|>3\,\mathrm{SEM}$, can be assessed from the mean and `-sem` files; the script does not write a separate pass/fail flag.

The Gillespie trajectories follow the nonlinear propensities, but the information conversion remains Gaussian. The comparison with the analytical expressions assesses moment-based LNA predictions, not the exact discrete mutual information or a separate test of the Gaussian approximation.

### 3.4 Reproducible seeds and optional acceleration

For zero-based motif index `m` within a class and replicate index `r`, the seeds are:

```text
FFL seed = salt + 1000*m + 2*r
OLE seed = FFL seed + 1
```

For optimized condition index `k` in `sorted(OPTIMIZED.items())`, `salt = 7919*(k+1)`. For representative conditions in the gate/class order described above, `salt = 104729*(k+1)`.

When Numba is importable, `@njit(cache=True)` compiles the simulation loops. The console reports whether compilation is enabled. Without Numba, the same algorithm executes as ordinary Python. Record the backend and software versions for a rerun; do not assume bit-for-bit agreement between the compiled and pure-Python random-number implementations. Numba may create cache files beside the script.

### 3.5 Output directories and exact filename conventions

When run from the repository root, the program creates:

```text
results/
`-- stochastic/
    |-- optimized/        # 48 data files
    `-- representative/   # 24 data files
```

These are generated directories, not directories already present in the source archive. Running the command from another working directory places `results/stochastic/` there instead.

For each condition, six files contain the means and SEMs of three quantities:

| Prefix | Stored quantity |
|---|---|
| `IMI` | Signed replicate-wise FFL-minus-OLE information difference. |
| `Itotal` | Gaussian total MI estimated from FFL trajectories. |
| `Ipath` | Gaussian joint-input MI estimated from OLE trajectories, used as the simulation comparison for analytical pathway MI. |

Optimized filenames have the form

```text
<quantity>-<class>-<GATE>-<Organism>.dat
<quantity>-<class>-<GATE>-<Organism>-sem.dat
```

with `<quantity>` equal to `IMI`, `Itotal`, or `Ipath`, `<class>` equal to `coherent` or `incoherent`, `<GATE>` equal to `AND` or `OR`, and `<Organism>` equal to **`Ecoli` or `Yeast`**.

For example:

```text
results/stochastic/optimized/IMI-coherent-AND-Ecoli.dat
results/stochastic/optimized/IMI-coherent-AND-Ecoli-sem.dat
results/stochastic/optimized/Itotal-coherent-AND-Ecoli.dat
results/stochastic/optimized/Itotal-coherent-AND-Ecoli-sem.dat
results/stochastic/optimized/Ipath-coherent-AND-Ecoli.dat
results/stochastic/optimized/Ipath-coherent-AND-Ecoli-sem.dat
```

Representative filenames omit the organism:

```text
results/stochastic/representative/IMI-coherent-AND.dat
results/stochastic/representative/IMI-coherent-AND-sem.dat
results/stochastic/representative/Itotal-coherent-AND.dat
results/stochastic/representative/Itotal-coherent-AND-sem.dat
results/stochastic/representative/Ipath-coherent-AND.dat
results/stochastic/representative/Ipath-coherent-AND-sem.dat
```

The same pattern applies to incoherent motifs and OR regulation. A complete run writes **72 `.dat` files**. Each is a **headerless two-column, tab-separated** table containing `motif` and `value`, with four rows in `C1`--`C4` or `I1`--`I4` order. Values are written to six decimal places. Mean and SEM files have matching row order. Use them together for points and error bars.

The program also prints all three means and SEMs to the console. It does not save individual trajectories, replicate-level values, moment matrices, or plots. The filename does not encode run length or replicate count; preserve a run log with the results. For example:

```bash
python "Numerics/Stochastic simulation/stochastic_imi_2.py" 100000 10 > stochastic_run.log 2>&1
```

## Interpreting analytical and stochastic comparisons

The representative notebooks and stochastic `REPRESENTATIVE` setting address the pre-optimization AND comparison and its OR counterpart. The optimized FFL notebooks, the embedded `OPTIMIZED` dictionary, and the OLE notebooks address the optimized information patterns and matched-reference comparisons. Use the full motif/organism/class/gate identifiers, rather than relying only on manuscript figure numbers, which may change between versions.

For numerical figure preparation, `Itotal` supplies simulation points for total MI, `Ipath` supplies the OLE-based pathway-information comparison, and `IMI` supplies the signed interference comparison. Their `-sem` files supply stochastic standard errors. In contrast, CEM error columns describe across-run sample SDs. These two kinds of error bars should not be labeled interchangeably.

The analytical notebook parameters and stochastic dictionaries are embedded values, not live links to CEM results. After changing any parameter, update all corresponding calculations explicitly. A CEM mean over information values is also not generally equal to the information evaluated at the mean parameter vector.

## Code-to-manuscript notation

### Operating parameters and pathway variables

| Manuscript quantity | CEM scripts | Stochastic script | Mathematica notation |
|---|---|---|---|
| $\beta_X,\beta_Y,\beta_Z$ | `beta_x`, `beta_y`, `beta_z` | `bx`, `by`, `bz` | `\[Beta]x`, `\[Beta]y`, `\[Beta]z` |
| $K_{XY},K_{XZ},K_{YZ}$ | `Kxy`, `Kxz`, `Kyz` | `kxy`, `kxz`, `kyz` | `Kxy`, `Kxz`, `Kyz` |
| Prescribed mean copy numbers | `xbar`, `ybar`, `zbar` | `xav`, `yav`, `zav` in `motif_arguments` | `xav`, `yav`, `zav` |
| $\alpha_X,\alpha_Y,\alpha_Z$ | `alpha_y`, `alpha_z`; input scale implicit | `ax`, `ay`, `az` | `\[Alpha]x`, `\[Alpha]y`, `\[Alpha]z` where used |
| $f'_{YX},f'_{ZX},f'_{ZY}$ | `fyxp`, `fzxp`, `fzyp` | Not required by the reaction simulation | `fyxp`, `fzxp`, `fzyp` |
| Regulatory signs $(X\to Y,X\to Z,Y\to Z)$ | Motif dictionaries | `TOPOLOGY` | Assigned in each motif block |
| Independent indirect-route input $\widetilde X$ | OLE not simulated by the CEM scripts | `xt` | Second input in the OLE notebooks |

### Information and fluctuation variables

| Manuscript quantity | Source names or output |
|---|---|
| $\eta_X^2$ | `eta_x` in CEM; `ex` in stochastic loops. |
| $\zeta_{XZ,d},\zeta_{XZ,\mathrm{ind}}$ | `eta_xz1`, `eta_xz2` in CEM; `\[Eta]xz1`, `\[Eta]xz2` in notebooks. |
| $\eta_{Z,0}^2$, $\eta_{Z,d}^2$, $\eta_{Z,\mathrm{ind}}^2$ | `eta_zi`, `eta_zd`, `eta_zind` in CEM. |
| $\eta_{Z,\mathrm{int}}^2$ | `eta_zsyn`; `syn` is a historical code identifier. |
| $\eta_{Z,\mathrm{path}}^2$ | `eta_p`. |
| $\eta_{Z|X,\mathrm{path}}^2$ | `denom_inner` in the IMI calculation. |
| $\eta_{Z|X,\mathrm{int}}^2$ | `numer_inner` in the IMI calculation. |
| $I(X;Z)$ | Total-MI control metric; representative `MI-...` and stochastic `Itotal-...` files. |
| $I_{\mathrm{path}}(X;Z)$ | Pathway-MI control metric; representative `PMI-...` and optimized `00_Ipath_...` files. |
| $I_{\mathrm{OLE}}(X,\widetilde X;Z)$ | OLE notebooks and stochastic `Ipath-...` files. |
| $I_{\mathrm{int}}(X;Z)$ | IMI control metric; representative and stochastic `IMI-...` files. |
| $S_{XZ}=2\zeta_{XZ,d}\zeta_{XZ,\mathrm{ind}}/\eta_X^2$ | `Df` in the optimized FFL notebooks. |
| Direct/indirect sensitivity quantities | `Ed`, `Ei` and `02_Ed_Ei_...` notebook exports. |

Names such as `imi_analytical_*`, `imis`, `avg_imi`, and `imi_trajectory` persist in the total/pathway comparison code. Their meaning is determined by the script's calculation and directory, not the inherited variable name.

## Reproducibility and reuse

Preserve the source files, software versions, parameter values, seeds, command-line arguments, and generated outputs used for a reported result. Running a program again in the same location overwrites files with the same names. Mathematica export tables, Python CEM tables, and stochastic tables use different header conventions, as specified above.

The supplied code implements the model and numerical protocols documented here. This README does not report a new full CEM optimization, a new long-trajectory simulation run, or a new execution of the Mathematica notebooks. Source inspection and command/output documentation are distinct from regenerating the manuscript's numerical results.

Please cite the accompanying manuscript when using these calculations. Questions about the implementation or reuse can be directed to the corresponding author. No separate license file is included in this source distribution.
