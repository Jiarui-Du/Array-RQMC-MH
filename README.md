# Array-RQMC-MH

This repository provides the MATLAB code, input data and saved results needed to reproduce the experiments, tables and figures in the Array-RQMC Metropolis--Hastings paper.

## Requirements

Tested with MATLAB R2024b Update 3. Statistics and Machine Learning Toolbox is required. Parallel Computing Toolbox is required for parallel execution; set `useParallel=false` or `UseParallel=false` to run serially. Open the repository root as the MATLAB Current Folder.

## Reproduce tables and figures

Plotting and summary scripts read saved results without running simulations.

| Entry | Output |
| --- | --- |
| `plot_1d` | One-dimensional variance figures |
| `plot_md` | Gaussian and Bayesian projection figures; select `gaussian`, `bayesian`, `ripley`, `pima`, or `all` |
| `plot_1d_horizon` | Fixed-particle-count horizon variance and scaled-variance figure |
| `vrf_1d` | One-dimensional VRF, slope intervals and TAE; optional CSV/TeX exports |
| `vrf_md` | Gaussian/Bayesian VRF, slope intervals and TAE |

Figures are saved as PDF and editable MATLAB FIG files under `figures/`. Summary scripts display results by default; set nonempty export paths to save CSV or LaTeX tables. An empty result path selects a complete saved result automatically; specify `matFile` to select a particular run.

## Run experiments

| Entry | Experiment |
| --- | --- |
| `run_1d` | Four one-dimensional targets with CDF, $x$ and $x^2$ tests |
| `run_gaussian_mwg` | Product Gaussian coordinatewise MWG with the projection test |
| `run_gaussian_joint` | Joint Gaussian Hilbert MH with the projection test |
| `run_bayesian` | Ripley and Pima logistic regression with the projection test |
| `run_1d_horizon` | Gaussian mixture/CDF with fixed $N$ and increasing $T$ |

Edit parameters at the start of each script. In `run_bayesian`, `dataset='all'` runs Ripley then Pima; choose either dataset to run it alone. The default settings are full experiments and can take many hours.

The main experiments use $N=2^{10},\ldots,2^{18}$, $R=200$, $T_0=8$, $T=128$, and master seed 2026. Estimates average transitions $T_0+1,\ldots,T$. MWG counts complete sweeps. The one-dimensional and MWG proposal standard deviation is $2.38$; joint Gaussian and Bayesian proposals use scale $2.38/\sqrt p$.

Gaussian dimensions are $p\in\{2,4,8\}$. Ripley and Pima have $p=3$ and $p=8$, respectively. The Gaussian projection is $f_p(\boldsymbol{x})=p^{-1/2}\sum_{j=1}^p x_j$. Bayesian chains start at the MLE and use

$$
f_{\mathrm{proj}}(\boldsymbol\beta)=\frac{1}{\sqrt p}\sum_{j=1}^{p}\frac{\beta_j-\widehat\beta_j}{\sqrt{\widehat\Sigma_{jj}}},
$$

where $\widehat\Sigma$ is the inverse observed information at the MLE.

Hilbert indices use $B_p=\min\{20,\lfloor64/p\rfloor\}$ bits per coordinate and exact `uint64` keys. Equal keys retain input order. Each transition uses a fresh scrambling.

The horizon entry defaults to $N\in\{2^{12},2^{14},2^{16}\}$, $T_0=128$, $T\in\{2^8,\ldots,2^{16}\}$, $R=50$, $\sigma=2.38$, and seed 2026. It studies variance against $L=T-T_0$. Saved MAT files contain the parameters used for each result. Set `resumeDir` to resume a horizon run.

## Statistical summaries and timing

Gaussian and Bayesian statistical stages use one combined `parfor` per $N$, followed by a checkpoint. One-dimensional tasks compute both independent methods per macro replication. Random seeds are assigned independently of worker scheduling. Progress reports elapsed time per grid point and total elapsed time.

Slope estimates use the largest five grid points. The $95\%$ percentile bootstrap intervals resample whole macro-replication estimates according to the random-stream pairing. Only time-average estimates are used in paper summaries.

Independent serial benchmarks use 20 complete single-test runs per method. The efficiency measure is

$$
\mathrm{TAE}=\mathrm{VRF}\,\frac{\operatorname{median}(t_{\mathrm{IID}})}{\operatorname{median}(t_{\mathrm{Array}})}.
$$

Timing includes random-stream setup, initialization, sorting, scrambling, Hilbert-index construction, transitions and test evaluation. Data/model preparation and file I/O are excluded. Statistical parallel wall times are not used for TAE. Timing results depend on hardware and system load.

The `runs/` directory contains only MAT files. Variances, VRFs, slope intervals and TAE tables are computed from these files; pre-exported CSV and LaTeX tables are not required.

## Directory layout

- `src/`: simulation, analysis and plotting functions.
- `data/`: model inputs, source datasets and preprocessing records.
- `runs/`: saved estimates and experiment metadata, organized by experiment.
- `figures/`: PDF and MATLAB FIG files.
- `tests/`: small validation functions.

For each Bayesian dataset, `data/<dataset>/author` contains the exact model inputs and `source` contains reconstruction data and row maps. See `data/README.txt` and `data/source_manifest.json` for sources and preprocessing. Result readers resolve local files rather than relying on paths recorded on the original computation machine.

## Small validation

```matlab
addpath('tests');
validate_public_package;
validate_1d_horizon;
```

These checks use temporary outputs and small grids ($N\le32$, $T\le16$, $R=4$), including two-worker consistency for the horizon experiment. They do not start a full experiment.
