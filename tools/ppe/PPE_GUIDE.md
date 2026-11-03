# Tuning CAM6 convection with a nudged-SCAM perturbed-parameter ensemble

**A guide to the science, the workflow, and the current state of the code**

> **DRAFT. Not yet reviewed.**
> This guide and the scripts in `tools/ppe/` were drafted by Claude (Anthropic's AI assistant, model Claude Opus 5.5) working with the project lead, and are being checked by the project lead. Until this notice is removed, treat the science summary, code descriptions and references as provisional: verify details against the code and the cited papers before relying on them, and report anything that looks wrong.

*Status as of 2026-11-03. Code lives in `tools/ppe/` on the `compassdev` branch of the CAM fork `jtruesdal/CAM-1`.*

---

## How to use this document

This document is written for someone joining the project who has not worked with SCAM, CAMIOP, nudging or the PPE scripts. It runs from the scientific motivation down to the commands, so it can be read top to bottom once and then used as a reference.

- **Sections 1–3** explain the goal, the theory behind the method, and the published work it builds on.
- **Section 4** gives the overall three-tier plan.
- **Section 5** describes what has been built and verified.
- **Sections 6–9** go through each tier in detail, including work that hasn't been built yet.
- **Sections 10–11** cover current status, bottlenecks and next steps.
- **The appendices** list shell variables, commands, pitfalls, code references and the literature.

Each step is marked **Done**, **In progress** or **Planned**, so it's clear what exists and what is still design.

---

## 1. Goal

The aim is to improve the parameter settings of CAM6's moist physics: deep convection (Zhang–McFarlane, "ZM"), the unified turbulence and shallow-convection scheme (CLUBB), and cloud microphysics (MG/PUMAS).

The usual way to tune a GCM is to run many multi-year climate simulations, which is expensive and mixes fast physics errors with slow feedbacks. Instead, we use thousands of cheap **single-column** simulations driven by the real atmosphere. In each one the model is held close to observations, and we measure how hard it has to be pushed to stay there. Parameter settings that need less pushing have physics that agrees better with the observed atmosphere. Promising settings are then checked in the full 3D model with progressively more expensive tests.

---

## 2. Scientific background

### 2.1 Single-column models and SCAM

A single-column model (SCM) runs a GCM's physics for one atmospheric column. The large-scale dynamics, mainly horizontal and vertical advection of heat, moisture and momentum, come from outside as prescribed forcing. Randall et al. (1996) describe SCMs as the link between field observations and climate models: they isolate the physics so it can be tested cheaply against observed cases.

SCAM is CAM's single-column configuration. SCAM6 and its uses are described in Gettelman et al. (2019). SCAM normally takes its forcing from an "IOP" file (intensive observing period), usually built from field-campaign observations such as ARM's.

### 2.2 CAMIOP: exact reproduction of a 3D column

CAM can also *write* an IOP file. With the `-camiop` configure option, a 3D run outputs, at every time step and for every column:

- the state (T, q and all tracers, u, v, surface pressure, ω)
- the dynamics tendencies of every prognostic variable (`divT3d`, `divq3d`, one `<tracer>_dten` per constituent, `divU3d`, `divV3d`)

This output is the **CAMIOP** file. When SCAM is run in "BFB mode" (`scm_cambfb_mode = .true.`, `scm_use_3dfrc = .true.`) from a CAMIOP file, at the same start time and with the same initial file, it reproduces the 3D model's column **bit-for-bit (BFB)**. Every value is identical, to the last bit, for as long as the run lasts.

That makes SCAM a microscope on the 3D model. Any column of a global simulation can be extracted and rerun alone, and we know exactly what the 3D model did there.

How the forcing is built and applied:

- **In the 3D run** (`src/dynamics/se/advect_tend.F90`, called from `stepon.F90`), the tendency is computed as a residual: (change of T over the dynamics step) minus (the physics forcing handed to the dycore). It therefore contains everything that isn't physics: advection, vertical remapping and hyperviscosity. It is computed within each time step, so it doesn't depend on earlier steps and survives model restarts.
- **In SCAM** (`advance_iop_forcing` in `src/dynamics/se/apply_iop_forcing.F90`), the tendency is added as a **fixed time series**: `T_new = T_old + dt·(physics tendency + divT3d)`. SCAM doesn't recompute advection from ω in this mode. T, q and the other tracers are **prognostic**: they start from the initial file (`ncdata`) and evolve under physics, the fixed forcing and nudging. The IOP file's `t` and `q` are never written back into the state; they're used only for the `TDIFF`/`QDIFF` diagnostics, which compare SCAM with the 3D column. Only the wind components are *replaced* by the 3D values every step (`u = uobs`, `v = vobs`).

### 2.3 Nudging

**Nudging** adds a relaxation term that pulls the model state toward a target state, here ERA5 reanalysis:

  `Nudge(X) = −(X − X_ERA5) / τ`

`τ` is the relaxation timescale. Nudging is widely used in CAM to constrain the large-scale circulation so different model versions can be compared under the same weather; see, for example, Kooperman et al. (2012) and Zhang et al. (2014).

In CAM, nudging is implemented in `src/physics/cam/nudging.F90` as an ordinary physics tendency at the end of the physics step. It also writes diagnostics:

- `Nudge_T`, `Nudge_Q` (and `Nudge_U`, `Nudge_V`): the applied tendencies
- `Target_T`, `Target_Q`, …: the ERA5 values being nudged toward

Our settings: U, V, T and Q are nudged with coefficient 0.1667 to hourly targets, which works out to **τ ≈ 6 hours**. Surface pressure is not nudged.

Because nudging is a *physics* tendency, it is **not** included in the CAMIOP dynamics forcing; it's subtracted out with the rest of the physics. A nudged SCAM run therefore recomputes nudging from its own state. With default parameters it reproduces the 3D nudging tendency exactly, because it is BFB.

### 2.4 Nudging tendencies as a measure of physics error

This is the core idea. With CAMIOP forcing and nudging, the SCAM temperature budget is

  `dT/dt = Phys(T; θ) + F(t) + Nudge(T)`

- `Phys(T; θ)` is the physics tendency, which depends on the parameters θ.
- `F(t)` is the fixed CAMIOP dynamics forcing, which doesn't depend on θ.
- `Nudge(T)` is the relaxation to ERA5.

Averaged over a long period, the left side is close to zero because the state stays near ERA5. So

  `⟨Nudge⟩ ≈ −⟨Phys(θ)⟩ − ⟨F⟩`

The time-mean nudging tendency is exactly what is needed to make up for the mismatch between what the physics does and what the observed atmosphere requires, given the large-scale forcing. If a parameter change makes the physics closer to the observed heating and moistening, the mean nudging tendency shrinks. With `F` fixed, **the change in mean nudging tendency between two parameter settings is approximately minus the change in their mean physics tendency.**

This is a single-column version of an established idea from numerical weather prediction:

- **Klinker and Sardeshmukh (1992)** diagnosed model errors from the systematic tendencies the model develops in the first time step of forecasts started from analyses.
- **Rodwell and Palmer (2007)** proposed the **initial-tendency** diagnostic. In a data-assimilation cycle, the mean analysis increment (the correction needed to bring the model back to observations) shows the model's systematic error. Process tendencies in the first hours of forecasts show which parameterization causes it.

Our continuous nudging tendency plays the role of the analysis increment. Nudging is a continuous, weak data assimilation, and the mean correction it makes is the systematic error. We call it the **SCAM analog of the Rodwell–Palmer initial-tendency diagnostic**.

### 2.5 What this metric can and can't tell us

**Strengths**

- **It attributes error to physics.** Under nudging, the state stays close to the real atmosphere, and the large-scale forcing comes from a 3D model that is itself nudged to ERA5. Systematic nudging is therefore mostly due to the physics.
- **It resolves the vertical profile.** The mean nudging profile shows *where* the physics heats or moistens too much or too little.
- **It's cheap and exact for the baseline.** The default member reproduces the 3D model, so any parameter-induced change is real and not an artifact of the single-column setup.

**Limitations**

1. **Fixed forcing means no circulation feedback.** In SCAM the dynamics tendencies are prescribed. A parameter change that would, in 3D, change the convective heating and hence the large-scale vertical motion can't do so here. Even vertical advection of the column's own temperature anomalies doesn't respond; nudging is the only restoring force. The SCAM stage is reliable for the *vertical structure, timing and partitioning* of heating and moistening, but can't test convection–circulation coupling. **This is why the 3D tiers are essential.**
2. **The result depends on τ.** The size of nudging tendencies depends on the relaxation timescale, so τ must be identical across the whole ensemble.
3. **Winds give no information in SCAM.** u and v are reset from the CAMIOP file every step, so `Nudge_U` and `Nudge_V` don't reflect the physics. Only T and Q are used.
4. **ERA5 isn't truth.** ERA5 has its own model-dependent errors, especially in humidity, the boundary layer and the upper troposphere. Some nudging can never be removed by any CAM parameter setting. Metrics are therefore judged *relative to the default member*, not against zero.
5. **Compensating errors are possible.** A parameter set can reduce T/Q nudging by having, say, radiative heating compensate for wrong convective heating while clouds get worse. **Independent observations** (Section 7.6) are needed to catch this.

### 2.6 Short-range errors and climate errors: CAPT and Transpose-AMIP

Under the Cloud-Associated Parameterizations Testbed (**CAPT**; Phillips et al. 2004), a climate model is run as a weather model: many short hindcasts are started from analyses, and the forecast errors are compared with observations. Because the large-scale state is close to the real atmosphere in the first days, errors that grow quickly can be traced to the fast physics. Xie et al. (2012) showed that many CAM systematic errors appear within the first few days of such hindcasts. The **Transpose-AMIP** project ran this kind of experiment across many models (Williams et al. 2013), and Ma et al. (2014) found that many CMIP5 climate errors are already present in short-range forecasts.

Two implications for this project:

- Fast-physics errors *can* be diagnosed in short, constrained runs.
- CAPT-style hindcasts are the natural 3D check on SCAM results, because they restore the circulation feedback that SCAM lacks.

### 2.7 Perturbed-parameter ensembles, emulators and history matching

A **perturbed-parameter ensemble** (PPE) runs a model many times with parameters sampled across their plausible ranges. Eidhammer et al. (2024) describe a large CAM6 PPE; our parameter ranges come from their Table 1.

Model runs are too expensive to search parameter space directly, so a statistical **emulator** is trained on the PPE: a fast surrogate, usually a Gaussian-process regression, that predicts a metric for any parameter setting together with its uncertainty. A common rule of thumb is about 10 runs per parameter for the first training design (Loeppky et al. 2009).

**History matching** (Williamson et al. 2013 applied it to climate models) uses the emulator to *rule out* parameter regions, rather than looking for a single "best" point. For each candidate θ and each metric, it computes an **implausibility**

  `I(θ) = |z − E[f(θ)]| / sqrt(Var_obs + Var_discrepancy + Var_emulator(θ))`

- `z` is the observed target.
- `E[f(θ)]` is the emulator's prediction.
- The denominator combines observational uncertainty, model discrepancy (error the model can't remove at any parameter setting) and emulator uncertainty.

Points with `I > 3` for any metric are ruled out. What remains is the **not-ruled-out-yet (NROY)** space. New simulations are then run inside it ("waves"), the emulators are refit, and the space shrinks further.

The **HighTune** project did exactly this for the LMDZ model:

- **Couvreux et al. (2021)** used single-column cases compared against LES to history-match boundary-layer and convection parameters.
- **Hourdin et al. (2021)** carried the reduced parameter space into global calibration.

Our design follows the same idea, with two differences: our single-column cases are *real-atmosphere columns forced by a nudged 3D model*, not idealized LES cases, and our main metric is the nudging tendency, supplemented by observations. For the broader practice of model tuning, see Hourdin et al. (2017).

---

## 3. Why this design

Putting Sections 2.4–2.7 together:

- **SCAM with CAMIOP plus nudging** gives thousands of cheap, exactly reproducible, physics-focused experiments at real places and times.
- **Nudging tendencies** give a vertically resolved, physically interpretable error signal that doesn't need observations at every column.
- **ARM sites and satellite data** guard against compensating errors.
- **Emulators and history matching** make the most of a limited number of runs and give honest uncertainty, not a single over-tuned answer.
- **Short-range 3D tests (CAPT)** check what SCAM can't, namely circulation feedback, before committing to expensive climate runs.

---

## 4. The three-tier plan

| Tier | What runs | Cost per run | What it decides |
|---|---|---|---|
| **0. Foundation** | Verification tests of the setup itself | Small | Whether SCAM results can be trusted, and how big a difference is meaningful (noise floor) |
| **1. SCAM PPE** | Thousands of nudged SCAM runs: many parameter sets × many columns × many monthly segments | Cheap (single column) | Emulators fitted; history matching removes implausible regions; NROY space |
| **2. 3D short-range** | Tens of candidate parameter sets in nudged 3D runs and CAPT hindcasts | Moderate | Whether SCAM improvements survive with circulation feedback |
| **3. Climate** | A few candidates in multi-year AMIP runs, with radiation-balance retuning | Expensive | Final parameter choice |

Each tier passes a smaller, better-constrained set of candidates to the next.

---

## 5. What exists today

### 5.1 Model configuration

| Item | Setting |
|---|---|
| Model | CAM6 physics, `compassdev` branch |
| Compset | `F2000climo_SROF`: F2000climo with a stub runoff model, because the standard runoff model fails in SCAM. Defined in `cime_config/config_compsets.xml`. |
| Grid | `ne30_ne30_mg17`. Physics must be on the np4 grid (not pg3), because the ERA5 nudging files are on ne30np4. |
| Build | Debug (`DEBUG = TRUE`), Intel compiler, Derecho |
| Nudging data | Hourly ERA5 on ne30np4, 32 levels: `/glade/campaign/cesm/stormspeed/ERA5_nudging_input_1h/2002-2008/ne30np4/L32/` |
| Nudging settings | U, V, T, Q coefficients 0.1667 (τ ≈ 6 h); PS not nudged; window 2002-01-01 to 2008-12-10 |
| 3D CAMIOP testmod | `cime_config/testdefs/testmods_dirs/cam/scm_prep_ERAnudge_global`, a copy of `scm_prep_todd_ERAnudge` with the `Nudge_hwin_*` lines removed. The original nudges everywhere *except* a box around Manus; without the window settings, nudging is global. |
| SCAM BFB settings | `scm_cambfb_mode = .true.`, `scm_use_3dfrc = .true.`, `scm_relaxation = .false.`, `scm_use_obs_uv = .true.` |
| Column | Manus. In NUOPC builds the column location is set by `scol_lat`/`scol_lon` in `nuopc.runconfig`, currently −1.8378, 147.8292. These are the ne30 column center nearest the site, not the ARM site coordinates. |
| History output | Tape 2 (`h1i` files): `TDIFF`, `QDIFF` (SCAM minus the 3D column, used for BFB checks). Tape 3 (`h2i` files): `Nudge_*`, `Target_*`, `num_a*`, cloud and surface fields, every step, 120 steps per file (`mfilt = 1,120,120,1,1,1`). |

**Key locations**

- **Code checkout (`$CAMROOT`):** `/glade/campaign/cgd/amp/jet/collections/cam6_4_207_compassdev.102726`
- **Base SCAM case:** `/glade/campaign/cgd/amp/jet/cases/CAM.F2000climo_SROF.ne30_ne30_mg17.intel.CAMIOP.SCAM.101a`. Its `run/` directory and `bld/cesm.exe` are the template and executable for every PPE member.
- **PPE areas (scratch):** `/glade/derecho/scratch/jet/scam_ppe/manus_5d` (test ensemble) and `/glade/derecho/scratch/jet/scam_ppe/manus_spin` (spin-up run)

### 5.2 What has been verified

| Check | Result |
|---|---|
| CICE single-column fix (degrees vs. radians for TLAT/TLON in `ice_mesh_mod`) | **Done.** Without it, F2000climo SCAM is not BFB with 3D. |
| 3D nudged CAMIOP run, 5 days from 2002-09-01, ne30, debug | **Done** |
| Nudged SCAM at Manus is BFB with the 3D column | **Verified.** `TDIFF = QDIFF = 0` at every step; `Nudge_T` and `Nudge_Q` nonzero, so nudging is actually on. |
| SCAM reads global ne30 nudging files correctly | **Verified** by the BFB result. In single-column mode, `src/control/ncdio_atm.F90` reads one column of a global file through a special `physgrid_scm` decomposition. |
| A PPE member run outside CIME (copied namelists, shared executable) is BFB with the CIME-run base case | **Verified.** `m0000` matched with `TDIFF = 0` and an identical `Nudge_T`. |
| Parameter changes take effect | **Verified.** `m0001` diverged by 0.47 K and 2.2e-4 kg/kg in 5 steps. |
| Several members in one PBS job | **Verified** with two |

### 5.3 The workflow scripts (`tools/ppe/`)

**Design choice: no CIME case per member.** A CIME clone per member would mean thousands of case directories and per-member namelist generation. Instead, one base SCAM case is built and run once with CIME. Each member is a plain directory holding copies of the base run directory's input files, with only two files changed:

- `atm_in`: the parameter values
- `nuopc.runconfig`: output paths and run length

All members run the same `cesm.exe`.

**A member directory** (e.g. `<ppe_root>/m0001/`) contains:

- **Inputs copied from the base run directory:** `atm_in`, `drv_in`, `drv_flds_in`, `lnd_in`, `ice_in`, `docn_in`, `docn.streams.xml`, `fd.yaml`, `nuopc.runconfig`, `nuopc.runseq`
- **Generated:** `ppe_member.json` (the member's parameters), `timing/`, `init_generated_files/`
- **Outputs:** CAM, CLM and CICE history files, logs, `cesm.log`, `time.txt` (memory and wall time from `/usr/bin/time -v`), and a `DONE` or `FAILED` marker

**The config file** (`ppe_config.toml`, copied from `ppe_config.example.toml`) has these sections:

| Section / key | Meaning |
|---|---|
| `case`, `base_rundir`, `exe`, `env_file` | Base SCAM case, its run directory, the executable, and `.env_mach_specific.sh` (the case's module environment) |
| `ppe_root` | Where member directories, `members.csv`, job scripts and metrics go |
| `[sampling] n_members`, `seed` | Size and random seed of the Latin hypercube |
| `[run] account, queue, job_priority, walltime, max_parallel` | PBS settings. `job_priority = "premium"` for tests. `max_parallel` is how many members run at once on one node (24, set by memory). |
| `[run] ncpus, mem` (optional) | Request part of a node, e.g. for small tests on the `develop` queue |
| `[run] stop` (optional) | Override the run length, e.g. `{ stop_option = "ndays", stop_n = 5 }` |
| `[params.<name>] range, scale` | Sampling range and `linear`/`log` scale for each perturbed namelist variable |
| `[analysis] spinup_hours, p_range_hpa` | Hours dropped at the start, and the pressure band for scalar metrics |

**The scripts**

| Script | Purpose |
|---|---|
| `ppe_sample.py CONFIG [--force]` | Reads each parameter's default from the base `atm_in` and writes `members.csv`: `m0000` = defaults, `m0001…` = Latin hypercube draws. Refuses to overwrite without `--force`. |
| `ppe_setup.py CONFIG [--members LIST] [--force]` | Creates member directories: copies the inputs, rewrites the 10 output-path (`diro`) lines in `nuopc.runconfig`, applies `[run] stop` if set, and patches parameter values into `atm_in`. `m0000`'s `atm_in` is byte-identical to the base. Only changes variables that already exist in `atm_in`; matching is exact, so `clubb_gamma_coef` never touches `clubb_gamma_coefb`. |
| `ppe_run.py CONFIG [--members LIST] --name NAME` | Writes, but does not submit, a PBS script `<ppe_root>/NAME.pbs`. The script runs members `max_parallel` at a time on one node, each with `mpiexec -n 1` in its own directory. The job log is `<ppe_root>/NAME.pbs.log` and ends with `done: N failed: M`. |
| `ppe_analyze.py CONFIG` | For every member with `DONE`: computes the nudging metrics (Section 7.5) and writes `<ppe_root>/ppe_metrics.nc` plus a summary table relative to `m0000`. **Tested on mock data only so far.** |
| `ppe_spinup.py MEMBER_DIR [--p-range 100 1000]` | Prints 6-hourly averages and plots (`spinup.png`) the nudging tendencies and aerosol numbers against time since start. Used to choose the spin-up discard (Section 6). |

**Parameters currently sampled (11)**

- **ZM:** `zmconv_c0_ocn` (rain conversion efficiency over ocean), `zmconv_ke` (rain evaporation efficiency), `zmconv_tau` (CAPE consumption timescale), `zmconv_dmpdz` (entrainment rate), `zmconv_capelmt` (CAPE threshold)
- **CLUBB:** `clubb_gamma_coef`, `clubb_c14`, `clubb_c11`, `clubb_c11b`, `clubb_beta`
- **MG:** `micro_mg_autocon_fact` (autoconversion)

The ranges were checked against Eidhammer et al. (2024), Table 1. Three (`zmconv_c0_ocn`, `zmconv_dmpdz` and `zmconv_capelmt`) were changed to the table's values; the others already agreed with it. Land-specific parameters (`zmconv_c0_lnd`, `zmconv_ke_lnd`) should be added when land columns are included.

---

## 6. Tier 0: verifying the foundation

Each check below has to pass before results from the next tier can be trusted.

| # | Check | Why | Method | Pass criterion | Status |
|---|---|---|---|---|---|
| 0.1 | Nudged BFB | The default member must reproduce 3D exactly | SCAM vs. 3D column | `TDIFF = QDIFF = 0`, `Nudge_T ≠ 0` | **Done** |
| 0.2 | Member mechanism | Running outside CIME must change nothing | `m0000` vs. base case | Identical output | **Done** |
| 0.3 | Full-day run and analysis on real files | Memory growth, wall time, first real use of `ppe_analyze.py` | 1-day run of `m0000`, `m0001` | Memory flat; analysis table produced | **In progress** (job `test2`) |
| 0.4 | Spin-up length | How many days to discard at the start of each segment | 5-day default run; later two 3D segments started a week apart | Nudging and aerosol series stop trending | **In progress** (job `spin5`) |
| 0.5 | Noise floor | Smallest meaningful difference between two members | Default-parameter members with tiny initial perturbations (below) | Spread of each metric | **Planned** |
| 0.6 | Failure handling | Extreme parameters may crash, hang or produce NaNs | 3–4 members at corners of the ranges | Each case recorded correctly | **Planned** |
| 0.7 | Scale | Many members reading the same files at once | 24 members on one node | Per-member time close to the single-run time | **Planned** |
| 0.8 | Other columns | Same base case, different location | Patch `scol_lat`/`scol_lon` | `TDIFF = 0` against the 3D column at the new site | **Planned** |
| 0.9 | Debug vs. optimized | Production must use optimized builds | Default member, both builds | Metric difference within the noise floor | **Planned** |

**Noise floor method (0.5).** In this nonlinear system a perturbation at roundoff level grows until it saturates at the level of the system's internal variability. Differences between PPE members are only meaningful if they exceed that spread. Run a few default-parameter members whose initial temperature is perturbed at the 1e-14 level, then compute the spread of each metric across them.

CAM's namelist variable `pertlim` (group `cam_initfiles_nl`) multiplies the initial temperature by `1 + r`, with `r` uniform in `±pertlim`. In SE SCAM it's applied while the initial file is read (`read_inidat`, `src/dynamics/se/dyn_comp.F90`). In CAMIOP mode T is prognostic and is never reset from the IOP file, so the perturbation carries through the run and `pertlim` should work for SCAM. That still needs one test: a member with `pertlim = 1e-14` should show small nonzero `TDIFF`.

Two practical points:

- The random pattern comes from Fortran's `random_number` with its default seed, so ensemble members should use different `pertlim` values (e.g. 1e-14, 2e-14, 3e-14…) to get distinct perturbations.
- `pertlim` isn't in the base `atm_in` by default, so add `pertlim = 0.0` to the base case's `user_nl_cam` once. That lets `ppe_setup.py` patch it per member.

**Failure handling (0.6) needs new script features:**

- a per-member timeout, because a hung member currently just runs into the job walltime and leaves no marker
- NaN detection in `ppe_analyze.py`
- recording crashes as their own category, because for history matching a crash means "implausible region", not "missing data"
- a `--failed-only` option to rerun failed members

**Changing the column (0.8).** The land domain file named in `nuopc.runconfig` (`single_column_lnd_domainfile`) is global, so `scol_lat`/`scol_lon` are the only site-specific settings. Set them to the exact ne30 column center nearest the site, read from the CAMIOP file's `lat`/`lon`. Then CAM, CLM, CICE and the data ocean all select the same point.

---

## 7. Tier 1: the SCAM PPE in detail

### 7.1 Producing the CAMIOP forcing: overlapping monthly segments

A SCAM run is BFB with its 3D partner only if it starts **exactly where the 3D run started (its first time step) and from the same initial files**. Four ways of covering 2002–2008 were considered.

**Option 1: one continuous 3D run, with SCAM started partway through it.** *Rejected.* In CAMIOP mode T, q and the tracers are prognostic, and their initial values come from the initial file named in `ncdata`, not from the CAMIOP file. The CAMIOP `t` and `q` are used only for the `TDIFF`/`QDIFF` comparison. A mid-file start would therefore begin from the initial file's state, not the 3D model's state on that day. The physics "memory" would also start cold: CLUBB's turbulence moments, radiative heating rates and the other physics-buffer fields that persist between steps.

**Option 2: warm start from 3D restart and initial files.** *Set aside for now.* SCAM would read its column of the physics buffer from a 3D restart file and start CAM's clock as a restart run. Scoped at about 150–250 lines of code. Even then it would only be nearly BFB, because the coupler's first-step surface fluxes would differ.

**Option 3: one continuous SCAM run that writes its own restart files, with ensemble segments restarted from them.** *Set aside.* Exact SCAM restart hasn't been tested, and a serial multi-year SCAM run isn't affordable at current speed (Section 10).

**Option 4: independent, overlapping monthly 3D segments.** *Chosen.* Each month gets its own short 3D run that starts a few spin-up days before the month and ends at the month's end. Each segment is an ordinary nudged startup run with its own CAMIOP file. SCAM runs for that segment start at the segment's first step, which is exactly the configuration verified to be BFB.

Example with a spin-up of `S` days:

```
3D segment "2002-09":  start 2002-09-01 minus S days  ->  end 2002-10-01
3D segment "2002-10":  start 2002-10-01 minus S days  ->  end 2002-11-01
...
Analysis window of each segment = its calendar month only.
Consecutive segments overlap by S days; analysis windows don't overlap
and together cover 2002-2008 without gaps.
```

Properties of Option 4:

- **Exact.** Each SCAM run reproduces its 3D segment bit-for-bit, with no new code.
- **Parallel.** All 84 monthly 3D segments can run at once, and so can all SCAM runs.
- **No SCAM restarts and no continuous multi-year file are needed.**
- **The cost is the spin-up days.** Every segment starts from the same climatological initial state, which the nudging pulls onto ERA5. Those first `S` days are run, both in 3D and in every SCAM member, and then discarded.
- **Segments differ from a continuous run, and that doesn't matter.** After spin-up the state is set by ERA5, not by the initial file, and the PPE metric is a time average over the analysis window.

**How long must `S` be?** The discard window must cover two adjustments, which happen at the same time, so it's the longer of the two, not their sum:

1. **Adjustment of the 3D segment from its initial state.** T, q, u and v lock onto ERA5 within about 1–2 days (τ ≈ 6 h); clouds and turbulence within hours. **Aerosols (MAM4) and land state are not nudged.** Aerosols may take a week or two; land soil moisture takes weeks to months.
2. **Adjustment of a perturbed SCAM member to its new parameters.** Fast physics, probably a day or two.

`S` must be measured, not guessed:

- **Column test (in progress, job `spin5`):** a 5-day default run analyzed with `ppe_spinup.py`. When the nudging tendencies stop trending and the `num_a*` aerosol ratios stop drifting, spin-up is over.
- **Decisive test (planned):** two 3D segments started about a week apart, for example 2002-08-25 and 2002-09-01, compared day by day. Once they agree to within noise, both have forgotten their initial state. That time is the required `S`. It automatically includes aerosols and any remote influence of land on the forcing.

**Land columns.** Every segment starts CLM from the same climatological initial file (`finidat`), so soil moisture can be wrong for weeks or more, biasing surface fluxes at land columns such as Darwin. Ocean columns (Manus, Nauru) are largely unaffected. Possible remedies for land columns, all still to be designed: longer spin-up, or segment-specific CLM initial files taken from an earlier continuous run.

**Building a 3D segment.** This is planned; no script exists yet. Each segment is a CIME case using the same compset, grid and nudging as the verified pair, with three changes:

- `RUN_STARTDATE` = month start minus `S`
- run length = month length plus `S`
- an optimized build for production

The same initial files (`ncdata`, `finidat`) must be used for every segment, and SCAM must use exactly the same ones. A small generator script that creates and submits all segments from a template should be written.

**Data volume.** CAMIOP output is large. It writes about 41 three-dimensional fields (including one tendency per tracer; CAM6/MAM4 has about 33 tracers) for every column, every time step, in double precision:

- about 500 MB per step at ne30
- about 24 GB per simulated day
- about 60 TB for 2002–2008

SCAM needs only the columns we choose. Options:

- **(a)** write daily files and extract the chosen columns with `ncks -d ncol,...`, deleting the global files as the run proceeds
- **(b)** single precision (`NDENS = 2`), acceptable once BFB isn't required
- **(c)** a code change so CAMIOP writes only selected columns; SCAM's IOP reader would also need to accept that format

Option (a) needs no code. It must first be checked that SCAM gives the same answer from a column-subset file as from the full file (BFB test).

### 7.2 Choosing columns

- **ARM Tropical Western Pacific sites (2002–2008):** Manus (deep convection over the warm pool), Nauru (on the edge of the warm pool, where convection is often suppressed), and Darwin (coastal land, monsoon and break regimes). Observations exist at these sites (Section 7.6).
- **Additional regimes without ARM data:** trade-wind cumulus, subtropical stratocumulus, and mid-latitude columns. The nudging-tendency metric doesn't need local observations, so these are nearly free once column extraction works. They guard against tuning that only fits the western Pacific.
- **Practical rule:** use each site's nearest ne30 column center for `scol_lat`/`scol_lon` (Section 6, check 0.8).

### 7.3 Ensemble design

- **First wave:** a Latin hypercube over the parameters, about 10 members per parameter (about 110 for 11 parameters), plus the default member. With 84 monthly segments and N columns, each member is run for every segment × column.
- **Later waves:** sample only inside the NROY space from history matching (Section 7.8), refit the emulators, repeat.
- **Cost:** per member and column, (30 + `S`) days per segment. Reducing the per-day cost (Section 10) is what makes this affordable.

### 7.4 Running members across segments and columns

The current scripts handle one base case, one segment and one column. To cover many segments and columns, `ppe_setup.py` needs a general per-member override mechanism, both for `atm_in` variables and for `nuopc.runconfig` entries:

| Variable | Where | Purpose |
|---|---|---|
| `iopfile` | `atm_in` | The segment's CAMIOP file |
| `ncdata` | `atm_in` | Initial file (must match the 3D segment) |
| `pertlim` | `atm_in` | Noise-floor members |
| `start_ymd`, `start_tod` | `nuopc.runconfig` | Segment start, equal to the 3D segment's start |
| `stop_option`, `stop_n` | `nuopc.runconfig` | Segment length (already supported) |
| `scol_lat`, `scol_lon` | `nuopc.runconfig` | Column location |

A natural layout is `<ppe_root>/<column>/<segment>/<member>/`, with `ppe_analyze.py` extended to collect over segments and columns.

### 7.5 Metrics from nudging tendencies

For each member, column and segment, `ppe_analyze.py` uses `Nudge_T` (converted to K/day) and `Nudge_Q` (converted to g/kg/day). It drops the first `spinup_hours` and computes, for each model level k:

- **Time-mean profile** `mean_k = ⟨Nudge⟩_t`: the systematic correction, which is the Rodwell–Palmer-style error signal.
- **RMS profile** `rms_k = sqrt(⟨Nudge²⟩_t)`: systematic plus random correction.

These are combined into mass-weighted scalars over a pressure band (default 100–1000 hPa). The layer weights `w_k` are proportional to the pressure thickness Δp_k, computed from the hybrid coefficients (`hyai`, `hybi`) and mean surface pressure:

- **`bias`** = sqrt( Σ w_k · mean_k² / Σ w_k ): the size of the systematic error
- **`tot`** = Σ w_k · rms_k / Σ w_k: total correction

Each scalar is also given **relative to `m0000`**: a value of 0.8 means 20% less nudging than the default. The output file `ppe_metrics.nc` stores the profiles (member × level), the scalars, the relative scalars and the parameter values.

**Still to do:**

- aggregate over segments and columns
- carry the noise floor as a per-metric uncertainty
- consider separating levels or regimes (e.g. lower and upper troposphere) as separate metrics, so improvement in one layer can't hide degradation in another
- consider using the leading principal components of the mean profile as emulator targets

### 7.6 Independent observations

Nudging metrics alone can reward compensating errors, so each candidate should also be scored against data the nudging doesn't use.

| Data | Variables | Notes |
|---|---|---|
| ARM Manus, Nauru, Darwin (2002–2008) | Surface precipitation; broadband surface radiative fluxes; cloud-fraction profiles from ARSCL (radar–lidar); liquid water path and precipitable water from microwave radiometers; MICROBASE cloud microphysics retrievals | Point measurements versus a ~100 km ne30 column: compare time averages and statistics (diurnal cycle, distributions), not individual times. |
| Satellite precipitation (e.g. GPM IMERG) | Precipitation over the grid box | Matches the column's spatial scale better |
| CERES (e.g. EBAF or SYN products) | Top-of-atmosphere and surface radiative fluxes, cloud radiative effect | Key protection against radiation–convection compensation |

**Choosing comparable quantities.** Nudged T and q are close to ERA5 by construction, so they aren't useful independent metrics. Precipitation, cloud, liquid water path and radiation are the informative quantities. Observational uncertainty and representativeness error go into the `Var_obs` term of the implausibility (Section 2.7).

### 7.7 Emulators

**Planned.** For each scalar metric (per column or regime, aggregated over segments), fit a Gaussian-process emulator of metric versus parameters. Python options include scikit-learn's `GaussianProcessRegressor` or GPy.

- **Inputs:** parameters scaled to [0, 1], in log space for log-sampled parameters.
- **Validation:** leave-one-out or held-out members. The emulator's predicted uncertainty should cover the true values at roughly the stated rate.
- **Profiles:** reduce with principal components and emulate the leading coefficients.

### 7.8 History matching

**Planned.** For each candidate θ (millions of emulator evaluations are cheap) and each metric, compute the implausibility `I(θ)` from Section 2.7:

- **For nudging metrics**, the target is zero systematic correction. `Var_discrepancy` must include the part that can't be removed, which is estimated from ERA5 errors and the best members of the ensemble, and `Var_obs` is replaced by the noise floor.
- **For observational metrics**, `z` is the observed value and `Var_obs` includes measurement and representativeness uncertainty.
- **Rule out** θ where the maximum (or second-highest) implausibility across metrics exceeds 3.
- **Waves:** sample new members inside NROY, run them in SCAM, refit, repeat until NROY stops shrinking or the emulator uncertainty is small compared with the other terms.
- **Output:** a set of plausible parameter vectors, not one optimum. A representative handful goes to Tier 2.

---

## 8. Tier 2: 3D short-range validation

**Planned.** This tier tests the candidates with circulation feedback included.

- **Nudged 3D runs** with each candidate: compare the global, 3D nudging tendencies with what SCAM predicted. This directly tests whether the SCAM improvement holds in 3D and shows where it fails.
- **CAPT-style hindcasts** (Phillips et al. 2004): from the nudged default run, save daily initial files, then run many free 3–5 day hindcasts per candidate. Score them against ERA5 (temperature, humidity and wind errors by lead time) and against satellite precipitation and radiation. Check whether day-1 to day-5 errors look like the known climate biases (Xie et al. 2012; Ma et al. 2014). For full fidelity, a candidate's hindcasts should start from a nudged spin-up with its own parameters, so clouds and soil moisture match its physics.

---

## 9. Tier 3: climate validation

**Planned.** Run a few candidates in multi-year AMIP simulations. Short-range metrics don't constrain the top-of-atmosphere energy balance or slow feedbacks, so each candidate needs a radiation-balance retune, usually with a small set of cloud and microphysics parameters (see Hourdin et al. 2017), and evaluation against standard climate diagnostics.

---

## 10. Current status and bottlenecks

**Jobs in progress**

- `test2` (`manus_5d`): `m0000` and `m0001`, 1 day. Measures memory over a full day and wall time, and is the first real run of `ppe_analyze.py`.
- `spin5` (`manus_spin`): `m0000`, 5 days. Spin-up measurement with `ppe_spinup.py`.

**Bottlenecks**

1. **Speed.** About 30 s per time step for one column (1430 s per simulated day, 98% in the atmosphere), plus about 210 s of initialization. That's about 2 h per member for 5 days, or about 40 days of wall time for 7 years. One column of CAM6 physics should take well under a second per step, so most of the time is overhead. The suspects:
   - reading the global CAMIOP file every step
   - single-column reads of the global ERA5 files
   - the debug build
   - possibly SE SCAM carrying the whole ne30 element structure (`scm_dyn_grid_indicies` searches the full `elem` array for the column)

   **First step:** profile with the timers written to the case's `timing/` directory.
2. **Memory.** 8.5 GB per member, probably from setting up the full ne30 grid (spectral-element grid, coupler mesh, I/O decompositions). Only about 24 members fit on a 235 GB Derecho node, and node-hours are charged whole, so each single-core member effectively costs about 5 cores.
3. **CAMIOP volume.** About 60 TB for the full period if written globally (Section 7.1).
4. **Spin-up length.** Not yet known (Section 7.1). The first-step nudging tendencies are large (about 48 K/day in T and 17 g/kg/day in q) because the initial state is climatological.

---

## 11. Next steps, in order

1. **Finish Tier 0 checks 0.3 and 0.4** (`test2`, `spin5`), and fix `ppe_analyze.py` if real CAM files expose problems (time decoding, the shapes of `hyai`/`P0`).
2. **Two-segment 3D spin-up test** to set `S` (Section 7.1).
3. **Noise floor** with `pertlim` (check 0.5); confirm that `pertlim` acts in CAMIOP-mode SCAM.
4. **Failure handling** (check 0.6) and the related script features.
5. **Scale test** with 24 members (check 0.7).
6. **Second site, Nauru,** by patching `scol_lat`/`scol_lon` (check 0.8); then Darwin, watching `LANDFRAC`.
7. **Per-member overrides in `ppe_setup.py`** (Section 7.4) and multi-segment, multi-column analysis.
8. **Performance:** profile, optimized build (with check 0.9), column-subset CAMIOP files, the memory footprint.
9. **Segment generator** for 3D CAMIOP production, and column extraction; then production for 2002–2008.
10. **Observational metrics** (ARM, IMERG, CERES), emulators and history matching (Sections 7.6–7.8).
11. **Tiers 2 and 3.**

---

## Appendix A: Shell environment and variables

The commands are written for **tcsh**, the login shell on Derecho, with Python from the `npl-2025b` conda environment. Variables are set with `set NAME = value` and used as `$NAME`. They are typed per session and aren't stored anywhere.

| Variable | Value (example) | Meaning |
|---|---|---|
| `CAMROOT` | `/glade/campaign/cgd/amp/jet/collections/cam6_4_207_compassdev.102726` | Code checkout (`compassdev`) |
| `PP` | `$CAMROOT/tools/ppe` | The PPE scripts |
| `SC` | `/glade/campaign/cgd/amp/jet/cases/CAM.F2000climo_SROF.ne30_ne30_mg17.intel.CAMIOP.SCAM.101a` | Base SCAM case directory |
| `B` (or `RD`) | `$SC/run` | Base case run directory, the template for members |
| `PR` | `/glade/derecho/scratch/$USER/scam_ppe/manus_5d` | PPE root of the test ensemble |
| `C` | `$PR/ppe_config.toml` | That ensemble's config |
| `P2` | `/glade/derecho/scratch/$USER/scam_ppe/manus_spin` | PPE root of the spin-up run |
| `C2` | `$P2/ppe_config.toml` | Its config |
| `TM` | `$CAMROOT/cime_config/testdefs/testmods_dirs/cam` | CAM testmods |

**A typical session**

```
set CAMROOT = /glade/campaign/cgd/amp/jet/collections/cam6_4_207_compassdev.102726
set PP = $CAMROOT/tools/ppe
set PR = /glade/derecho/scratch/$USER/scam_ppe/manus_5d
set C = $PR/ppe_config.toml

python3 $PP/ppe_sample.py $C                          # members.csv
python3 $PP/ppe_setup.py  $C --members m0000,m0001    # member directories
python3 $PP/ppe_run.py    $C --members m0000,m0001 --name test2
qsub $PR/test2.pbs                                    # you submit; scripts never do
tail -1 $PR/test2.pbs.log                             # "done: 2 failed: 0"
python3 $PP/ppe_analyze.py $C                         # ppe_metrics.nc + table
python3 $PP/ppe_spinup.py $PR/m0000                   # spin-up table + plot
```

**Starting a new PPE area from an existing config:** copy the config, edit `ppe_root` and anything else, then run the steps above. In tcsh, use single quotes for TOML strings inside `sed`, e.g. `sed -i "s#^ppe_root = .*#ppe_root = '$P2'#" $C2`.

## Appendix B: Pitfalls already encountered

- **Nudging can be silently off.** If the run starts before the nudging window (`Nudge_Beg_*`), nudging waits and never acts in a short run. A first "BFB" test started on 0001-09-01 and was BFB only because nudging was off in both runs. **Always confirm `Nudge_T` is nonzero.**
- **History file sizes.** With `mfilt = 1` for a tape, CAM writes one file per output time. Members copy the base `atm_in`, so set history options in the base case before creating members.
- **tcsh quoting.** `\"` inside double quotes doesn't work. With `<< 'EOF'`, the terminator line must literally be `'EOF'`, quotes included. Use unquoted `EOF` with bodies that contain no `$`, `!` or backticks.
- **Long single-line commands** overflow the terminal; split them with intermediate variables.
- **Job priority** can't be changed on a queued job: `qdel` and resubmit. Premium priority starts small jobs in about a minute.
- **No builds on login nodes:** run `create_test` and `case.build` through `qcmd`.
- **Concurrent jobs must not share member directories.** Two jobs writing the same `mNNNN/` would corrupt each other's output; use separate `ppe_root`s.

## Appendix C: Code references (`compassdev`)

| File | What to look at |
|---|---|
| `src/dynamics/se/apply_iop_forcing.F90` | `advance_iop_forcing`: fixed-tendency forcing; winds set from the IOP file |
| `src/dynamics/se/advect_tend.F90`, `stepon.F90` | How CAMIOP tendencies are computed (dynamics residual, within each step) |
| `src/dynamics/se/se_single_column_mod.F90` | `scm_setinitial`: in CAMIOP mode the state isn't set from the IOP file (T and q are prognostic; IOP `t`/`q` feed only `TDIFF`/`QDIFF`); `scm_dyn_grid_indicies` |
| `src/dynamics/se/dyn_comp.F90` | `read_inidat` (initial state from `ncdata`; `pertlim`) |
| `src/control/scamMod.F90` | `setiopupdate_init` (IOP time index chosen by date); fields read from the IOP file |
| `src/physics/cam/nudging.F90` | Nudging window (`nudging_init`); `Nudge_*`/`Target_*` diagnostics |
| `src/control/ncdio_atm.F90` | `physgrid_scm` single-column reads of global files |
| `src/cpl/nuopc/atm_comp_nuopc.F90` | Column location from `scol_lat`/`scol_lon` |
| `src/control/history_scam.F90` | Which fields CAMIOP writes, and on which tape |
| `cime_config/SystemTests/sct.py` | The `SCT` system test (3D → CAMIOP → SCAM BFB) |

## Appendix D: References

*Citations are from memory; check the details (year, journal, pages) before using them in a publication.*

- Couvreux, F., and coauthors, 2021: Process-based climate model development harnessing machine learning: I. A calibration tool for parameterization improvement. *J. Adv. Model. Earth Syst.*, 13, e2020MS002217.
- Eidhammer, T., and coauthors, 2024: An extensible perturbed parameter ensemble for the Community Atmosphere Model version 6. *Geosci. Model Dev.*, 17, 7835–7853. (Parameter ranges: Table 1.)
- Gettelman, A., J. E. Truesdale, J. T. Bacmeister, P. M. Caldwell, R. B. Neale, P. A. Bogenschutz, and I. R. Simpson, 2019: The Single Column Atmosphere Model version 6 (SCAM6): Not a scam but a tool for model evaluation and development. *J. Adv. Model. Earth Syst.*, 11, 1381–1401.
- Hourdin, F., and coauthors, 2017: The art and science of climate model tuning. *Bull. Amer. Meteor. Soc.*, 98, 589–602.
- Hourdin, F., and coauthors, 2021: Process-based climate model development harnessing machine learning: II. Model calibration from single column to global. *J. Adv. Model. Earth Syst.*, 13, e2020MS002225.
- Klinker, E., and P. D. Sardeshmukh, 1992: The diagnosis of mechanical dissipation in the atmosphere from large-scale balance requirements. *J. Atmos. Sci.*, 49, 608–627.
- Kooperman, G. J., M. S. Pritchard, S. J. Ghan, M. Wang, R. C. J. Somerville, and L. M. Russell, 2012: Constraining the influence of natural variability to improve estimates of global aerosol indirect effects in a nudged version of the Community Atmosphere Model 5. *J. Geophys. Res.*, 117, D23204.
- Loeppky, J. L., J. Sacks, and W. J. Welch, 2009: Choosing the sample size of a computer experiment: A practical guide. *Technometrics*, 51, 366–376.
- Ma, H.-Y., and coauthors, 2014: On the correspondence between mean forecast errors and climate errors in CMIP5 models. *J. Climate*, 27, 1781–1798.
- Phillips, T. J., and coauthors, 2004: Evaluating parameterizations in general circulation models: Climate simulation meets weather prediction. *Bull. Amer. Meteor. Soc.*, 85, 1903–1915.
- Randall, D. A., K.-M. Xu, R. J. C. Somerville, and S. Iacobellis, 1996: Single-column models and cloud ensemble models as links between observations and climate models. *J. Climate*, 9, 1683–1697.
- Rodwell, M. J., and T. N. Palmer, 2007: Using numerical weather prediction to assess climate models. *Quart. J. Roy. Meteor. Soc.*, 133, 129–146.
- Williams, K. D., and coauthors, 2013: The Transpose-AMIP II experiment and its application to the understanding of Southern Ocean cloud biases in climate models. *J. Climate*, 26, 3258–3274.
- Williamson, D., M. Goldstein, L. Allison, A. Blaker, P. Challenor, L. Jackson, and K. Yamazaki, 2013: History matching for exploring and reducing climate model parameter space using observations and a large perturbed physics ensemble. *Climate Dyn.*, 41, 1703–1729.
- Xie, S., H.-Y. Ma, J. S. Boyle, S. A. Klein, and Y. Zhang, 2012: On the correspondence between short- and long-time-scale systematic errors in CAM4/CAM5 for the Year of Tropical Convection. *J. Climate*, 25, 7937–7955.
- Zhang, K., H. Wan, X. Liu, S. J. Ghan, G. J. Kooperman, P.-L. Ma, P. J. Rasch, D. Neubauer, and U. Lohmann, 2014: Technical Note: On the use of nudging for aerosol–climate model intercomparison studies. *Atmos. Chem. Phys.*, 14, 8631–8645.
