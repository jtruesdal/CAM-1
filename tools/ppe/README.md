# SCAM PPE tools

See `PPE_GUIDE.md` for the science, the full workflow and next steps
(draft, under review).

Perturbed-parameter ensemble (PPE) workflow for nudged SCAM driven by a
CAMIOP file. One base SCAM case is built and run once with CIME; every
member is a plain directory with copied namelists, run with the shared
`cesm.exe`. No CIME clone per member.

1. `cp ppe_config.example.toml <ppe_root>/ppe_config.toml` and edit paths,
   `[run] account`, and the parameter ranges.
2. `ppe_sample.py cfg` writes `<ppe_root>/members.csv`: `m0000` holds the
   base-case values from the base `atm_in`, followed by Latin hypercube draws.
3. `ppe_setup.py cfg` creates `<ppe_root>/mNNNN/`, copying the base run
   files, rewriting the `diro` paths in `nuopc.runconfig`, optionally setting
   the run length, and patching `atm_in`. `m0000` gets an unmodified `atm_in`
   and is the BFB check against the base case.
4. `ppe_run.py cfg` writes a PBS script that runs up to `max_parallel`
   members at once on one node. It does not submit; run `qsub` yourself.
   Each member leaves `cesm.log`, `time.txt` (max RSS), and `DONE` or
   `FAILED`. `[run] machine` selects the batch header and MPI launcher
   (`derecho` or `izumi`).
5. `ppe_analyze.py cfg` reads `Nudge_T`/`Nudge_Q` from every member with a
   `DONE` marker, drops `[analysis] spinup_hours`, and writes
   `<ppe_root>/ppe_metrics.nc`: time-mean and RMS profiles per member plus
   mass-weighted scalars over `p_range_hpa` (`bias_*` from the mean profile,
   `tot_*` from the RMS profile), each also relative to `m0000`, alongside
   the parameter values. Needs all of a member's steps in its history files
   (set a large `mfilt` for the tape holding the nudging fields).

Only variables already present in the base `atm_in` can be perturbed. Add
any others to the base case's `user_nl_cam` and rerun `preview_namelists`
(or rerun the case) first.
