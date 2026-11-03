#!/usr/bin/env python3
"""Write a PBS job script that runs PPE members from their own directories.

Each member is one single-task cesm.exe run in its member directory. The job
runs up to run.max_parallel members at once on one node. The script is only
written, not submitted: submit it yourself with qsub.

[run] machine selects the batch header and MPI launcher:
  derecho (default)  PBS Pro: -A account, select=..., job_priority;
                     launcher "mpiexec -n 1 --cpu-bind none" (Cray PALS)
  izumi              PBS/Torque: -l nodes=1:ppn=N, no account or priority;
                     launcher "mpiexec -n 1" with CPU binding disabled for
                     mvapich2 and openmpi so concurrent members don't all
                     land on core 0

usage: ppe_run.py ppe_config.toml [--members m0000,m0001] [--name NAME]
"""
import argparse
import sys

from ppe_common import load_config, read_members

HEADERS = {
    "derecho": """#PBS -N {jobname}
#PBS -A {account}
#PBS -q {queue}
#PBS -l select=1:ncpus={ncpus}:mpiprocs={ncpus}:ompthreads=1{mem}
#PBS -l walltime={walltime}
#PBS -l job_priority={job_priority}
#PBS -j oe
#PBS -o {root}/{jobname}.pbs.log
""",
    "izumi": """#PBS -N {jobname}
#PBS -q {queue}
#PBS -l nodes=1:ppn={ncpus}
#PBS -l walltime={walltime}
#PBS -S /bin/bash
#PBS -j oe
#PBS -o {root}/{jobname}.pbs.log
""",
}

LAUNCH = {
    "derecho": "mpiexec -n 1 --cpu-bind none",
    "izumi": "mpiexec -n 1",
}

DEFAULTS = {  # queue, ncpus, max_parallel
    "derecho": ("main", 128, 24),
    "izumi": ("short", 48, 8),
}

BODY = """
source {env_file}
export OMP_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0
export OMPI_MCA_hwloc_base_binding_policy=none
EXE={exe}
MAXP={max_parallel}

run_member() {{
    cd "$1" || return 1
    rm -f DONE FAILED
    if {launch} /usr/bin/time -v -o time.txt \\
            "$EXE" > cesm.log 2>&1; then
        touch DONE
    else
        touch FAILED
    fi
}}

MEMBERS=(
{member_lines}
)

for d in "${{MEMBERS[@]}}"; do
    while [ "$(jobs -rp | wc -l)" -ge "$MAXP" ]; do wait -n; done
    run_member "$d" &
done
wait
echo "done: $(ls {root}/m*/DONE 2>/dev/null | wc -l)" \\
     "failed: $(ls {root}/m*/FAILED 2>/dev/null | wc -l)"
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config")
    ap.add_argument("--members", help="comma-separated subset to run")
    ap.add_argument("--name", default="ppe_run", help="job and file name")
    args = ap.parse_args()
    cfg = load_config(args.config)
    root, run = cfg["ppe_root"], cfg.get("run", {})

    machine = run.get("machine", "derecho")
    if machine not in HEADERS:
        sys.exit(f"config: [run] machine '{machine}' must be one of "
                 f"{', '.join(HEADERS)}")
    dq, dn, dp = DEFAULTS[machine]

    _, members = read_members(root / "members.csv")
    names = [m[0] for m in members]
    if args.members:
        want = args.members.split(",")
        bad = [m for m in want if m not in names]
        if bad:
            sys.exit(f"not in members.csv: {', '.join(bad)}")
        names = want
    unset = [m for m in names if not (root / m / "atm_in").is_file()]
    if unset:
        sys.exit(f"not set up yet (run ppe_setup.py): {', '.join(unset)}")

    required = ("account", "walltime") if machine == "derecho" \
        else ("walltime",)
    for key in required:
        if key not in run:
            sys.exit(f"config: [run] needs '{key}' on {machine}")
    priority = run.get("job_priority", "regular")
    if priority not in ("economy", "regular", "premium"):
        sys.exit(f"config: [run] job_priority '{priority}' must be "
                 "economy, regular or premium")
    if machine == "izumi" and "job_priority" in run:
        print("note: job_priority is ignored on izumi")
    mem = f":mem={run['mem']}" if "mem" in run and machine == "derecho" \
        else ""

    fields = dict(
        jobname=args.name, account=run.get("account", ""),
        queue=run.get("queue", dq), walltime=run["walltime"],
        job_priority=priority, ncpus=int(run.get("ncpus", dn)), mem=mem,
        root=root, env_file=cfg["env_file"], exe=cfg["exe"],
        max_parallel=int(run.get("max_parallel", dp)),
        launch=LAUNCH[machine],
        member_lines="\n".join(f'    "{root / m}"' for m in names))
    script = "#!/bin/bash\n" + HEADERS[machine].format(**fields) + \
        BODY.format(**fields)
    out = root / f"{args.name}.pbs"
    out.write_text(script)
    print(f"wrote {out} ({len(names)} members, {machine}); "
          f"submit with: qsub {out}")


if __name__ == "__main__":
    main()
