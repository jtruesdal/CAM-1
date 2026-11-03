#!/usr/bin/env python3
"""Create one run directory per PPE member from members.csv.

For each member: copy the base RUNDIR namelists/config files, point the
nuopc.runconfig diro lines at the member directory, optionally override the
run length, and patch the member's parameter values into atm_in. m0000 gets
unmodified copies (apart from the diro paths) and is the BFB reference.

usage: ppe_setup.py ppe_config.toml [--members m0000,m0003] [--force]
"""
import argparse
import json
import shutil
import sys

from ppe_common import (RUN_DIRS, RUN_FILES, load_config, patch_nml,
                        patch_runconfig, read_members)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config")
    ap.add_argument("--members", help="comma-separated subset to set up")
    ap.add_argument("--force", action="store_true",
                    help="recreate member directories that already exist")
    args = ap.parse_args()
    cfg = load_config(args.config)
    base, root = cfg["base_rundir"], cfg["ppe_root"]

    missing = [f for f in RUN_FILES if not (base / f).is_file()]
    if missing:
        sys.exit(f"base rundir {base} is missing: {', '.join(missing)}")

    names, members = read_members(root / "members.csv")
    if args.members:
        want = set(args.members.split(","))
        members = [m for m in members if m[0] in want]
        if len(members) != len(want):
            sys.exit("some requested members are not in members.csv")

    stop = cfg.get("run", {}).get("stop")
    base_atm = (base / "atm_in").read_text()
    base_rc = (base / "nuopc.runconfig").read_text()

    for name, vals in members:
        d = root / name
        if d.exists():
            if not args.force:
                print(f"{name}: exists, skipped (use --force)")
                continue
            shutil.rmtree(d)
        d.mkdir(parents=True)
        for sub in RUN_DIRS:
            (d / sub).mkdir(parents=True, exist_ok=True)
        for f in RUN_FILES:
            if f not in ("atm_in", "nuopc.runconfig"):
                shutil.copy2(base / f, d / f)

        rc, n = patch_runconfig(base_rc, base, d, stop)
        if n == 0:
            sys.exit("nuopc.runconfig has no base RUNDIR paths; check config")
        (d / "nuopc.runconfig").write_text(rc)

        atm = base_atm if name == "m0000" else patch_nml(base_atm, vals)
        (d / "atm_in").write_text(atm)

        (d / "ppe_member.json").write_text(json.dumps(
            {"member": name, "params": vals, "base_rundir": str(base),
             "stop": stop}, indent=1) + "\n")
        print(f"{name}: created ({n} rundir paths rewritten)")


if __name__ == "__main__":
    main()
