#!/usr/bin/env python3
"""Draw a Latin hypercube over the PPE parameters and write members.csv.

Member m0000 holds the base-case defaults read from the base atm_in; it is
run from an unmodified copy of the namelists and is the BFB reference.
Members m0001..mNNNN are the LHS draws.

usage: ppe_sample.py ppe_config.toml
"""
import argparse
import sys

import numpy as np

from ppe_common import (fortran_to_float, load_config, read_nml_value,
                        write_members)


def lhs(n, d, rng):
    """Plain Latin hypercube on [0,1)^d: one point per stratum per dim."""
    u = (rng.random((n, d)) + np.arange(n)[:, None]) / n
    for j in range(d):
        u[:, j] = u[rng.permutation(n), j]
    return u


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing members.csv")
    args = ap.parse_args()
    cfg = load_config(args.config)

    params = cfg["params"]
    names = list(params)
    atm_in = (cfg["base_rundir"] / "atm_in").read_text()

    defaults = []
    for name in names:
        raw = read_nml_value(atm_in, name)
        if raw is None:
            sys.exit(f"'{name}' is not in the base atm_in; add it to the "
                     "base case user_nl_cam and rerun preview_namelists")
        defaults.append(fortran_to_float(raw))
        lo, hi = params[name]["range"]
        if not lo < hi:
            sys.exit(f"{name}: range must be [low, high] with low < high")
        if params[name].get("scale", "linear") == "log" and lo * hi <= 0:
            sys.exit(f"{name}: log scale needs a range of one sign")
        if not min(lo, hi) <= defaults[-1] <= max(lo, hi):
            print(f"warning: {name} default {defaults[-1]} is outside "
                  f"range {params[name]['range']}")

    n = int(cfg["sampling"]["n_members"])
    rng = np.random.default_rng(int(cfg["sampling"]["seed"]))
    u = lhs(n, len(names), rng)

    rows = [defaults]
    for k in range(n):
        row = []
        for j, name in enumerate(names):
            lo, hi = params[name]["range"]
            if params[name].get("scale", "linear") == "log":
                s = np.sign(lo)
                a, b = np.log(abs(lo)), np.log(abs(hi))
                row.append(float(s * np.exp(a + u[k, j] * (b - a))))
            else:
                row.append(float(lo + u[k, j] * (hi - lo)))
        rows.append(row)

    out = cfg["ppe_root"] / "members.csv"
    if out.exists() and not args.force:
        sys.exit(f"{out} exists; use --force to overwrite")
    cfg["ppe_root"].mkdir(parents=True, exist_ok=True)
    write_members(out, names, rows)
    print(f"wrote {out}: m0000 (defaults) + {n} LHS members, "
          f"{len(names)} parameters")


if __name__ == "__main__":
    main()
