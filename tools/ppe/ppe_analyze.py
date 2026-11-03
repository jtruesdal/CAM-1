#!/usr/bin/env python3
"""Compute nudging-tendency metrics for every finished PPE member.

For each member with a DONE marker, read the per-step history fields
Nudge_T and Nudge_Q, drop a spin-up window, and compute:
  mean_<F>(lev)  time-mean profile (systematic error the nudging corrects)
  rms_<F>(lev)   time RMS profile
  bias_<F>       mass-weighted RMS of mean_<F> over the [p_top, p_bot] band
  tot_<F>        mass-weighted mean of rms_<F> over the same band
and the ratio of each scalar to m0000's. Units: K/day and g/kg/day.
Writes <ppe_root>/ppe_metrics.nc (dims member, lev) with the parameter
values from members.csv, and prints a summary table.

usage: ppe_analyze.py ppe_config.toml
"""
import argparse
import glob
import sys

import numpy as np
import xarray as xr

from ppe_common import load_config, read_members

FIELDS = {"Nudge_T": 86400.0, "Nudge_Q": 86400.0e3}  # -> K/day, g/kg/day
UNITS = {"Nudge_T": "K/day", "Nudge_Q": "g/kg/day"}


def load_member(d, spinup_h):
    """Concatenate all history files holding Nudge_T; drop spin-up."""
    parts = []
    for f in sorted(glob.glob(f"{d}/*.cam.h*.nc")):
        ds = xr.open_dataset(f, decode_times=True)
        if "Nudge_T" in ds:
            keep = list(FIELDS) + ["PS", "hyai", "hybi", "P0"]
            parts.append(ds[[v for v in keep if v in ds]])
    if not parts:
        return None
    ds = xr.concat(parts, dim="time", data_vars="minimal")
    ds = ds.sortby("time")
    t = ds["time"].values
    hours = np.array([(x - t[0]).total_seconds() / 3600.0 for x in t])
    ds = ds.isel(time=np.nonzero(hours >= spinup_h)[0])
    return ds.squeeze("ncol", drop=True) if "ncol" in ds.dims else ds


def layer_weights(ds, p_top, p_bot):
    """Time-mean layer pressure thickness, zero outside [p_top, p_bot] hPa."""
    p0 = float(ds["P0"]) if "P0" in ds else 1.0e5
    ps = float(ds["PS"].mean()) if "PS" in ds else 1.0e5
    hyai, hybi = ds["hyai"].values, ds["hybi"].values
    if hyai.ndim > 1:  # concatenated along time
        hyai, hybi = hyai[0], hybi[0]
    pint = hyai * p0 + hybi * ps
    dp = np.diff(pint)
    pmid = 0.5 * (pint[1:] + pint[:-1]) / 100.0
    return np.where((pmid >= p_top) & (pmid <= p_bot), dp, 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config")
    args = ap.parse_args()
    cfg = load_config(args.config)
    root = cfg["ppe_root"]
    an = cfg.get("analysis", {})
    spinup_h = float(an.get("spinup_hours", 24.0))
    p_top, p_bot = an.get("p_range_hpa", [100.0, 1000.0])

    names, members = read_members(root / "members.csv")
    results, done = {}, []
    for m, vals in members:
        d = root / m
        if not (d / "DONE").exists():
            continue
        ds = load_member(d, spinup_h)
        if ds is None or ds.sizes.get("time", 0) == 0:
            print(f"{m}: no Nudge_T samples after {spinup_h} h spin-up")
            continue
        w = layer_weights(ds, p_top, p_bot)
        r = {"lev": ds["lev"].values, "nt": ds.sizes["time"]}
        for f, scale in FIELDS.items():
            x = ds[f].transpose("time", "lev").values * scale
            mean, rms = x.mean(axis=0), np.sqrt((x ** 2).mean(axis=0))
            r[f"mean_{f}"], r[f"rms_{f}"] = mean, rms
            r[f"bias_{f}"] = float(np.sqrt((w * mean ** 2).sum() / w.sum()))
            r[f"tot_{f}"] = float((w * rms).sum() / w.sum())
        results[m] = r
        done.append((m, vals))

    if not done:
        sys.exit("no finished members with usable output")
    if "m0000" not in results:
        sys.exit("m0000 (defaults) is missing; it is the reference")
    nts = {results[m]["nt"] for m, _ in done}
    if len(nts) > 1:
        print(f"warning: members have different sample counts {sorted(nts)}")

    mem = [m for m, _ in done]
    lev = results["m0000"]["lev"]
    data = {}
    for f in FIELDS:
        for kind in ("mean", "rms"):
            data[f"{kind}_{f}"] = (
                ("member", "lev"),
                np.stack([results[m][f"{kind}_{f}"] for m in mem]),
                {"units": UNITS[f]})
        for kind in ("bias", "tot"):
            v = np.array([results[m][f"{kind}_{f}"] for m in mem])
            data[f"{kind}_{f}"] = ("member", v, {"units": UNITS[f]})
            data[f"{kind}_{f}_rel"] = ("member",
                                       v / results["m0000"][f"{kind}_{f}"])
    for n in names:
        data[n] = ("member", np.array([v[n] for _, v in done]))
    out = xr.Dataset(data, coords={"member": mem, "lev": lev},
                     attrs={"spinup_hours": spinup_h,
                            "p_range_hpa": [p_top, p_bot],
                            "samples_per_member": max(nts)})
    path = root / "ppe_metrics.nc"
    out.to_netcdf(path)

    print(f"{len(mem)} members, {max(nts)} samples each after "
          f"{spinup_h:g} h spin-up, {p_top:g}-{p_bot:g} hPa")
    print(f"{'member':8s} {'biasT':>8s} {'totT':>8s} {'biasQ':>8s} "
          f"{'totQ':>8s}   (relative to m0000)")
    for i, m in enumerate(mem):
        print(f"{m:8s} " + " ".join(
            f"{float(out[k + '_rel'][i]):8.3f}"
            for k in ("bias_Nudge_T", "tot_Nudge_T",
                      "bias_Nudge_Q", "tot_Nudge_Q")))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
