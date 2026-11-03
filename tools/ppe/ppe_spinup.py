#!/usr/bin/env python3
"""Plot how a member's nudging tendencies and aerosols evolve after startup.

Reads every history file in a member directory that holds Nudge_T and
writes <member_dir>/spinup.png with, against hours since the first sample:
  - mass-weighted mean and RMS of Nudge_T over [p_top, p_bot] (K/day)
  - the same for Nudge_Q (g/kg/day)
  - mass-weighted mean of num_a1..num_a4 relative to their first value
Also prints 6-hourly averages of those series so the spin-up can be read
without opening the plot.

usage: ppe_spinup.py MEMBER_DIR [--p-range 100 1000]
"""
import argparse
import glob
import sys

import numpy as np
import xarray as xr

AEROSOLS = ["num_a1", "num_a2", "num_a3", "num_a4"]


def load(d):
    parts = []
    for f in sorted(glob.glob(f"{d}/*.cam.h*.nc")):
        ds = xr.open_dataset(f)
        if "Nudge_T" in ds:
            keep = ["Nudge_T", "Nudge_Q", "PS", "hyai", "hybi", "P0"] + \
                   [a for a in AEROSOLS if a in ds]
            parts.append(ds[[v for v in keep if v in ds]])
    if not parts:
        sys.exit(f"no history files with Nudge_T in {d}")
    ds = xr.concat(parts, dim="time", data_vars="minimal").sortby("time")
    return ds.squeeze("ncol", drop=True) if "ncol" in ds.dims else ds


def weights(ds, p_top, p_bot):
    p0 = float(ds["P0"]) if "P0" in ds else 1.0e5
    hyai, hybi = ds["hyai"].values, ds["hybi"].values
    if hyai.ndim > 1:
        hyai, hybi = hyai[0], hybi[0]
    ps = ds["PS"].values.reshape(-1, 1) if "PS" in ds else np.full((1, 1),
                                                                   1.0e5)
    pint = hyai[None, :] * p0 + hybi[None, :] * ps          # (time, ilev)
    dp = np.diff(pint, axis=1)
    pmid = 0.5 * (pint[:, 1:] + pint[:, :-1]) / 100.0
    w = np.where((pmid >= p_top) & (pmid <= p_bot), dp, 0.0)
    return w / w.sum(axis=1, keepdims=True)                 # (time, lev)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("member_dir")
    ap.add_argument("--p-range", nargs=2, type=float, default=[100., 1000.])
    args = ap.parse_args()

    ds = load(args.member_dir)
    t = ds["time"].values
    hours = np.array([(x - t[0]).total_seconds() / 3600.0 for x in t])
    w = weights(ds, *args.p_range)

    series = {}
    for f, scale, unit in (("Nudge_T", 86400.0, "K/day"),
                           ("Nudge_Q", 86400.0e3, "g/kg/day")):
        x = ds[f].transpose("time", "lev").values * scale
        series[f"{f} mean ({unit})"] = (w * x).sum(axis=1)
        series[f"{f} rms ({unit})"] = np.sqrt((w * x ** 2).sum(axis=1))
    for a in AEROSOLS:
        if a in ds:
            x = (w * ds[a].transpose("time", "lev").values).sum(axis=1)
            series[f"{a} / initial"] = x / x[0]

    # 6-hourly text summary
    names = list(series)
    print(f"{len(hours)} samples, {hours[-1]:.1f} h, "
          f"{args.p_range[0]:g}-{args.p_range[1]:g} hPa mass-weighted")
    print(f"{'hours':>9s} " + " ".join(f"{n.split(' (')[0]:>16s}"
                                        for n in names))
    for h0 in np.arange(0.0, hours[-1] + 1e-6, 6.0):
        sel = (hours >= h0) & (hours < h0 + 6.0)
        if sel.any():
            print(f"{h0:4.0f}-{h0 + 6:<4.0f} " + " ".join(
                f"{series[n][sel].mean():16.4g}" for n in names))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipped plot")
        return
    groups = [[n for n in names if n.startswith("Nudge_T")],
              [n for n in names if n.startswith("Nudge_Q")],
              [n for n in names if n.startswith("num_a")]]
    groups = [g for g in groups if g]
    fig, axes = plt.subplots(len(groups), 1, figsize=(9, 3 * len(groups)),
                             sharex=True)
    for ax, g in zip(np.atleast_1d(axes), groups):
        for n in g:
            ax.plot(hours, series[n], label=n)
        ax.axhline(0.0 if "num_a" not in g[0] else 1.0, color="0.7", lw=0.8)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    np.atleast_1d(axes)[-1].set_xlabel("hours since first sample")
    fig.suptitle(args.member_dir)
    fig.tight_layout()
    out = f"{args.member_dir.rstrip('/')}/spinup.png"
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
