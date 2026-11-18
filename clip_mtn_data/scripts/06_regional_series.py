#!/usr/bin/env python
"""Regional mountain-mean pr and tas per run: all mountain cells and elevation bands, monthly.

Weights are mtn_frac * cell_area, so each cell counts by its mountain area; cells with
mtn_frac > 0 are all used (foothills included). Bands use the model orography of the run.
Output: results/regional_series/<domain>_monthly.csv, one row per run/band/month:
  pr_mm  = monthly total (mm), tas_C = monthly mean (deg C), ndays = days present.
"""
import argparse
import glob
import os
import sys

import cftime
import netCDF4 as nc
import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from mtn_common import EXTRACTED, MTN

EDGES = [0, 1000, 2000, 3000, 4000, 5000, np.inf]
BANDS = ["0-1000", "1000-2000", "2000-3000", "3000-4000", "4000-5000", ">5000"]


def weights(f, mask):
    """(ncell, 1 + nband) weight matrix over kept cells, from the run's own orography."""
    with nc.Dataset(f) as d:
        orog = d["orog"][:].filled(np.nan) if "orog" in d.variables else mask["orog"]
    keep = mask["frac"] > 0
    w = (mask["frac"] * mask["area"])[keep]
    z = np.maximum(orog[keep], 0)  # sub-sea-level cells go to the lowest band
    band = np.digitize(z, EDGES[1:-1])
    W = np.zeros((w.size, 1 + len(BANDS)))
    W[:, 0] = w
    W[np.arange(w.size), 1 + band] = w
    return keep, W


def series(files, var, mask):
    keep = W = None
    out = []
    for f in files:
        with nc.Dataset(f) as d:
            if W is None:
                keep, W = weights(f, mask)
            x = d[var][:].filled(np.nan)[:, keep].astype("f8")
            t = d["time"]
            dates = cftime.num2date(t[:], t.units, t.calendar)
        valid = ~np.isnan(x)
        num = np.nan_to_num(x) @ W
        den = valid @ W
        m = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        out.append(pd.DataFrame(m, columns=["all"] + BANDS).assign(year=[x.year for x in dates], month=[x.month for x in dates]))
    df = pd.concat(out, ignore_index=True)
    return df, (W[:, 1:] > 0).sum(0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--domain", required=True)
    a = p.parse_args()
    with nc.Dataset(f"{MTN}/masks/{a.domain}_mountain_mask.nc") as m:
        mask = dict(frac=m["mtn_frac"][:].filled(0), area=m["cell_area"][:].filled(0), orog=m["orog"][:].filled(np.nan))
    rows, cells = [], []
    for rdir in sorted(glob.glob(f"{EXTRACTED}/{a.domain}/*/*/day")):
        driver, exp = rdir.split("/")[-3:-1]
        per_var = {}
        for var in ("pr", "tas"):
            files = sorted(glob.glob(f"{rdir}/{var}/{var}_*_mtn_day_*.nc"))
            if not files:
                continue
            df, nband = series(files, var, mask)
            long = df.melt(id_vars=["year", "month"], var_name="band", value_name="v")
            g = long.groupby(["year", "month", "band"], sort=False)
            agg = g.v.agg(["sum", "mean", "count"]).reset_index()
            agg.loc[agg["count"] == 0, ["sum", "mean"]] = np.nan  # band without cells: NaN, not a zero total
            if var == "pr":
                agg[var] = agg["sum"] * 86400.0  # kg m-2 s-1 daily means -> mm per month
            else:
                agg[var] = agg["mean"] - 273.15
            agg = agg.rename(columns={"count": f"ndays_{var}"})[["year", "month", "band", var, f"ndays_{var}"]]
            per_var[var] = agg
            cells.append(dict(domain=a.domain, driver=driver, experiment=exp, var=var, **dict(zip(BANDS, map(int, nband)))))
        if not per_var:
            continue
        merged = None
        for v in per_var.values():
            merged = v if merged is None else merged.merge(v, on=["year", "month", "band"], how="outer")
        merged.insert(0, "experiment", exp); merged.insert(0, "driver", driver); merged.insert(0, "domain", a.domain)
        rows.append(merged)
        print(f"{a.domain} {driver} {exp}: {merged.year.min()}-{merged.year.max()} {len(merged)} rows", flush=True)
    if not rows:
        raise SystemExit(f"no extracted data for {a.domain}")
    out = pd.concat(rows, ignore_index=True).rename(columns={"pr": "pr_mm", "tas": "tas_C"})
    out.to_csv(f"{MTN}/results/regional_series/{a.domain}_monthly.csv", index=False, float_format="%.5g")
    pd.DataFrame(cells).to_csv(f"{MTN}/results/regional_series/{a.domain}_band_cells.csv", index=False)
    print("wrote", f"{MTN}/results/regional_series/{a.domain}_monthly.csv", flush=True)


if __name__ == "__main__":
    main()
