#!/usr/bin/env python
"""Extract mountain cells of one run and variable to yearly files on the native grid.

Cells with mtn_frac > 0 are kept, so stricter thresholds can be applied later without re-extracting.
Values are copied untouched: no regridding, no unit change, no lossy packing.

Usage: 04_extract.py --domain EUR-12 --driver ERA5 --exp evaluation --var pr [--years 2000 2000] [--source v20250415]
"""
import argparse
import glob
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import cftime
import netCDF4 as nc
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from mtn_common import EXTRACTED as OUTROOT, MTN, run_dir

DATE = re.compile(r"_(\d{4})\d{4}-(\d{4})\d{4}\.nc$")
BOUNDS = {"pr": (0.0, 0.1), "tas": (170.0, 335.0)}  # kg m-2 s-1 (0.1 = 8640 mm/day), K


def sources(run, freq, var, choice):
    """Map source label -> files. Refuses to guess when several labels hold data."""
    d = f"{run}/{freq}/{var}"
    groups = {}
    flat = sorted(glob.glob(f"{d}/*.nc"))
    if flat:
        groups["flat"] = flat
    for v in sorted(glob.glob(f"{d}/*/")):
        fs = sorted(glob.glob(f"{v}*.nc"))
        if fs:
            groups[os.path.basename(v.rstrip("/"))] = fs
    if not groups:
        raise FileNotFoundError(d)
    if choice:
        return choice, groups[choice]
    if len(groups) > 1 and "flat" in groups:
        # AFR-12: the flat monthly files and v20260910 hold identical values (checked 2026-11-07),
        # and only the flat set is complete (EC-Earth3-Veg historical v20260910 lacks 1971-2010).
        print(f"several sources {sorted(groups)}; using flat monthly files", flush=True)
        return "flat", groups["flat"]
    if len(groups) > 1:
        raise SystemExit(f"Several sources in {d}: { {k: len(v) for k, v in groups.items()} }. Pick one with --source.")
    return next(iter(groups.items()))


def run_orog(run, fallback):
    fs = sorted(glob.glob(f"{run}/fx/orog/*.nc") + glob.glob(f"{run}/fx/orog/*/*.nc"))
    if fs:
        return np.ma.filled(nc.Dataset(fs[0])["orog"][:].astype("f4"), np.nan).squeeze(), "run"
    return fallback, "domain mask"


def extract_year(a):
    year, files, A = a
    out = A["outdir"] + f"/{A['var']}_{A['domain']}_{A['driver']}_{A['exp']}_mtn_{A['freq']}_{year}.nc"
    if os.path.exists(out) and not A["overwrite"]:
        return dict(year=year, status="exists")
    with nc.Dataset(A["mask"]) as m:
        keep = m["mtn_frac"][:] > 0
        mlat, mlon = m["lat"][:], m["lon"][:]
        frac, gid, prim = m["mtn_frac"][:], m["gmba_id"][:], m["is_primary"][:]
        mask_orog = m["orog"][:].filled(np.nan) if "orog" in m.variables else None
    orog, orog_src = run_orog(A["run"], mask_orog)
    ny, nx = keep.shape
    chunks, tvals, tbnds, used, attrs = [], [], [], [], None
    units = cal = None
    for f in files:
        y0, y1 = map(int, DATE.search(f).groups())
        if y1 < year or y0 > year:
            continue
        with nc.Dataset(f) as d:
            la, lo = d["lat"][:], d["lon"][:]
            if la.shape != (ny, nx) or np.abs(la - mlat).max() > 1e-3 or np.abs((lo - mlon + 180) % 360 - 180).max() > 1e-3:
                raise RuntimeError(f"grid mismatch with mask: {f}")
            t = d["time"]
            if units is None:
                units, cal = t.units, t.calendar
                attrs = {k: d[A["var"]].getncattr(k) for k in d[A["var"]].ncattrs() if k not in ("_FillValue", "missing_value", "grid_mapping")}
            elif (t.units, t.calendar) != (units, cal):
                raise RuntimeError(f"time encoding changed in {f}")
            tv = t[:].astype("f8")
            yrs = np.array([x.year for x in cftime.num2date(tv, units, cal)])
            sel = np.where(yrs == year)[0]
            if not sel.size:
                continue
            v = d[A["var"]]
            v.set_auto_mask(True)
            data = np.ma.filled(v[sel[0]:sel[-1] + 1].astype("f4"), np.nan)
            data[:, ~keep] = np.nan
            chunks.append(data)
            tvals.append(tv[sel])
            if "time_bnds" in d.variables:
                tbnds.append(d["time_bnds"][sel[0]:sel[-1] + 1].astype("f8"))
            used.append(os.path.basename(f))
    if not chunks:
        return dict(year=year, status="no data")
    data = np.concatenate(chunks)
    tv = np.concatenate(tvals)
    order = np.argsort(tv)
    data, tv = data[order], tv[order]
    n_dup = int((np.diff(tv) == 0).sum())
    if n_dup:
        raise RuntimeError(f"{year}: {n_dup} duplicate time steps")
    tb = np.concatenate(tbnds)[order] if tbnds else None
    # Daily expected counts per calendar; other frequencies are reported but not judged.
    exp_days = {"noleap": 365, "365_day": 365, "360_day": 360}.get(cal, 366 if (year % 4 == 0 and (year % 100 or year % 400 == 0)) else 365)
    gaps = int((np.diff(tv) > 1.5 * np.median(np.diff(tv))).sum()) if len(tv) > 1 else 0
    vals = data[:, keep]
    nmiss = int(np.isnan(vals).sum())
    lo, hi = BOUNDS.get(A["var"], (-np.inf, np.inf))
    nout = int(((vals < lo) | (vals > hi)).sum())
    qc = dict(year=year, status="ok", nt=len(tv), expected=exp_days if A["freq"] == "day" else None, gaps=gaps, ncells=int(keep.sum()),
              missing=nmiss, out_of_range=nout, vmin=float(np.nanmin(vals)), vmax=float(np.nanmax(vals)), vmean=float(np.nanmean(vals)), files=len(used))

    tmp = out + f".tmp{os.getpid()}"
    with nc.Dataset(tmp, "w", format="NETCDF4") as o:
        o.createDimension("time", None); o.createDimension("rlat", ny); o.createDimension("rlon", nx); o.createDimension("bnds", 2)
        o.setncatts(dict(title=f"{A['var']} over GMBA mountains (mtn_frac>0), native {A['domain']} grid", domain_id=A["domain"], driving_source_id=A["driver"],
                         driving_experiment_id=A["exp"], source_id="RegCM5-0", frequency=A["freq"], source_label=A["source"],
                         source_files=",".join(used), note="Use 2D lat/lon for positions. Values unchanged from source."))
        t = o.createVariable("time", "f8", ("time",)); t[:] = tv; t.setncatts(dict(units=units, calendar=cal, standard_name="time", axis="T", bounds="time_bnds"))
        if tb is not None:
            o.createVariable("time_bnds", "f8", ("time", "bnds"))[:] = tb
        def stat(name, val, typ, **att):
            x = o.createVariable(name, typ, ("rlat", "rlon"), zlib=True, complevel=4); x[:] = val; x.setncatts(att)
        stat("lat", mlat, "f8", units="degrees_north", standard_name="latitude")
        stat("lon", mlon, "f8", units="degrees_east", standard_name="longitude")
        stat("mtn_frac", frac, "f4", units="1", long_name="fraction of cell inside GMBA mountains")
        stat("gmba_id", gid, "i4", long_name="dominant GMBA Basic unit")
        stat("is_primary", prim, "i1", long_name="dominant unit assigned to this domain and cell is interior")
        if orog is not None:
            stat("orog", orog, "f4", units="m", long_name="model orography", source=orog_src)
        v = o.createVariable(A["var"], "f4", ("time", "rlat", "rlon"), zlib=True, complevel=1, shuffle=True,
                             chunksizes=(min(31, len(tv)), min(64, ny), min(64, nx)), fill_value=1e20)
        v.setncatts(attrs); v.coordinates = "lat lon"
        v[:] = np.ma.masked_invalid(data)
    os.replace(tmp, out)
    qc["MB"] = round(os.path.getsize(out) / 1e6, 1)
    return qc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--domain", required=True); p.add_argument("--driver", required=True); p.add_argument("--exp", required=True)
    p.add_argument("--var", required=True); p.add_argument("--freq", default="day"); p.add_argument("--years", nargs=2, type=int)
    p.add_argument("--source"); p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    if a.driver.startswith("wrong_"):
        raise SystemExit("wrong_* runs are withdrawn by the data owner")
    run = run_dir(a.domain, f"{a.driver}/{a.exp}")
    label, files = sources(run, a.freq, a.var, a.source)
    span = sorted({y for f in files for y in range(int(DATE.search(f)[1]), int(DATE.search(f)[2]) + 1)})
    years = [y for y in span if not a.years or a.years[0] <= y <= a.years[1]]
    outdir = f"{OUTROOT}/{a.domain}/{a.driver}/{a.exp}/{a.freq}/{a.var}"
    os.makedirs(outdir, exist_ok=True)
    A = dict(domain=a.domain, driver=a.driver, exp=a.exp, var=a.var, freq=a.freq, outdir=outdir, run=run, overwrite=a.overwrite,
             mask=f"{MTN}/masks/{a.domain}_mountain_mask.nc", source=label)
    print(f"{a.domain} {a.driver} {a.exp} {a.var} {a.freq}: source={label}, {len(files)} files, years {years[0]}-{years[-1]}", flush=True)
    bad = 0
    qclog = f"{outdir}/qc.jsonl"
    with ProcessPoolExecutor(a.workers) as ex, open(qclog, "a") as q:
        for r in ex.map(extract_year, [(y, files, A) for y in years]):
            print(json.dumps(r), flush=True)
            q.write(json.dumps(r) + "\n"); q.flush()
            if r.get("status") == "ok" and (r["out_of_range"] or r["missing"] or (r["expected"] and r["nt"] != r["expected"])):
                bad += 1
    print(f"finished; years flagged by QC: {bad}", flush=True)


if __name__ == "__main__":
    main()
