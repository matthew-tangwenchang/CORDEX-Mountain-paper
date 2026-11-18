#!/usr/bin/env python
"""Compare an extracted yearly file with the raw source, independently of the extraction code path."""
import glob
import sys

import netCDF4 as nc
import numpy as np


def main(path, srcdir):
    o = nc.Dataset(path)
    var = o.getncattr("title").split()[0]
    files = o.source_files.split(",")
    raw, rt = [], []
    for f in files:
        d = nc.Dataset(glob.glob(f"{srcdir}/**/{f}", recursive=True)[0])
        raw.append(d[var][:].filled(np.nan)); rt.append(d["time"][:])
    raw, rt = np.concatenate(raw).astype("f4"), np.concatenate(rt)
    ot = o["time"][:]
    i = np.searchsorted(rt, ot)
    assert (rt[i] == ot).all(), "time mismatch"
    raw = raw[i]
    x = o[var][:].filled(np.nan)
    keep = o["mtn_frac"][:] > 0
    print(f"{var}: nt={len(ot)} kept cells {int(keep.sum())} of {keep.size}")
    print("outside mask all NaN:", bool(np.isnan(x[:, ~keep]).all()))
    print("inside mask identical to raw (NaN-aware):", bool(np.array_equal(x[:, keep], raw[:, keep], equal_nan=True)))
    print(f"source sample means kept/all: {np.nanmean(raw[:, keep]):.6g} / {np.nanmean(raw):.6g}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
