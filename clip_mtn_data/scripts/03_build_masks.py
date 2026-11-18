#!/usr/bin/env python
"""Build GMBA mountain masks on each CORDEX domain's native grid and resolve domain overlaps.

Cell positions come from the 2D lat/lon in the files. rlat/rlon are NOT used: they sit half a
cell away from lat/lon (the orography match confirms lat/lon), are corrupt in some SAM-12 ERA5
files, and in AUS/SAM do not agree with the crs pole.

Products:
  masks/<domain>_mountain_mask.nc   per cell: mtn_frac, gmba_id, edge_dist, is_primary, lat, lon, orog, sftlf
  masks/gmba_units_by_domain.csv    per GMBA Basic unit x domain: area inside domain / interior
  masks/gmba_unit_primary.csv       per GMBA Basic unit: primary domain and hierarchy names
  masks/overlap_area_by_domains.csv mountain area by the set of domains covering it
Also a 0.25 deg global overlap map in masks/, and rasters in $SCRATCH/mountains/rasters.
"""
import os
import sys

import geopandas as gpd
import netCDF4 as nc
import numpy as np
import pandas as pd
import rasterio.features
import rasterio.transform
import shapely

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from mtn_common import DEM, DOMAINS_12KM, GRID_RUN, MTN, OROG_RUN, SHP, first_file, geo2rot

BUFFER = 0  # no lateral-boundary buffer: agreed in the group meeting to use every cell of each domain
NY, NX, RES = 21600, 43200, 1 / 120
BLOCK = 720
A0 = (RES * 111.195) ** 2  # km2 of a 30s pixel at the equator
OUT = f"{MTN}/masks"
RAST = os.path.join(os.environ.get("SCRATCH", "/tmp"), "mountains", "rasters")
SHIFTS = (-0.5, 0.0, 0.5)


def wrap(a):
    return (a + 180) % 360 - 180


def load_grid(dom):
    """Fit a regular rotated grid to the file's 2D lat/lon, trying the crs pole and its antipode."""
    d = nc.Dataset(first_file(dom, GRID_RUN[dom], "day", "pr"))
    lat, lon = d["lat"][:].astype("f8").filled(np.nan), wrap(d["lon"][:].astype("f8").filled(np.nan))
    pa, po = float(d["crs"].grid_north_pole_latitude), float(d["crs"].grid_north_pole_longitude)
    ny, nx = lat.shape
    J, I = np.mgrid[0:ny, 0:nx]
    best = None
    for a, b in ((pa, po), (-pa, (po + 180) % 360)):
        xr, yr = geo2rot(lon, lat, a, b)
        rx, ry = wrap(xr - xr[0, 0]), yr - yr[0, 0]
        sx, sy = rx[0, -1] / (nx - 1), ry[-1, 0] / (ny - 1)
        res_x, res_y = rx - I * sx, ry - J * sy
        err = max(np.abs(res_x).max(), np.abs(res_y).max())
        if best is None or err < best["err"]:
            best = dict(err=err, plat=a, plon=b, x0=xr[0, 0] + res_x.mean(), y0=yr[0, 0] + res_y.mean(), sx=sx, sy=sy)
    print(f"{dom}: {ny}x{nx} pole=({best['plat']:.2f},{best['plon']:.2f}) step=({best['sx']:.5f},{best['sy']:.5f}) "
          f"fit residual max={best['err']:.5f} deg", flush=True)
    if best["err"] > 0.01:
        raise RuntimeError(f"{dom}: lat/lon are not a regular rotated grid (residual {best['err']:.4f} deg)")
    return dict(lat=lat, lon=lon, ny=ny, nx=nx, **best)


def load_fx(dom, var):
    try:
        v = nc.Dataset(first_file(dom, OROG_RUN[dom], "fx", var))[var]
        return np.ma.filled(v[:].astype("f8"), np.nan).squeeze()
    except (FileNotFoundError, IndexError):
        return None


def build_rasters():
    """Rasterize GMBA at the ETOPO grid, rows ordered south to north like the DEM file."""
    os.makedirs(RAST, exist_ok=True)
    fid, fany, flu = (f"{RAST}/gmba_basic_idx.npy", f"{RAST}/gmba_any.npy", f"{RAST}/gmba_basic_lookup.csv")
    g = gpd.read_file(SHP)
    bad = ~g.geometry.is_valid
    print(f"repairing {int(bad.sum())} invalid geometries", flush=True)
    g.loc[bad, "geometry"] = shapely.make_valid(g.geometry[bad].values)
    if not all(map(os.path.exists, (fid, fany, flu))):
        tr = rasterio.transform.from_origin(-180, 90, RES, RES)
        basic = g[g.MapUnit == "Basic"].reset_index(drop=True)
        # uint16 index: 0 = no unit, n = row n-1 of the lookup.
        print("rasterizing Basic units", flush=True)
        idx = rasterio.features.rasterize(
            ((geom, i + 1) for i, geom in enumerate(basic.geometry)), out_shape=(NY, NX), transform=tr, fill=0, dtype="uint16")
        print("rasterizing full union", flush=True)
        anym = rasterio.features.rasterize(
            ((geom, 1) for geom in g.geometry), out_shape=(NY, NX), transform=tr, fill=0, dtype="uint8")
        np.save(fid, idx[::-1])
        np.save(fany, anym[::-1])
        cols = ["GMBA_V2_ID", "MapName", "Hier_Lvl", "Area", "Elev_Low", "Elev_High", "Path"] + [f"Level_{i:02d}" for i in range(1, 6)]
        basic[cols].to_csv(flu, index_label="idx0")
    return np.load(fid), np.load(fany), pd.read_csv(flu)


def main():
    os.makedirs(OUT, exist_ok=True)
    idx, anym, lut = build_rasters()
    dem = nc.Dataset(DEM)
    dlat, dlon = np.asarray(dem["lat"][:]), np.asarray(dem["lon"][:])
    LON1 = dlon

    # Self-test of raster orientation and hierarchy: if these fail, nothing below is valid.
    for name, x, y in [("Everest", 86.925, 27.988), ("Kilimanjaro", 37.356, -3.067), ("Aconcagua", -70.011, -32.653), ("Paris", 2.35, 48.85), ("Amazon", -60, -3)]:
        i, j = np.abs(dlat - y).argmin(), np.abs(dlon - x).argmin()
        k = int(idx[i, j])
        print(f"raster test {name:12s} any={anym[i, j]} z={float(dem['z'][i, j]):7.0f} unit={lut.MapName[k - 1] if k else '-'}", flush=True)
    cos_all = np.cos(np.deg2rad(dlat))
    pix = A0 * cos_all[:, None]
    print(f"raster area: union {(anym * pix).sum()/1e6:.2f} M km2; pixels in union without Basic unit "
          f"{((anym > 0) & (idx == 0)).sum()/1e6:.2f} M px, area {(((anym > 0) & (idx == 0)) * pix).sum()/1e6:.3f} M km2", flush=True)
    del pix

    bits = np.zeros((NY, NX), dtype="uint8")  # bit k = pixel lies in the interior of domain k
    unit_rows, grids = [], {}
    for bit, dom in enumerate(DOMAINS_12KM):
        G = load_grid(dom)
        ny, nx = G["ny"], G["nx"]
        JJ, II = np.mgrid[0:ny, 0:nx]
        edge = np.minimum.reduce([JJ, ny - 1 - JJ, II, nx - 1 - II]).astype("i4")
        orog = load_fx(dom, "orog")
        sftlf = load_fx(dom, "sftlf")
        area = np.zeros(ny * nx)
        mtn = np.zeros(ny * nx)
        pairs = []  # (cell * 20000 + unit index, area) accumulated over blocks
        sums = {(a, b): (np.zeros(ny * nx), np.zeros(ny * nx)) for a in SHIFTS for b in SHIFTS}
        for i0 in range(0, NY, BLOCK):
            i1 = min(i0 + BLOCK, NY)
            LON, LAT = np.meshgrid(LON1, dlat[i0:i1])
            x, y = geo2rot(LON, LAT, G["plat"], G["plon"])
            del LON, LAT
            fi = wrap(x - G["x0"]) / G["sx"]
            fj = (y - G["y0"]) / G["sy"]
            del x, y
            ii, jj = np.rint(fi).astype("i4"), np.rint(fj).astype("i4")
            ok = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny)
            near = (fi >= -1) & (fi < nx) & (fj >= -1) & (fj < ny)
            if not near.any():
                continue
            rows = np.broadcast_to(cos_all[i0:i1, None], ok.shape)
            if ok.any():
                cell = jj[ok].astype("i8") * nx + ii[ok]
                w = A0 * rows[ok]
                area += np.bincount(cell, weights=w, minlength=ny * nx)
                am = anym[i0:i1][ok] > 0
                mtn += np.bincount(cell[am], weights=w[am], minlength=ny * nx)
                u = idx[i0:i1][ok]
                has = u > 0
                if has.any():
                    key, inv = np.unique(cell[has] * 20000 + u[has], return_inverse=True)
                    pairs.append((key, np.bincount(inv, weights=w[has])))
                interior = edge.ravel()[cell] >= BUFFER
                mask_px = np.zeros(ok.shape, bool)
                mask_px[ok] = interior & am
                bits[i0:i1][mask_px] |= 1 << bit
            # Registration test: DEM mean on the grid shifted by sub-cell offsets.
            z = np.maximum(np.ma.filled(dem["z"][i0:i1, :], 0.0), 0.0)
            for (sy, sx), (s, n) in sums.items():
                a, b = np.rint(fi[near] - sx).astype("i4"), np.rint(fj[near] - sy).astype("i4")
                k = (a >= 0) & (a < nx) & (b >= 0) & (b < ny)
                c = b[k].astype("i8") * nx + a[k]
                n += np.bincount(c, minlength=ny * nx)
                s += np.bincount(c, weights=z[near][k], minlength=ny * nx)
            print(f"  {dom} rows {i1}/{NY}", end="\r", flush=True)
        print(flush=True)

        area2d = area.reshape(ny, nx)
        frac = np.divide(mtn, area, out=np.zeros_like(mtn), where=area > 0).reshape(ny, nx)

        # Sum per (cell, unit) across blocks, then take each cell's dominant unit.
        K = np.concatenate([p[0] for p in pairs])
        W = np.concatenate([p[1] for p in pairs])
        key, inv = np.unique(K, return_inverse=True)
        wsum = np.bincount(inv, weights=W)
        cell_u, unit_u = key // 20000, (key % 20000).astype("i4")
        order = np.lexsort((wsum, cell_u))
        last = np.r_[cell_u[order][1:] != cell_u[order][:-1], True]
        dom_idx = np.zeros(ny * nx, dtype="i4")
        dom_idx[cell_u[order][last]] = unit_u[order][last]
        gid = np.full(ny * nx, 0, dtype="i4")
        gid[dom_idx > 0] = lut.GMBA_V2_ID.values[dom_idx[dom_idx > 0] - 1]

        df = pd.DataFrame(dict(unit=unit_u, a=wsum, e=edge.ravel()[cell_u]))
        df["inner"] = df.a.where(df.e >= BUFFER, 0.0)
        df["ew"] = df.a * df.e
        t = df.groupby("unit").agg(area_km2=("a", "sum"), interior_km2=("inner", "sum"), ew=("ew", "sum")).reset_index()
        t["mean_edge_dist"] = t.ew / t.area_km2
        t["domain"] = dom
        t["GMBA_V2_ID"] = lut.GMBA_V2_ID.values[t.unit.values - 1]
        unit_rows.append(t.drop(columns=["ew", "unit"]))

        # Registration score: expected best at shift (0,0) relative to 2D lat/lon.
        if orog is not None:
            gy, gx = np.gradient(orog)
            rough = (np.hypot(gx, gy) > 100) & (orog > 500) & (area2d > 0)
            line = []
            for k, (s, n) in sums.items():
                m = (s / np.maximum(n, 1)).reshape(ny, nx)
                line.append((np.sqrt(np.mean((m[rough] - orog[rough]) ** 2)), k))
            line.sort()
            zero = dict((k, v) for v, k in line)[(0.0, 0.0)]
            print(f"{dom} registration (n={int(rough.sum())}): best shift (rlat,rlon cells)={line[0][1]} RMSE={line[0][0]:.1f} m; "
                  f"RMSE at 2D lat/lon position={zero:.1f} m; worst={line[-1][0]:.1f} m", flush=True)
        print(f"{dom}: cells mtn_frac>0 {int((frac > 0).sum())}, >=0.5 {int((frac >= 0.5).sum())} of {ny * nx}; "
              f"mountain area {(frac * area2d).sum()/1e6:.3f} M km2, of it interior {(frac * area2d)[edge >= BUFFER].sum()/1e6:.3f}", flush=True)
        grids[dom] = dict(G=G, frac=frac, gid=gid.reshape(ny, nx), edge=edge, orog=orog, sftlf=sftlf, area=area2d)

    # Overlaps: mountain area by which domains' interiors cover it.
    pix = A0 * cos_all
    codes = {}
    for i0 in range(0, NY, BLOCK):
        b = bits[i0:i0 + BLOCK]; am = anym[i0:i0 + BLOCK] > 0
        w = np.broadcast_to(pix[i0:i0 + BLOCK, None], b.shape)[am]
        for c, s in zip(*(lambda bb: (np.unique(bb, return_inverse=True)[0], np.bincount(np.unique(bb, return_inverse=True)[1], weights=w)))(b[am])):
            codes[int(c)] = codes.get(int(c), 0.0) + s
    rows = [dict(domains="+".join(d for k, d in enumerate(DOMAINS_12KM) if c >> k & 1) or "NONE", area_Mkm2=s / 1e6) for c, s in codes.items()]
    ov = pd.DataFrame(rows).sort_values("area_Mkm2", ascending=False)
    ov["pct"] = 100 * ov.area_Mkm2 / ov.area_Mkm2.sum()
    ov.to_csv(f"{OUT}/overlap_area_by_domains.csv", index=False)
    print(ov.round(3).to_string(index=False), flush=True)

    # Primary domain per Basic unit: most interior area, near-ties go to the unit further from boundaries.
    u = pd.concat(unit_rows, ignore_index=True)
    u.to_csv(f"{OUT}/gmba_units_by_domain.csv", index=False)
    prim = []
    for gidv, h in u.groupby("GMBA_V2_ID"):
        h = h.sort_values("interior_km2", ascending=False)
        top = h[h.interior_km2 >= 0.99 * h.interior_km2.iloc[0]].sort_values("mean_edge_dist", ascending=False).iloc[0]
        prim.append(dict(GMBA_V2_ID=gidv, primary_domain=top.domain if top.interior_km2 > 0 else "NONE",
                         primary_interior_km2=top.interior_km2, n_domains=int((h.interior_km2 > 0).sum())))
    prim = lut.merge(pd.DataFrame(prim), on="GMBA_V2_ID", how="left")
    prim["primary_domain"] = prim.primary_domain.fillna("NONE")
    prim["primary_frac"] = prim.primary_interior_km2.fillna(0) / prim.Area
    prim.to_csv(f"{OUT}/gmba_unit_primary.csv", index=False)
    print("Basic units per primary domain:\n", prim.groupby("primary_domain").agg(n=("Area", "size"), Mkm2=("Area", lambda a: a.sum() / 1e6)).round(3).to_string(), flush=True)
    print("largest units with no interior coverage or <50% in their primary domain:", flush=True)
    print(prim[prim.primary_frac < 0.5].sort_values("Area", ascending=False)[["MapName", "Level_02", "Area", "primary_domain", "primary_frac"]].head(25).round(2).to_string(index=False), flush=True)

    pmap = dict(zip(prim.GMBA_V2_ID, prim.primary_domain))
    for dom, D in grids.items():
        G = D["G"]
        ny, nx = G["ny"], G["nx"]
        lookup = np.vectorize(lambda v: pmap.get(int(v), "NONE") == dom)
        primary = (D["gid"] > 0) & (D["edge"] >= BUFFER) & lookup(D["gid"])
        with nc.Dataset(f"{OUT}/{dom}_mountain_mask.nc", "w", format="NETCDF4") as o:
            o.createDimension("rlat", ny); o.createDimension("rlon", nx)
            o.setncatts(dict(title=f"GMBA v2.0 standard mountain mask on {dom} native grid", source_dem="ETOPO 2022 30s surface (raster reference only)",
                             interior_buffer_cells=BUFFER, fitted_pole=f"{G['plat']},{G['plon']}",
                             note="Cell positions are the 2D lat/lon. rlat/rlon in the source files are offset by half a cell and are not used."))
            def add(name, data, typ, **att):
                v = o.createVariable(name, typ, ("rlat", "rlon"), zlib=True, complevel=4)
                v[:] = data; v.setncatts(att)
            add("lat", G["lat"], "f8", units="degrees_north", standard_name="latitude")
            add("lon", G["lon"], "f8", units="degrees_east", standard_name="longitude")
            add("mtn_frac", D["frac"], "f4", long_name="fraction of cell area inside GMBA mountains", units="1")
            add("gmba_id", D["gid"], "i4", long_name="GMBA_V2_ID of the Basic unit with most area in the cell (0 = none)")
            add("edge_dist", D["edge"], "i4", long_name="cells to the nearest lateral boundary")
            add("is_primary", primary.astype("i1"), "i1", long_name="dominant unit assigned to this domain and cell is interior")
            add("cell_area", D["area"], "f4", units="km2", long_name="cell area summed from 30s pixels")
            if D["orog"] is not None:
                add("orog", D["orog"], "f4", units="m", long_name="model orography", source_run=OROG_RUN[dom])
            if D["sftlf"] is not None:
                add("sftlf", D["sftlf"], "f4", units="%", long_name="model land fraction", source_run=OROG_RUN[dom])
        print(f"wrote {dom}: primary cells {int(primary.sum())}", flush=True)

    # Coarse map for the notebook: 0.25 deg, mountain pixels only, per-domain interior fraction.
    f = 30
    with nc.Dataset(f"{OUT}/overlap_map_025deg.nc", "w") as o:
        o.createDimension("lat", NY // f); o.createDimension("lon", NX // f)
        o.createVariable("lat", "f4", ("lat"))[:] = dlat[f // 2::f][:NY // f]
        o.createVariable("lon", "f4", ("lon"))[:] = dlon[f // 2::f][:NX // f]
        am = (anym > 0).astype("f4")
        o.createVariable("mtn_frac", "f4", ("lat", "lon"), zlib=True)[:] = am.reshape(NY // f, f, NX // f, f).mean((1, 3))
        for k, dom in enumerate(DOMAINS_12KM):
            cov = ((bits >> k) & 1).astype("f4")
            o.createVariable(f"cov_{dom.replace('-', '_')}", "f4", ("lat", "lon"), zlib=True)[:] = cov.reshape(NY // f, f, NX // f, f).mean((1, 3))
    print("done", flush=True)


if __name__ == "__main__":
    main()
