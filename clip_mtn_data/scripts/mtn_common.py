"""Shared helpers for the global-mountain CORDEX-CMIP6 workflow."""
import glob
import os

import numpy as np

# Paths. Defaults are the ones used on Leonardo (CINECA); override them with environment variables elsewhere.
ROOT = os.environ.get("CORDEX_ROOT", "/leonardo_work/ICT26_ESP/CORDEX-CMIP6/DD")  # CORDEX-CMIP6 DD archive
MTN = os.environ.get("MTN_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # this folder
DEM = os.environ.get("MTN_DEM", f"{MTN}/DEM_and_shp/ETOPO_2022_v1_30s_N90W180_surface.nc")
SHP = os.environ.get("MTN_SHP", f"{MTN}/DEM_and_shp/GMBA_Inventory_v2.0_standard/GMBA_Inventory_v2.0_standard.shp")
# Extracted mountain data (large): $MTN_OUT, else $SCRATCH/mountains/extracted, else <MTN>/extracted.
EXTRACTED = os.environ.get("MTN_OUT") or (os.path.join(os.environ["SCRATCH"], "mountains", "extracted")
                                          if os.environ.get("SCRATCH") else f"{MTN}/extracted")
DOMAINS_12KM = ["AFR-12", "AUS-12", "EAS-12", "EUR-12", "NAM-12", "SAM-12"]
# A run whose rlat/rlon are intact and one holding fx/orog, per domain.
GRID_RUN = {
    "AFR-12": "ERA5/evaluation", "AUS-12": "MPI-ESM1-2-HR/historical", "EAS-12": "ERA5/evaluation",
    "EUR-12": "ERA5/evaluation", "NAM-12": "ERA5/evaluation", "SAM-12": "MPI-ESM1-2-HR/historical",
}
OROG_RUN = {
    "AFR-12": "ERA5/evaluation", "AUS-12": "MPI-ESM1-2-HR/historical", "EAS-12": "ERA5/evaluation",
    "EUR-12": "ERA5/evaluation", "NAM-12": "ERA5/evaluation", "SAM-12": "ERA5/evaluation",
}


def run_dir(domain, run):
    return f"{ROOT}/{domain}/ICTP/{run}/r1i1p1f1/RegCM5-0/v1-r1"


def first_file(domain, run, freq, var):
    d = f"{run_dir(domain, run)}/{freq}/{var}"
    fs = sorted(glob.glob(f"{d}/*.nc") + glob.glob(f"{d}/*/*.nc"))
    if not fs:
        raise FileNotFoundError(d)
    return fs[0]


def _basis(plat, plon):
    return np.deg2rad(plat), np.deg2rad(plon) - np.pi


def rot2geo(rlon, rlat, plat, plon):
    """Rotated (deg) -> geographic (deg, lon in [-180, 180)). Pole is the grid north pole."""
    pa, l0 = _basis(plat, plon)
    rlo, rla = np.deg2rad(rlon), np.deg2rad(rlat)
    xp, yp, zp = np.cos(rla) * np.cos(rlo), np.cos(rla) * np.sin(rlo), np.sin(rla)
    x = xp * np.sin(pa) - zp * np.cos(pa)
    z = zp * np.sin(pa) + xp * np.cos(pa)
    lon = np.arctan2(yp, x) + l0
    return (np.rad2deg(lon) + 180) % 360 - 180, np.rad2deg(np.arcsin(np.clip(z, -1, 1)))


def geo2rot(lon, lat, plat, plon):
    """Geographic (deg) -> rotated (deg, rlon in [-180, 180))."""
    pa, l0 = _basis(plat, plon)
    lo, la = np.deg2rad(lon) - l0, np.deg2rad(lat)
    x, y, z = np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)
    xp = x * np.sin(pa) + z * np.cos(pa)
    zp = z * np.sin(pa) - x * np.cos(pa)
    return np.rad2deg(np.arctan2(y, xp)), np.rad2deg(np.arcsin(np.clip(zp, -1, 1)))
