# Clipping CORDEX-CMIP6 12 km output to the world's mountains

Code and masks to extract the mountain part of the CORDEX-CMIP6 dynamical downscaling runs made with
**RegCM5-0 (ICTP)** on the six 12 km domains **AFR-12, AUS-12, EAS-12, EUR-12, NAM-12, SAM-12**.
Mountains are defined by the GMBA Mountain Inventory v2.0 (standard) as in: https://www.earthenv.org/mountains .

The masks are already built and included here (`masks/`). If you have the same grids with RegCM5.0 output on
your own system, you do not need the DEM or the GMBA shapefile: start directly from step 2 below.


## Contents here

Path: What it is 
------
`masks/<domain>_mountain_mask.nc`: Per-cell mountain fraction (`mtn_frac`), GMBA unit, `is_primary`, cell area, model orography, on the native grid of each domain
`masks/*.csv`, `masks/overlap_map_025deg.nc`: GMBA units per domain, primary-domain assignment, overlap between domains
`masks/README.md`: **Full description of the mask variables and of the method.** Read this first. 
`scripts/`: Python and Slurm scripts (see the workflow below) 
`notebooks/01_mask_check.ipynb`: Maps and tables of the masks 
`notebooks/02_regional_timeseries.ipynb`: Annual and seasonal mountain-mean `pr` and `tas`, by domain and by per 1000 m elevation band 


## Requirements

Python 3.10 or later with the packages in `requirements.txt`:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Configuration

The defaults are the paths on Leonardo (CINECA). Elsewhere, set these environment variables:

Variable: Meaning | Default
|---|---|---|
`CORDEX_ROOT`: Root of the your output path | in Leonardo is: `/leonardo_work/ICT26_ESP/CORDEX-CMIP6/DD`
`MTN_OUT`: Where extracted files are written | in Leonardo is: `$SCRATCH/mountains/extracted`
`MTN_DIR`: This folder | the folder that contains `scripts/` 
`PY`: Python used by the Slurm scripts | `python3` 

An example of the file: `EUR-12/ICTP/MPI-ESM1-2-HR/ssp370/r1i1p1f1/RegCM5-0/v1-r1/day/pr/v20250415/*.nc`.
Daily file names must end in `_YYYYMMDD-YYYYMMDD.nc` (CORDEX convention), which is how the years are found. Revise these code If name is different.

## Workflow (If you have the same grids with RegCM5.0 output, you can skip step 1 and start directly from step 2.)

**1. Build the masks** (optional: they are already in `masks/`)

```bash
sbatch scripts/03_build_masks.sh
```

This needs the files in `DEM_and_shp/` and one `pr` file plus `fx/orog` and `fx/sftlf` per domain.

**2. Extract daily `pr` and `tas` over the mountains**

(a) All six domains: one job per domain, each followed by the regional series (step 4). 

```bash
bash scripts/submit_all.sh --no-mask
```

Keep `--no-mask`. Without it, the script first rebuilds the masks (step 1), which needs the files in
`DEM_and_shp/` and overwrites the masks shipped in `masks/`.

(b) One domain only (all its runs, then its regional series), for example to test the setup:

```bash
sbatch scripts/04_extract_domain.sh EUR-12
```

(c) One run and one variable, interactively or inside your own job:

```bash
python scripts/04_extract.py --domain EUR-12 --driver ERA5 --exp evaluation --var pr [--years 1990 2000]
```

What the extraction does:
- It keeps every cell with `mtn_frac > 0` and sets all other cells to missing. 
- Output is one file per year on the native grid:
  `$MTN_OUT/<domain>/<driver>/<exp>/day/<var>/<var>_<domain>_<driver>_<exp>_mtn_day_<year>.nc`.
  Each file also contains `lat`, `lon`, `mtn_frac`, `gmba_id`, `is_primary` and `orog`.
- Quality checks for each year (days against the calendar, gaps, missing values, physical range) are
  appended to `qc.jsonl` in the same folder.
- It can be restarted: years already written are skipped, so re-running after new data arrive adds only the new years.
- The script stops if the grid of a source file does not match the mask (tolerance 0.001°).

Runs at 12 km with all of 1950–2100 take tens of GB per domain for `pr` and `tas` together
(AFR-12 ≈ 40 GB, EAS-12 ≈ 44 GB, EUR-12 ≈ 28 GB).

**3. Check an extracted file against the raw data**

```bash
python scripts/05_verify_extract.py <extracted_file.nc> <raw_run_folder>
```

It confirms that cells outside the mask are missing and that cells inside are identical to the source.

**4. Regional series**

```bash
python scripts/06_regional_series.py --domain EUR-12
```

Writes `results/regional_series/<domain>_monthly.csv`: monthly totals of `pr` (mm) and monthly means of `tas` (°C),
averaged over the mountain cells with weights `mtn_frac × cell_area`. The series are given for all cells and for every 1000 m
bands of model orography. The notebook `02_regional_timeseries.ipynb` plots them.

## Using the masks with other data

The masks match the RegCM5-0 grids of these six domains only. For those may have different grids, rebuild the masks with `scripts/03_build_masks.py`: point
`GRID_RUN`/`OROG_RUN` in `scripts/mtn_common.py` to a run of that model.

To apply a mask yourself:

```python
import xarray as xr
m = xr.open_dataset("masks/EUR-12_mountain_mask.nc")
pr = xr.open_dataset("pr_EUR-12_..._day_20000101-20001231.nc")["pr"]
pr_mtn = pr.where(m["mtn_frac"].values > 0)            # same (rlat, rlon) grid
```

Use the 2D `lat`/`lon` for positions, not `rlat`/`rlon`, which are offset by half a cell in the source
files (see `masks/README.md`).

