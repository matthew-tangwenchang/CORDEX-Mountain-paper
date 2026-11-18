# Mountain masks for CORDEX-CMIP6 12 km domains (RegCM5-0)

The masks mark which grid cells of each 12 km CORDEX-CMIP6 domain lie in mountains, as defined by the
GMBA Mountain Inventory v2.0 (standard). They are built by `scripts/03_build_masks.py`. The scripts
in `scripts/04_extract*.py/.sh` use them to extract daily `pr` and `tas`.

Domains: AFR-12, AUS-12, EAS-12, EUR-12, NAM-12, SAM-12.

## Inputs

| Input | File |
|---|---|
| Mountain definition | `DEM_and_shp/GMBA_Inventory_v2.0_standard/GMBA_Inventory_v2.0_standard.shp` (8327 polygons, WGS84) |
| Reference DEM | `DEM_and_shp/ETOPO_2022_v1_30s_N90W180_surface.nc` (30 arc-second) |
| Model grids | 2D `lat`/`lon` of the daily `pr` files of each domain |
| Model orography and land fraction | `fx/orog`, `fx/sftlf` of the domain |

The ETOPO DEM is used only to check grid registration. It is not used to define mountains or elevation
bands. Elevation analyses use the model orography in RegCM `orog`, because `tas` and `pr` were simulated on it.

## Files

| File | Content |
|---|---|
| `<domain>_mountain_mask.nc` | Per-cell mask variables on the native grid of the domain (see below) |
| `gmba_unit_primary.csv` | One row per GMBA Basic unit: name, hierarchy (`Level_01`–`Level_05`), area, assigned `primary_domain` (`NONE` = not covered) |
| `gmba_units_by_domain.csv` | Area of each GMBA Basic unit inside each domain |
| `overlap_area_by_domains.csv` | Mountain area (M km²) by the set of domains that cover it |
| `overlap_map_025deg.nc` | 0.25° global map: GMBA mountain fraction and coverage fraction per domain (for plotting) |

## Variables in `<domain>_mountain_mask.nc`

All variables are 2D on `(rlat, rlon)`, the same grid as the model output of the domain.

| Variable | Units | Meaning |
|---|---|---|
| `lat`, `lon` | degrees | Cell centre coordinates. **Use these for location and plotting.** `lon` is in [-180, 180). |
| `mtn_frac` | 1 | Fraction of the cell area inside GMBA mountains (0–1). The denominator is the whole cell area, including any sea. |
| `gmba_id` | – | `GMBA_V2_ID` of the GMBA Basic unit with the largest area in the cell; 0 = no mountain. Names and hierarchy are in `gmba_unit_primary.csv`. |
| `edge_dist` | cells | Distance to the nearest lateral boundary of the domain (0 = outermost row/column). Kept for reference only (see "No boundary buffer"). |
| `is_primary` | 0/1 | 1 if the cell's `gmba_id` unit is assigned to this domain. Used to avoid double counting where domains overlap (see below). |
| `cell_area` | km² | Cell area, summed from the 30 arc-second pixels that fall into the cell. |
| `orog` | m | Model surface height (`fx/orog`). |
| `sftlf` | % | Model land fraction (`fx/sftlf`). |

## How the mask is computed (helped by claude code)

1. **Rasterise GMBA at 30 arc-seconds (~1 km).** The 8327 polygons are rasterised on the ETOPO grid, with
   a pixel counted as inside when its centre is inside. Two layers are made:
   - the union of all polygons, which says whether a pixel is mountain, and
   - the index of the GMBA **Basic** unit containing the pixel.

   GMBA polygons are nested (for example Eurasia > Central Asia > Himalaya > ...), so the polygon areas cannot be summed: their sum is 142 M km², the union only 24.8 M km². The Basic units cover the union with almost no overlap (24.85 M km²), so they are used for unit identity. Before rasterising,22 polygons with ring self-intersections are repaired with `shapely.make_valid` (in memory only).
2. **Locate each 30 arc-second pixel in the model grid.** The rotated-pole grid is fitted to the 2D
   `lat`/`lon` of the model files (fit residual < 0.0001°). Each pixel is converted to rotated coordinates and assigned to the nearest cell centre.
3. **Aggregate per model cell** with pixel areas weighted by cos(latitude):
   `mtn_frac` = mountain area / cell area; `gmba_id` = Basic unit with the most area in the cell.

### Why `rlat`/`rlon` are not used

In the source files, `rlat`/`rlon` are offset by **half a grid cell** from the 2D `lat`/`lon` in every domain. Some SAM-12 ERA5 files contain corrupt `rlat`/`rlon` values. In AUS-12 and SAM-12, the `crs` pole attributes are the antipode of the pole that matches `lat`/`lon`. The 2D `lat`/`lon` are the true cell positions: ETOPO aggregated to the `lat`/`lon` positions matches the model `orog` best in all six domains (best sub-cell shift = 0; RMSE 140–213 m). In the four domains where `rlat`/`rlon` were tested directly (AFR, EAS, EUR, NAM), placing the data at `rlat`/`rlon` positions gives RMSE 206–260 m. The extracted files keep `rlat`/`rlon` as dimension values only.

## Overlapping domains and `is_primary`

Some mountains lie in two domains (mainly AFR-12/EUR-12 and AUS-12/EAS-12). Each domain keeps all of its cells. For statistics that combine domains, each GMBA Basic unit is assigned to one **primary domain**:
- the domain where the unit has the largest area,
- or, if two domains are within 1 % of each other, the domain where the unit lies farther from the boundary.

A cell has `is_primary = 1` when `gmba_id > 0` and the primary domain of its `gmba_id` unit is this
domain. Usage:
- **Global or multi-domain statistics:** use only `is_primary == 1`, so each mountain is counted once.
- **Single-domain analyses:** use all cells and ignore `is_primary`.
- **Cells covered by two domains:** these can be used to compare the same model in two domain setups.

GMBA mountains outside all six domains (Arctic islands, Svalbard, Novaya Zemlya, Chukotka, southern Mexico and Central America, Sumatra, among others) are listed with `primary_domain = NONE` in `gmba_unit_primary.csv`.
