# Input data (not stored on GitHub)

These two datasets are needed only to **rebuild** the masks (`scripts/03_build_masks.py`). They are too large
for GitHub (1.6 GB and 330 MB), so download them into this folder with the names below.

```
DEM_and_shp/
├── ETOPO_2022_v1_30s_N90W180_surface.nc
└── GMBA_Inventory_v2.0_standard/
    ├── GMBA_Inventory_v2.0_standard.shp
    ├── GMBA_Inventory_v2.0_standard.shx
    ├── GMBA_Inventory_v2.0_standard.dbf
    ├── GMBA_Inventory_v2.0_standard.prj
    └── GMBA_Inventory_v2.0_standard.cpg
```

## GMBA Mountain Inventory v2.0, standard

Download "GMBA Inventory v2.0 standard" from https://www.earthenv.org/mountains and unzip it here.

Snethlage, M. A., et al. (2022). A hierarchical inventory of the world's mountains for global comparative
mountain science. *Scientific Data* 9, 149. https://doi.org/10.1038/s41597-022-01256-y

## ETOPO 2022 v1, 30 arc-second, surface elevation (netCDF)

```bash
wget https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc
```

Product page: https://www.ncei.noaa.gov/products/etopo-global-relief-model (https://doi.org/10.25921/fd45-gt74).

The paths can also be set with the environment variables `MTN_DEM` and `MTN_SHP`.
