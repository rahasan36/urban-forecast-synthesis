## Analysis

The folder analysis holds the code that turns the source forecast maps into the numbers of the paper and its Supplementary Information. The scripts read the rasters in data/rasters, in the layout of the Zenodo deposit, and write their results to analysis/output/<folder>. Both folders are ignored by git. The environment variables URBAN_FORECAST_RASTERS and URBAN_FORECAST_OUTPUT point the scripts to other copies. analysis/common.py holds the forecasts, scenarios, years, regions, land-cover classes and grid shared by all scripts, and analysis/landcover_class_codes.json groups the ESA CCI classes into the 12 nonurban classes.

| Folder | Contents |
|---|---|
| harmonize | Projection of the inputs onto the 1-km equal-area grid (EPSG:6933). harmonize_forecasts.py reprojects the 81 source maps listed in source_index.csv (forecast, scenario, year, RCP, file under data/rasters/source, SHA-256), stores fractional forecasts as threshold ranks and builds the common valid domain (common_valid.tif, grid.json). landcover_baseline.py projects the ESA CCI / C3S land-cover maps of 2016 and 2020, region_raster.py rasterizes the 20 AR6 reporting regions and country_raster.py the Natural Earth countries. |
| measure | measure_maps.py: urban area, overlap with the 2016 land cover and boundary contact by region and class of every map, with WGS84 ellipsoidal side lengths (Tables S2, S4 to S7). configuration.py: patch count, mean patch area and edge density (Table S3, Fig. S2). summaries.py: changes from 2020 to 2100 and ensemble contact means, including the FLUS weighting (Tables S4, S6, S8). zone_contact.py: boundary and cropland and forest contact of each forecast by country (Fig. 3) or region (Table S9, Fig. S6). blocks.py: area and boundary of each forecast in 100-km blocks (Fig. 2, Figs. S3 and S4). checks.py: accounting identities and gross and net changes. |
| agreement | ensemble.py: nine-forecast membership rasters, union and unanimous areas, pair intersections, quantity and allocation, leave-one-out ensembles, inherited and new disagreement, shared boundary segments and 50-km block sums (Tables S12 and S13, Figs. S7 and S13 to S15). intersections_latitude.py: exact intersections, persistence and latitude profiles (Figs. S8 and S9). agreement_stats.py and test_agreement_stats.py: the pooled pairwise Jaccard index and the amount and placement statistics (Eqs. S6 to S10) with their synthetic checks. aggregate.py and zone_scale.py: pairwise tallies by 20-km block, country and region and the zone sums (Fig. 4). amount_placement.py: the amount and placement split by zone size, distance and region. pairwise_table.py: Table S12 and its note. concordance.py: Kendall's W of regional ranks (Table S9). |
| clustering | block_sums.py: urban counts and union boundary in 50-km blocks. block_regions.py: region of each block. gistar.py: the Gi* statistic (Eq. S12) and its synthetic check. primary_clusters.py: disagreement intensity (Eq. S11) and Gi* over a circular neighborhood of radius four blocks for blocks with at least 10 km2 of union urban land. strict_clusters.py: block support and the stricter 9 x 9 specification with at least 25 km2, at least nine neighbors and Benjamini-Yekutieli control (Eq. S13). region_summary.py: regional shares under both specifications and the Spearman correlation of their rankings (Fig. S12). |
| benchmarks | observed_2020.py and builtup_2020.py: intersections of each forecast's 2020 map with observed urban land (class 190 of the 2020 land cover) and with the GHS-BUILT-S 2020 built-up fraction. benchmark_tables.py: Tables S10 and S11. |
| area_reconciliation | full_domain_areas.py: urban area of every source map and base-year map on its full source domain and on the common domain, on a global 1-km EPSG:6933 grid. gao_land_weighted.py: the land-weighted area of Gao 2021. area_comparison.py: the comparison of full-source and common-domain areas (Methods; the full comparison is in the deposit). |

Run the scripts from the repository folder in this order. The harmonization needs the source maps of the nine forecasts, which are distributed by their authors (Table S1) and are not redistributed here; all later steps start from the harmonized rasters on Zenodo.

    python analysis/harmonize/landcover_baseline.py 2016
    python analysis/harmonize/landcover_baseline.py 2020
    python analysis/harmonize/region_raster.py
    python analysis/harmonize/country_raster.py
    python analysis/harmonize/harmonize_forecasts.py
    python analysis/measure/measure_maps.py
    python analysis/measure/configuration.py
    python analysis/measure/summaries.py
    python analysis/agreement/ensemble.py
    python analysis/measure/checks.py
    python analysis/agreement/intersections_latitude.py
    python analysis/measure/zone_contact.py country SSP1 2100
    python analysis/measure/zone_contact.py region SSP1 2100      (and the other eight scenario-year pairs)
    python analysis/measure/blocks.py SSP1 2100                   (and the other eight scenario-year pairs)
    python analysis/agreement/aggregate.py
    python analysis/agreement/zone_scale.py
    python analysis/agreement/amount_placement.py
    python analysis/agreement/pairwise_table.py
    python analysis/agreement/concordance.py
    python analysis/clustering/block_sums.py SSP1 2020
    python analysis/clustering/block_sums.py SSP1 2100
    python analysis/clustering/block_regions.py
    python analysis/clustering/primary_clusters.py
    python analysis/clustering/strict_clusters.py
    python analysis/clustering/region_summary.py
    python analysis/benchmarks/observed_2020.py
    python analysis/benchmarks/builtup_2020.py
    python analysis/benchmarks/benchmark_tables.py
    python analysis/area_reconciliation/full_domain_areas.py
    python analysis/area_reconciliation/area_comparison.py
    python analysis/agreement/test_agreement_stats.py
    python analysis/clustering/gistar.py

The global grid has 34,736 x 13,367 cells. measure_maps.py, ensemble.py and configuration.py read whole rasters and need about 4 GB of memory or more; aggregate.py and zone_scale.py take --part 0 and --part 1 to process half of the rows at a time, and --combine to join the parts.

## Zenodo deposit

The scripts expect the deposit unzipped into data/rasters:

| Path in data/rasters | Contents | Size |
|---|---|---|
| harmonized/<forecast>_SSP<s>_<year>.tif and .json | The 81 harmonized maps (forecasts Li2021, Chen2020, Chen2022, Gao2021, Hou2022, Koomen2UP (2UP 2023), He2023, Zhang2023 and Zhao2025 (FU3D 2025); SSP1, SSP3, SSP5; 2020, 2050, 2100), uint8 0/1/255, or threshold ranks 0 to 3 for Gao2021 and He2023, with the value histogram of each map | 255 MB |
| harmonized/ensemble_SSP<s>_<year>.tif | The nine membership rasters, uint16 bitmask (bit 0 Li 2021 to bit 8 FU3D 2025), 65535 invalid | 77 MB |
| common_valid.tif, grid.json | Common valid domain (130,444,801 cells) and grid description | 3 MB |
| baseline/ESA_CCI_2016_lccs_cea1km.tif | 2016 land cover on the analysis grid | 45 MB |
| baseline/ESA_CCI_2020_lccs_cea1km.tif | 2020 land cover on the analysis grid | 34 MB |
| baseline/AR6_20regions_cea1km.tif, AR6_20regions_codes.json | Reporting regions | 3 MB |
| baseline/countries_cea1km.tif, country_codes.json | Natural Earth countries | 4 MB |

The source maps are read from data/rasters/source/<forecast>/... at the paths in harmonize/source_index.csv and area_reconciliation/baseline_sources.csv, the land-cover netCDF files from data/rasters/source/ESA_CCI, the AR6 regions from data/rasters/source/AR6_regions and the built-up layer from data/rasters/source/GHSL. These third-party files keep their own licenses and are downloaded from their publishers.
