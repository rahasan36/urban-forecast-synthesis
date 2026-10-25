# Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

This repository holds the code and data for the paper above, which compares nine global urban land forecasts under SSP1, SSP3 and SSP5 for 2020, 2050 and 2100 on a common 1-km equal-area grid. The paper measures how far the forecasts agree on how much land becomes urban, where that land lies and the boundaries it forms with cropland and forest. The repository includes the Python code of the analysis, the scripts that draw the figures of the paper and its Supplementary Information with the figures they draw, the inputs of the figure scripts, and the tables of the data deposit.

The outputs of the analysis are on Zenodo (DOI to be added at deposit). Neither this repository nor the Zenodo deposit redistributes the forecast maps, which are available from their authors.

## Forecasts

The file names and tables use one key for each forecast. All nine are distributed by their authors at the links below, are described in the papers listed under References, and are not redistributed here.

| Forecast in the paper | Key in files | Data |
|---|---|---|
| Li 2021<sup>1</sup> | Li2021 | https://doi.org/10.6084/m9.figshare.9696218 |
| Chen 2020<sup>2</sup> | Chen2020 | https://doi.org/10.1594/PANGAEA.905890 |
| Chen 2022<sup>3</sup> | Chen2022 | https://doi.org/10.5281/zenodo.4584775 |
| Gao 2021<sup>4</sup> | Gao2021 | https://doi.org/10.7927/1z4r-ez63 |
| Hou 2022<sup>5</sup> | Hou2022 | https://doi.org/10.6084/m9.figshare.20088368 |
| 2UP 2023<sup>6</sup> | Koomen2UP | https://doi.org/10.5281/zenodo.15396282 |
| He 2023<sup>7</sup> | He2023 | https://doi.org/10.6084/m9.figshare.20391117.v4 |
| Zhang 2023<sup>8</sup> | Zhang2023 | https://doi.org/10.6084/m9.figshare.23542860 |
| FU3D 2025<sup>9</sup> | Zhao2025 | https://doi.org/10.6084/m9.figshare.26795932.v2 |

Gao 2021 and He 2023 map urban fractions, which enter the analysis as urban at a fraction of 0.40, with 0.30 and 0.50 as alternatives.

## Figures

The folder figures holds the script behind each figure of the paper and its Supplementary Information, one folder per figure, with the figure it draws as PNG (300 dpi) and PDF. The folder data/figures holds the inputs the scripts read. Run a script from the repository folder without arguments, for example

    python figures/Fig2_block_dispersion/make_fig2_block_dispersion.py

Each script writes the figure to its own folder, or to the folder given by --out. The input arguments, listed by --help, point a script to other copies of the inputs. The scripts do not overwrite an existing figure and add _r1, _r2 and so on to the name of a new one. The maps read the Natural Earth layers in data/figures/naturalearth, so the scripts need no network access.

| Figure | Content | Script (in figures/) | Inputs (in data/figures) |
|---|---|---|---|
| Figure 1 | Variation among forecasts in urban boundary length relative to urban area | Fig1_boundary_dispersion/make_fig1_boundary_dispersion.py | results/global_boundary_area.csv |
| Figure 2 | Block-scale variation in boundary length relative to urban area | Fig2_block_dispersion/make_fig2_block_dispersion.py | blocks100km, naturalearth |
| Figure 3 | Cropland and forest shares of the projected urban boundary in 2100, SSP1 | Fig3_boundary_contact/make_fig3_boundary_contact.py | countries, regions, boundary_composition, naturalearth |
| Figure 4 | Pairwise agreement among forecasts on the location of urban land, SSP1 | Fig4_pairwise_agreement/make_fig4_pairwise_agreement.py | pairwise_agreement, countries, naturalearth |
| Supplementary Fig. 1 | Projected urban expansion and overlap with other land covers | SFig01_landcover_conversion/make_sfig01_landcover_conversion.py | results/global_metrics.csv, results/landcover_endpoint_changes.csv |
| Supplementary Fig. 2 | Urban landscape configuration across nine forecasts | SFig02_configuration/make_sfig02_configuration.py | results/configuration_metrics.csv |
| Supplementary Fig. 3 | Boundary density of each forecast relative to the block median | SFig03_block_signatures/make_sfig03_block_signatures.py | blocks100km |
| Supplementary Fig. 4 | Block size, date, and the dispersion of boundary length and urban area | SFig04_block_threshold/make_sfig04_block_threshold.py | blocks100km |
| Supplementary Fig. 5 | Urban boundaries adjoining cropland and forest across scenarios | SFig05_contact_over_time/make_sfig05_contact_over_time.py | results/landcover_metrics.csv |
| Supplementary Fig. 6 | Regional change in cropland and forest contact and its spread among forecasts | SFig06_regional_contact_growth/make_sfig06_regional_contact_growth.py | results/landcover_endpoint_changes.csv |
| Supplementary Fig. 7 | Agreement and disagreement on urban location in 2100, SSP1 | SFig07_agreement_origin/make_sfig07_agreement_origin.py | agreement_blocks, regions |
| Supplementary Fig. 8 | Exact forecast intersections in the nine-forecast ensemble, SSP1 | SFig08_intersections/make_sfig08_intersections.py | intersections |
| Supplementary Fig. 9 | Latitudinal distribution of unanimous urban agreement and forecast variation | SFig09_latitude/make_sfig09_latitude.py | latitude |
| Supplementary Fig. 10 | Urban union in 2020, 2050 and 2100, SSP1, SSP3 and SSP5 | SFig10_urban_union_hexagons/make_sfig10_urban_union_hexagons.py | blocks5km |
| Supplementary Fig. 11 | Timing of disagreement among the nine forecasts, SSP1, SSP3 and SSP5 | SFig11_disagreement_timing/make_sfig11_disagreement_timing.py | blocks5km, regions, naturalearth |
| Supplementary Fig. 12 | Local clustering of urban support under two specifications, SSP1, 2100 | SFig12_support_clusters/make_sfig12_support_clusters.py | support_clusters, naturalearth |
| Supplementary Fig. 13 | Urban membership in six metropolitan windows, SSP1, 2020 and 2100 | SFig13_metropolitan_windows/make_sfig13_metropolitan_windows.py | metropolitan |
| Supplementary Fig. 14 | Urban membership in seven large metropolitan regions, SSP1, 2020 and 2100 | SFig14_metropolitan_established/make_sfig14_metropolitan_established.py | metropolitan |
| Supplementary Fig. 15 | Urban membership in six lower-income urbanizing regions, SSP1, 2020 and 2100 | SFig15_metropolitan_growth/make_sfig15_metropolitan_growth.py | metropolitan |

The same list is in figures/figure_scripts.csv.

## Figure inputs

| Folder in data/figures | Contents |
|---|---|
| results | Measurements of every forecast map, written by the measure scripts: global urban area and boundary length (global_metrics.csv, global_boundary_area.csv), overlap and boundary contact with each land-cover class by region (landcover_metrics.csv), their change from 2020 to 2100 (landcover_endpoint_changes.csv) and patch area and edge density (configuration_metrics.csv). |
| blocks100km | Urban area and boundary length of each forecast in 100 × 100 km blocks for every pathway and year, and the region of each block. |
| countries | Natural Earth countries on the analysis grid with their codes, and the cropland and forest contact of each forecast by country in 2100 under SSP1. |
| regions | Display outlines of the 19 reporting regions and country outlines for the maps. |
| boundary_composition | Median cropland and forest shares of the urban boundary by region, SSP1, 2100. |
| pairwise_agreement | Counts of forecasts classifying each cell as urban, tallied in 20-km blocks for 2100 and globally for 2020, with the pair intersections and urban areas by region and country, written by analysis/figure_inputs/prepare_pairwise_agreement_inputs.py. |
| agreement_blocks | Union, unanimous, majority and single-forecast urban land and inherited and new disagreement in 50-km blocks, SSP1, 2100, with the grid description. |
| intersections | Area of each of the 511 nonempty forecast combinations, SSP1. |
| latitude | Unanimous urban land and the share of each forecast in the band union, in 2° latitude bands. |
| blocks5km | Union urban land and the year disagreement first appears, in 5-km blocks, for the three pathways, and the valid cells of each block. |
| support_clusters | Gi* z scores of disagreement intensity under the primary specification and the cluster classes of the stricter specification in 50-km blocks, with the cropland boundary of each block and the regional summary. |
| metropolitan | Windows of the membership rasters around the metropolitan regions of Supplementary Figs 13 to 15, 2020 and 2100, SSP1, cut by analysis/figure_inputs/extract_metropolitan_windows.py. |
| naturalearth | Natural Earth 1:110m land, coastline and land boundaries and 1:50m admin-0 countries. |

## Analysis

The folder analysis holds the code that turns the forecast maps into the numbers of the paper and its Supplementary Information. The scripts read the rasters in data/rasters (not provided, available from the data sources listed) and write their results to analysis/output/<folder>. The environment variables URBAN_FORECAST_RASTERS and URBAN_FORECAST_OUTPUT point the scripts to other copies. analysis/common.py holds the forecasts, scenarios, years, regions, land-cover classes and grid shared by all scripts, and analysis/landcover_class_codes.json groups the land-cover classes into the 12 nonurban classes.

| Folder | Contents |
|---|---|
| harmonize | Projection of the inputs onto the 1-km equal-area grid (EPSG:6933). harmonize_forecasts.py reprojects the 81 source maps listed in source_index.csv (forecast, scenario, year, RCP, file under data/rasters/source and SHA-256 checksum), stores fractional forecasts as threshold ranks and builds the common valid domain (common_valid.tif, grid.json). landcover_baseline.py projects the land-cover maps of 2016 and 2020, region_raster.py rasterizes the 20 reporting regions and country_raster.py the Natural Earth countries. |
| measure | measure_maps.py: urban area, overlap with the 2016 land cover and boundary contact by region and class of every map, with WGS84 ellipsoidal side lengths (Tables S2 and S4 to S7). configuration.py: patch count, mean patch area and edge density (Table S3, Supplementary Fig. 2). summaries.py: changes from 2020 to 2100 and ensemble contact means, including the weighting of the three FLUS forecasts (Tables S4, S6 and S8). zone_contact.py: boundary and cropland and forest contact of each forecast by country (Figure 3) or region (Table S9, Supplementary Fig. 6). blocks.py: area and boundary of each forecast in 100-km blocks (Figure 2, Supplementary Figs 3 and 4). checks.py: accounting identities and gross and net changes. |
| agreement | ensemble.py: nine-forecast membership rasters, union and unanimous areas, pair intersections, quantity and allocation, leave-one-out ensembles, inherited and new disagreement and 50-km block sums (Tables S12 and S13, Supplementary Figs 7 and 13 to 15). intersections_latitude.py: exact intersections and latitude profiles (Supplementary Figs 8 and 9). agreement_stats.py and test_agreement_stats.py: the pooled pairwise Jaccard index and the amount and placement statistics (Eqs. S6 to S10) with their checks on synthetic grids. aggregate.py and zone_scale.py: pairwise tallies by 20-km block, country and region (Figure 4). amount_placement.py: the amount and placement split by zone. pairwise_table.py: Table S12 and its note. concordance.py: Kendall's coefficient of concordance of the regional ranks (Table S9). |
| clustering | block_sums.py: urban counts and union boundary in 50-km blocks. block_regions.py: region of each block. gistar.py: the Gi* statistic (Eq. S12) and its check on a synthetic grid. primary_clusters.py: disagreement intensity (Eq. S11) and Gi* over a circular neighborhood of radius four blocks, for blocks with at least 10 km² of union urban land. strict_clusters.py: the stricter specification over a 9 × 9 block neighborhood, with at least 25 km² of union urban land, at least nine neighbors and Benjamini-Yekutieli control (Eq. S13). region_summary.py: regional shares under both specifications and the rank correlation between them (Supplementary Fig. 12). |
| benchmarks | observed_2020.py and builtup_2020.py: correspondence of each forecast's 2020 map with observed urban land (class 190 of the 2020 land cover) and with the 2020 built-up surface fraction. benchmark_tables.py: Tables S10 and S11. |
| area_reconciliation | full_domain_areas.py: urban area of every analyzed and base-year map on its full source domain and on the common domain. gao_land_weighted.py: the land-weighted area of Gao 2021. area_comparison.py: the comparison of the two areas reported in the Methods. |
| figure_inputs | prepare_pairwise_agreement_inputs.py and extract_metropolitan_windows.py write the inputs in data/figures/pairwise_agreement and data/figures/metropolitan. |

deposit_tables.py writes the tables in data/tables from the outputs of area_comparison.py, summaries.py and intersections_latitude.py.

Run the scripts from the repository folder in the order below. The harmonization needs the source maps of the nine forecasts, which their authors distribute, and the land-cover, region and built-up surface layers listed in the Sources section at the end of this README. The outputs of every step are in the Zenodo deposit, so the table steps can also start from there.

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
    python analysis/measure/zone_contact.py region SSP1 2100      (and the other eight pathway and year pairs)
    python analysis/measure/blocks.py SSP1 2100                   (and the other eight pathway and year pairs)
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
    python analysis/deposit_tables.py
    python analysis/agreement/test_agreement_stats.py
    python analysis/clustering/gistar.py

The global grid has 34,736 × 13,367 cells. measure_maps.py, ensemble.py and configuration.py read whole rasters and need about 4 GB of memory or more. aggregate.py and zone_scale.py take --part 0 and --part 1 to process half of the rows at a time, and --combine to join the parts.

## Data tables

| File in data/tables | Contents |
|---|---|
| urban_area_by_domain.csv | Urban area of each of the 81 analyzed maps on its full source domain and on the common analysis domain, the area outside the common domain and its share, and the total reported by the authors where one was located. |
| urban_area_base_year_maps.csv | The same for the base-year map distributed with each forecast. FU3D 2025 distributes no base-year extent map. |
| regional_contact_change.csv | For each region and pathway, the ratio of 2100 to 2020 contact length with cropland and with forest (median, minimum and maximum across the nine forecasts, and the number of forecasts with growth), and the ratio of the longest to the shortest forecast contact length in 2100. |
| exact_intersections.csv | Area of each of the 511 nonempty forecast combinations for every pathway and year. |

The file names and checksums of the analyzed source maps are in analysis/harmonize/source_index.csv.

## Zenodo deposit

The Zenodo deposit holds analysis_output.zip, the outputs of the analysis in the layout the scripts write, and analysis_grid.zip, the common valid domain, the grid description and the region code table. Unzip analysis_output.zip into analysis/output and analysis_grid.zip into data/rasters. pairwise_table.py, concordance.py, region_summary.py, benchmark_tables.py, area_comparison.py, summaries.py, strict_clusters.py and deposit_tables.py then run without the forecast maps.

The folder data/rasters is not part of the repository. Download each forecast from the source given in Table S1 of the paper and in the source_url column of analysis/harmonize/source_index.csv, and place its files under data/rasters/source/<key> at the paths listed in source_index.csv (the analyzed maps) and in analysis/area_reconciliation/baseline_sources.csv (the base-year maps), for example

    data/rasters/source/Li2021/gUrban_2020_SSP1.tif
    data/rasters/source/Chen2020/Global urban product_SSPs_2015_2100/SSP1/global_SSP1_2020.tif

The SHA-256 checksums in source_index.csv identify the analyzed files. The land-cover files go to data/rasters/source/ESA_CCI, the region polygons to data/rasters/source/AR6_regions and the built-up surface layer to data/rasters/source/GHSL. The harmonization scripts then write the 1-km maps to data/rasters/harmonized and data/rasters/baseline.

## Environment

The scripts were run with Python 3.10, numpy 2.2, pandas 2.3, scipy 1.15, matplotlib 3.10, rasterio 1.4 (GDAL 3.10), pyproj 3.7, shapely 2.1, cartopy 0.25, pyogrio 0.13, openpyxl 3.1 and scikit-image. environment.yml creates a conda environment with these versions (`conda env create -f environment.yml`), and requirements.txt lists the packages for pip.

## Sources

Land cover from the Copernicus Climate Change Service, land cover classification gridded maps from 1992 to present derived from satellite observations, version 2.1.1 (https://doi.org/10.24381/cds.006f2c9a). Built-up surface from GHS-BUILT-S R2023A, European Commission, Joint Research Centre (Pesaresi and Politis 2023, https://doi.org/10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA). Reporting regions from the IPCC AR6 grouping of countries into 20 regions. Country boundaries and land from Natural Earth (https://www.naturalearthdata.com), which is in the public domain. The urban forecasts are cited in the paper and are available from their authors. The data in data/figures and on Zenodo keep the licenses of their sources.

## References

<sup>1</sup> Li, X., Zhou, Y., Hejazi, M., Wise, M., Vernon, C., Iyer, G., & Chen, W. (2021). Global urban growth between 1870 and 2100 from integrated high resolution mapped data and urban dynamic modeling. *Communications Earth & Environment, 2*(1), 201. https://doi.org/10.1038/s43247-021-00273-w

<sup>2</sup> Chen, G., Li, X., Liu, X., Chen, Y., Liang, X., Leng, J., . . . Huang, K. (2020). Global projections of future urban land expansion under shared socioeconomic pathways. *Nature Communications, 11*(1), 537. https://doi.org/10.1038/s41467-020-14386-x

<sup>3</sup> Chen, G., Li, X., & Liu, X. (2022). Global land projection based on plant functional types with a 1-km resolution under socio-climatic scenarios. *Scientific Data, 9*(1), 125. https://doi.org/10.1038/s41597-022-01208-6

<sup>4</sup> Gao, J., & Pesaresi, M. (2021). Downscaling SSP-consistent global spatial urban land projections from 1/8-degree to 1-km resolution 2000–2100. *Scientific Data, 8*(1), 281. https://doi.org/10.1038/s41597-021-01052-0

<sup>5</sup> Hou, H., Zhou, B.-B., Pei, F., Hu, G., Su, Z., Zeng, Y., . . . Li, X. (2022). Future land use/land cover change has nontrivial and potentially dominant impact on global gross primary productivity. *Earth's Future, 10*(9), e2021EF002628. https://doi.org/10.1029/2021EF002628

<sup>6</sup> Koomen, E., van Bemmel, M. S., van Huijstee, J., Andrée, B. P. J., Ferdinand, P. A., & van Rijn, F. J. A. (2023). An integrated global model of local urban development and population change. *Computers, Environment and Urban Systems, 100*, 101935. https://doi.org/10.1016/j.compenvurbsys.2022.101935

<sup>7</sup> He, W., Li, X., Zhou, Y., Shi, Z., Yu, G., Hu, T., . . . Gong, P. (2023). Global urban fractional changes at a 1 km resolution throughout 2100 under eight scenarios of Shared Socioeconomic Pathways (SSPs) and Representative Concentration Pathways (RCPs). *Earth System Science Data, 15*(8), 3623–3639. https://doi.org/10.5194/essd-15-3623-2023

<sup>8</sup> Zhang, T., Cheng, C., & Wu, X. (2023). Mapping the spatial heterogeneity of global land use and land cover from 2020 to 2100 at a 1 km resolution. *Scientific Data, 10*(1), 748. https://doi.org/10.1038/s41597-023-02637-7

<sup>9</sup> Zhao, Q., Meng, Q., Gao, L., & Zhu, M. (2025). FU3D: The first global projections of future urban three-dimensional (3D) expansion for the 21st century under shared socioeconomic pathways. *Scientific Data, 12*(1), 1555. https://doi.org/10.1038/s41597-025-05711-4
