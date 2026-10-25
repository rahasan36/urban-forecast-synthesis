"""Divergent urban boundaries and high spatial uncertainty in global urban land forecasts

Author: Rakibul Ahasan
"""
from pathlib import Path
import json
import os
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

REPO = Path(__file__).resolve().parents[1]
RASTERS = Path(os.environ.get("URBAN_FORECAST_RASTERS", REPO / "data" / "rasters"))
SOURCE = RASTERS / "source"
BASELINE = RASTERS / "baseline"
HARMONIZED = RASTERS / "harmonized"
OUTPUT = Path(os.environ.get("URBAN_FORECAST_OUTPUT", REPO / "analysis" / "output"))

PRODUCTS = ["Li2021", "Chen2020", "Chen2022", "Gao2021", "Hou2022", "Koomen2UP", "He2023", "Zhang2023", "Zhao2025"]
LABELS = ["Li 2021", "Chen 2020", "Chen 2022", "Gao 2021", "Hou 2022", "2UP 2023", "He 2023", "Zhang 2023", "FU3D 2025"]
LABEL = dict(zip(PRODUCTS, LABELS))
SSPS = [1, 3, 5]
YEARS = [2020, 2050, 2100]
FRACTIONAL = ["Gao2021", "He2023"]
URBAN_CODE = {"Li2021": 1, "Chen2020": 2, "Chen2022": 6, "Hou2022": 13, "Koomen2UP": 1, "Zhang2023": 4, "Zhao2025": 1}
THRESHOLD_RANK = {0.3: 1, 0.4: 2, 0.5: 3}
PRIMARY_THRESHOLD = 0.4

WIDTH, HEIGHT = 34736, 13367
X0, Y0 = -17367529.34, 7289708.89
CELL = 1000.0
EARTH_RADIUS = 6378137.0
COS30 = np.cos(np.radians(30.0))

REGIONS = {1: "Asia-Pacific Developed", 2: "Australia & New Zealand", 3: "Caribbean", 4: "Developing Pacific",
           5: "Eastern Africa", 6: "Eastern Asia", 7: "Eurasia", 8: "Greenland, Bermuda + others",
           9: "India & Sri Lanka", 10: "Meso America", 11: "Middle East", 12: "North Africa",
           13: "Northern and Western Europe", 14: "Southern Asia (rest)", 15: "South America",
           16: "South-East Asia", 17: "Southern and Eastern Europe", 18: "Southern and Middle Africa",
           19: "USA & Canada", 20: "Western Africa"}
REGIONS_MAIN = [i for i in REGIONS if i != 8]

CLASS_CODES = json.loads((REPO / "analysis" / "landcover_class_codes.json").read_text(encoding="utf-8-sig"))
CLASS_CODES = {("Flooded vegetation" if k == "Flooded forest" else k): v for k, v in CLASS_CODES.items()}
CLASSES = list(CLASS_CODES) + ["Baseline urban"]
URBAN_LANDCOVER_CODE = 190

LANDCOVER_2016 = BASELINE / "ESA_CCI_2016_lccs_cea1km.tif"
LANDCOVER_2020 = BASELINE / "ESA_CCI_2020_lccs_cea1km.tif"
REGION_RASTER = BASELINE / "AR6_20regions_cea1km.tif"
REGION_CODES = BASELINE / "AR6_20regions_codes.json"
COUNTRY_RASTER = BASELINE / "countries_cea1km.tif"
COUNTRY_CODES = BASELINE / "country_codes.json"
COMMON_VALID = RASTERS / "common_valid.tif"
GRID = RASTERS / "grid.json"


def forecast_path(product, ssp, year):
    return HARMONIZED / f"{product}_SSP{ssp}_{year}.tif"


def ensemble_path(ssp, year):
    return HARMONIZED / f"ensemble_SSP{ssp}_{year}.tif"


def urban_mask(a, product, threshold=PRIMARY_THRESHOLD):
    if product in FRACTIONAL:
        return (a >= THRESHOLD_RANK[threshold]) & (a != 255)
    return (a == 1) & (a != 255)


def urban_rank(product):
    return THRESHOLD_RANK[PRIMARY_THRESHOLD] if product in FRACTIONAL else 1


def output_dir(stage):
    d = OUTPUT / stage
    d.mkdir(parents=True, exist_ok=True)
    return d

