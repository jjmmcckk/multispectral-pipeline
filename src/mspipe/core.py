"""Multispectral processing: band mapping, grid alignment, windowed index calculation, classes and zones."""
from __future__ import annotations

import os
import re
import tempfile
from dataclasses import dataclass, field

import numpy as np
import rasterio
import rasterio.shutil
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

# Vegetation index formulas. Inputs are float arrays scaled to 0 to 1.
INDICES = {
    "NDVI": (("R", "NIR"), lambda b: _nd(b["NIR"], b["R"])),
    "GNDVI": (("G", "NIR"), lambda b: _nd(b["NIR"], b["G"])),
    "NDWI": (("G", "NIR"), lambda b: _nd(b["G"], b["NIR"])),          # McFeeters: open water above 0
    "NDRE": (("RE", "NIR"), lambda b: _nd(b["NIR"], b["RE"])),        # needs a red edge band
    "SAVI": (("R", "NIR"), lambda b: 1.5 * (b["NIR"] - b["R"]) / np.maximum(b["NIR"] + b["R"] + 0.5, 1e-6)),
}

# Classes from NDVI, with open water found first (see water_mask)
CLASSES = {1: "Open water", 2: "Bare or senescent", 3: "Sparse vegetation", 4: "Moderate vegetation", 5: "Dense vegetation"}
NDVI_BREAKS = (0.2, 0.4, 0.6)       # sparse from 0.2, moderate from 0.4, dense from 0.6


def _nd(a, b):
    return np.where((a + b) > 0, (a - b) / np.maximum(a + b, 1e-6), np.nan)


def parse_bands(spec: str) -> dict[str, int]:
    """'R=1,G=2,B=3,NIR=4' -> {'R': 1, ...} (1 based band numbers)."""
    out = {}
    for part in spec.split(","):
        name, num = part.split("=")
        out[name.strip().upper()] = int(num)
    if "NIR" not in out or "R" not in out:
        raise ValueError("band map needs at least R and NIR")
    return out


def date_from_name(path: str) -> str:
    m = re.search(r"(19|20)\d{6}", os.path.basename(path))
    return m.group(0) if m else os.path.splitext(os.path.basename(path))[0]


def otsu_threshold(values: np.ndarray, bins: int = 256) -> float:
    """Threshold that best splits values into two groups (Otsu's method)."""
    v = values[np.isfinite(values)]
    hist, edges = np.histogram(v, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers) / np.maximum(w0, 1)
    m1 = (np.sum(hist * centers) - np.cumsum(hist * centers)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return float(centers[np.argmax(between)])


def water_mask(b: dict, ndvi: np.ndarray, method: str, nir_threshold: float | None) -> np.ndarray:
    """Open water. 'ndwi' (green above NIR) suits calibrated reflectance. 'nir' suits uncalibrated
    8 bit imagery, where a fixed NDWI threshold also catches dry soil: water is where NIR is darker
    than a threshold found for each scene (Otsu) and NDVI is low."""
    if method == "ndwi":
        return _nd(b["G"], b["NIR"]) > 0
    return (b["NIR"] < nir_threshold) & (ndvi < 0.1)


def classify(ndvi: np.ndarray, water: np.ndarray | None) -> np.ndarray:
    cls = np.zeros(ndvi.shape, dtype=np.uint8)
    valid = ~np.isnan(ndvi)
    cls[valid] = 2
    cls[valid & (ndvi >= NDVI_BREAKS[0])] = 3
    cls[valid & (ndvi >= NDVI_BREAKS[1])] = 4
    cls[valid & (ndvi >= NDVI_BREAKS[2])] = 5
    if water is not None:
        cls[valid & water] = 1
    return cls


@dataclass
class Grid:
    crs: object
    transform: object
    width: int
    height: int

    @property
    def pixel_area_m2(self) -> float:
        return abs(self.transform.a * self.transform.e)


def reference_grid(path: str) -> Grid:
    with rasterio.open(path) as src:
        return Grid(src.crs, src.transform, src.width, src.height)


@dataclass
class DateResult:
    date: str
    path: str
    index_paths: dict[str, str]
    class_path: str
    class_area_ha: dict[str, float]
    ndvi_hist: np.ndarray
    stats: dict = field(default_factory=dict)


def _write_cog(tmp_path: str, out_path: str):
    with rasterio.open(tmp_path) as src:
        rasterio.shutil.copy(src, out_path, driver="COG", COMPRESS="DEFLATE", OVERVIEWS="AUTO", BLOCKSIZE=256)


def process_image(path: str, bands: dict[str, int], grid: Grid, out_dir: str, scale: float = 255.0,
                  indices=("NDVI", "GNDVI", "NDWI", "NDRE", "SAVI"), block: int = 512, water: str = "nir") -> DateResult:
    """Align one image to the reference grid and compute indices block by block, so memory use stays
    flat however large the image is. Writes one float32 COG per index and a class COG. Indices whose
    bands are missing from the band map (NDRE without a red edge band) are skipped."""
    date = date_from_name(path)
    indices = [i for i in indices if all(b in bands for b in INDICES[i][0])]
    os.makedirs(out_dir, exist_ok=True)
    tmpdir = tempfile.mkdtemp()
    profile = dict(driver="GTiff", width=grid.width, height=grid.height, crs=grid.crs, transform=grid.transform,
                   tiled=True, blockxsize=256, blockysize=256, compress="deflate")
    idx_tmp = {i: os.path.join(tmpdir, f"{i}.tif") for i in indices}
    cls_tmp = os.path.join(tmpdir, "classes.tif")
    writers = {i: rasterio.open(p, "w", count=1, dtype="float32", nodata=np.nan, **profile) for i, p in idx_tmp.items()}
    cls_w = rasterio.open(cls_tmp, "w", count=1, dtype="uint8", nodata=0, **profile)
    hist = np.zeros(200, dtype=np.int64)
    counts = np.zeros(len(CLASSES) + 1, dtype=np.int64)
    ndvi_vals = []

    with rasterio.open(path) as src, WarpedVRT(src, crs=grid.crs, transform=grid.transform, width=grid.width,
                                               height=grid.height, resampling=Resampling.bilinear) as vrt:
        nodata = src.nodata
        nir_threshold = None
        if water == "nir":       # scene threshold from a reduced resolution read of the NIR band
            f = max(1, int(np.ceil(max(grid.width, grid.height) / 1024)))
            nir = vrt.read(bands["NIR"], out_shape=(grid.height // f, grid.width // f)).astype("float32")
            nir_threshold = otsu_threshold(nir[nir > 0] / scale)
        for row in range(0, grid.height, block):
            for col in range(0, grid.width, block):
                w = Window(col, row, min(block, grid.width - col), min(block, grid.height - row))
                raw = {name: vrt.read(num, window=w).astype("float32") for name, num in bands.items()}
                mask = np.zeros(raw["NIR"].shape, dtype=bool)
                for a in raw.values():
                    mask |= (a == nodata) if nodata is not None else (a == 0)
                b = {k: v / scale for k, v in raw.items()}
                out = {}
                for i in indices:
                    val = INDICES[i][1](b).astype("float32")
                    val[mask] = np.nan
                    out[i] = val
                    writers[i].write(val, 1, window=w)
                ndvi = out["NDVI"]
                cls = classify(ndvi, water_mask(b, ndvi, water, nir_threshold))
                cls[mask] = 0
                cls_w.write(cls, 1, window=w)
                hist += np.histogram(ndvi[~np.isnan(ndvi)], bins=200, range=(-1, 1))[0]
                counts += np.bincount(cls.ravel(), minlength=len(CLASSES) + 1)
                land = (cls >= 2)
                if land.any():
                    ndvi_vals.append(ndvi[land][:: max(1, land.sum() // 20000)])
    for wtr in writers.values():
        wtr.close()
    cls_w.close()

    index_paths = {}
    for i, p in idx_tmp.items():
        index_paths[i] = os.path.join(out_dir, f"{date}_{i.lower()}.tif")
        _write_cog(p, index_paths[i])
    class_path = os.path.join(out_dir, f"{date}_classes.tif")
    _write_cog(cls_tmp, class_path)

    ha = grid.pixel_area_m2 / 10000
    area = {CLASSES[k]: round(float(counts[k] * ha), 2) for k in CLASSES}
    land_ndvi = np.concatenate(ndvi_vals) if ndvi_vals else np.array([np.nan])
    land_px = counts[2:].sum()
    stats = {"median_ndvi_land": float(np.nanmedian(land_ndvi)),
             "p90_ndvi_land": float(np.nanpercentile(land_ndvi, 90)),
             "water_ha": area["Open water"],
             "dense_veg_pct_of_land": float(100 * counts[5] / land_px) if land_px else 0.0,
             "valid_ha": round(float(counts[1:].sum() * ha), 2),
             "water_method": water, "nir_threshold": round(nir_threshold, 3) if nir_threshold is not None else None}
    return DateResult(date, path, index_paths, class_path, area, hist, stats)


def zone_stats(ndvi_path: str, class_path: str, zone_m: float) -> list[dict]:
    """Square management zones over the grid, with NDVI statistics and the share of each class."""
    with rasterio.open(ndvi_path) as n, rasterio.open(class_path) as c:
        ndvi, cls, t = n.read(1), c.read(1), n.transform
    px = int(round(zone_m / abs(t.a)))
    rows = []
    for r0 in range(0, ndvi.shape[0] - px + 1, px):
        for c0 in range(0, ndvi.shape[1] - px + 1, px):
            v = ndvi[r0:r0 + px, c0:c0 + px]
            k = cls[r0:r0 + px, c0:c0 + px]
            land = (k >= 2) & ~np.isnan(v)
            x0, y0 = t.c + c0 * t.a, t.f + r0 * t.e          # north up grid
            x1, y1 = t.c + (c0 + px) * t.a, t.f + (r0 + px) * t.e
            row = {"zone": f"r{r0 // px:02d}c{c0 // px:02d}", "bounds": (x0, y1, x1, y0),
                   "water_pct": round(100 * float((k == 1).mean()), 1),
                   "median_ndvi": round(float(np.median(v[land])), 3) if land.any() else None,
                   "dense_pct": round(100 * float((k == 5).mean()), 1)}
            rows.append(row)
    return rows
