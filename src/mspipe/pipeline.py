"""Run the pipeline over one or more dates: align, compute indices, summarise, compare dates, build zones."""
from __future__ import annotations

import csv
import json
import os
from datetime import datetime

import numpy as np
import pyproj
import rasterio
import rasterio.shutil
from shapely.geometry import Polygon, mapping

from .core import CLASSES, DateResult, Grid, date_from_name, process_image, reference_grid, zone_stats


def _pretty(date: str) -> str:
    try:
        return datetime.strptime(date, "%Y%m%d").strftime("%d %b %Y")
    except ValueError:
        return date


def ndvi_change(a: DateResult, b: DateResult, out_dir: str, threshold: float = 0.1) -> dict:
    """NDVI difference (later minus earlier) on land pixels of both dates, written as a COG."""
    with rasterio.open(a.index_paths["NDVI"]) as ra, rasterio.open(b.index_paths["NDVI"]) as rb, \
            rasterio.open(a.class_path) as ca, rasterio.open(b.class_path) as cb:
        na, nb = ra.read(1), rb.read(1)
        land = (ca.read(1) >= 2) & (cb.read(1) >= 2)
        profile = ra.profile
        px_ha = abs(ra.transform.a * ra.transform.e) / 10000
    diff = np.where(land, nb - na, np.nan).astype("float32")
    path = os.path.join(out_dir, f"change_ndvi_{a.date}_{b.date}.tif")
    tmp = path + ".tmp.tif"
    profile.update(driver="GTiff", dtype="float32", nodata=np.nan, tiled=True, blockxsize=256, blockysize=256, compress="deflate")
    with rasterio.open(tmp, "w", **profile) as w:
        w.write(diff, 1)
    with rasterio.open(tmp) as src:
        rasterio.shutil.copy(src, path, driver="COG", COMPRESS="DEFLATE", OVERVIEWS="AUTO")
    os.remove(tmp)
    d = diff[~np.isnan(diff)]
    return {"from": a.date, "to": b.date, "path": path, "mean_change": float(d.mean()) if d.size else None,
            "gain_ha": round(float((d > threshold).sum() * px_ha), 1), "loss_ha": round(float((d < -threshold).sum() * px_ha), 1),
            "threshold": threshold}


def zones_geojson(rows: list[dict], grid: Grid, path: str):
    to_ll = pyproj.Transformer.from_crs(grid.crs, 4326, always_xy=True)
    feats = []
    for r in rows:
        x0, y0, x1, y1 = r["bounds"]
        lon, lat = to_ll.transform([x0, x1, x1, x0], [y0, y0, y1, y1])
        geom = Polygon([(round(a, 7), round(b, 7)) for a, b in zip(lon, lat)])
        props = {k: v for k, v in r.items() if k != "bounds"}
        feats.append({"type": "Feature", "properties": props, "geometry": mapping(geom)})
    with open(path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)


def run(images: list[str], bands: dict[str, int], out_dir: str, scale: float = 255.0, zone_m: float = 200.0,
        change: tuple[str, str] | None = None, block: int = 512, water: str = "nir") -> dict:
    images = sorted(images, key=date_from_name)
    os.makedirs(out_dir, exist_ok=True)
    rasters = os.path.join(out_dir, "rasters")
    grid = reference_grid(images[0])
    results = [process_image(p, bands, grid, rasters, scale, block=block, water=water) for p in images]

    series = []
    for r in results:
        row = {"date": r.date, "label": _pretty(r.date), **{k: round(v, 3) if isinstance(v, float) else v for k, v in r.stats.items()}}
        row.update({f"{name} (ha)": ha for name, ha in r.class_area_ha.items()})
        series.append(row)
    with open(os.path.join(out_dir, "timeseries.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(series[0]))
        w.writeheader()
        w.writerows(series)

    by_date = {r.date: r for r in results}
    chg = None
    if len(results) > 1:
        a, b = (by_date[change[0]], by_date[change[1]]) if change else (results[0], results[-1])
        chg = ndvi_change(a, b, rasters)

    latest = results[-1]
    zones = zone_stats(latest.index_paths["NDVI"], latest.class_path, zone_m)
    zones_geojson(zones, grid, os.path.join(out_dir, "zones.geojson"))

    summary = {
        "images": [os.path.basename(p) for p in images], "bands": bands, "scale": scale, "water_method": water,
        "grid": {"crs": grid.crs.to_string(), "pixel_m": abs(grid.transform.a), "width": grid.width, "height": grid.height,
                 "reference": os.path.basename(images[0])},
        "series": series, "change": {k: v for k, v in (chg or {}).items() if k != "path"} or None,
        "zones": {"size_m": zone_m, "count": len(zones), "date": latest.date},
        "classes": CLASSES,
    }
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    return {"summary": summary, "results": results, "grid": grid, "change": chg, "zones": zones}
