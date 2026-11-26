"""Command line interface: mspipe IMAGE [IMAGE ...] --bands R=1,G=2,B=3,NIR=4 --out results/"""
from __future__ import annotations

import argparse
import sys

from .core import parse_bands
from .pipeline import run
from .report import write_report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="mspipe", description="Vegetation and water indices, classes, change and zones from multispectral imagery.")
    p.add_argument("images", nargs="+", help="Multispectral GeoTIFFs; dates are read from YYYYMMDD in the file names")
    p.add_argument("--bands", default="R=1,G=2,B=3,NIR=4", help="Band numbers, e.g. R=3,G=2,B=1,RE=4,NIR=5")
    p.add_argument("--scale", type=float, default=255.0, help="Value that maps to 1.0 (255 for 8 bit, 1 for reflectance)")
    p.add_argument("--out", default="ms_results", help="Output folder")
    p.add_argument("--zone", type=float, default=200.0, help="Management zone size in metres")
    p.add_argument("--change", nargs=2, metavar=("FROM", "TO"), help="Dates (YYYYMMDD) to compare; default first and last")
    p.add_argument("--block", type=int, default=512, help="Processing block size in pixels")
    p.add_argument("--water", choices=["nir", "ndwi"], default="nir",
                   help="Open water rule: 'nir' (scene threshold, for uncalibrated imagery) or 'ndwi' (for reflectance)")
    p.add_argument("--title", default="Multispectral pipeline report")
    p.add_argument("--notes", default="", help="Text added under the report title")
    a = p.parse_args(argv)

    bands = parse_bands(a.bands)
    result = run(a.images, bands, a.out, a.scale, a.zone, tuple(a.change) if a.change else None, a.block, a.water)
    write_report(result, bands, a.out, a.title, a.notes)
    for r in result["summary"]["series"]:
        print(f"{r['label']}: median NDVI {r['median_ndvi_land']:.2f}, open water {r['water_ha']} ha, dense vegetation {r['dense_veg_pct_of_land']:.1f}% of land")
    c = result["summary"]["change"]
    if c:
        print(f"Change {c['from']} to {c['to']}: mean {c['mean_change']:+.3f}, gain {c['gain_ha']} ha, loss {c['loss_ha']} ha")
    print(f"Report written to {a.out}/report.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
