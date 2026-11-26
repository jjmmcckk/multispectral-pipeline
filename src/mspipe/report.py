"""HTML report: latest date maps, NDVI and water through time, change between two dates, zones."""
from __future__ import annotations

import base64
import html
import io
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.vrt import WarpedVRT  # noqa: E402

from .core import CLASSES  # noqa: E402

CLASS_COLORS = ["#ffffff", "#2b7bba", "#c8b28a", "#e3e38a", "#8cc65a", "#2e7d32"]


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _read_aligned(path, grid, bands, max_px=1200):
    f = max(1, int(np.ceil(max(grid.width, grid.height) / max_px)))
    shape = (grid.height // f, grid.width // f)
    with rasterio.open(path) as src, WarpedVRT(src, crs=grid.crs, transform=grid.transform, width=grid.width,
                                               height=grid.height, resampling=Resampling.bilinear) as vrt:
        return {k: vrt.read(v, out_shape=shape).astype(float) for k, v in bands.items()}


def _stretch(a):
    lo, hi = np.nanpercentile(a, [2, 98])
    return np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1)


def latest_maps(run: dict, bands: dict) -> str:
    r = run["results"][-1]
    b = _read_aligned(r.path, run["grid"], {k: bands[k] for k in ("R", "G", "B", "NIR") if k in bands})
    with rasterio.open(r.index_paths["NDVI"]) as n, rasterio.open(r.class_path) as c:
        ndvi = n.read(1, out_shape=b["R"].shape)
        cls = c.read(1, out_shape=b["R"].shape, resampling=Resampling.nearest)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.6))
    if "B" in b:
        ax[0, 0].imshow(np.dstack([_stretch(b["R"]), _stretch(b["G"]), _stretch(b["B"])]))
    ax[0, 0].set_title("True colour")
    ax[0, 1].imshow(np.dstack([_stretch(b["NIR"]), _stretch(b["R"]), _stretch(b["G"])]))
    ax[0, 1].set_title("False colour (NIR, red, green): vegetation shows red")
    im = ax[1, 0].imshow(ndvi, cmap="RdYlGn", vmin=-0.2, vmax=0.8)
    ax[1, 0].set_title("NDVI")
    fig.colorbar(im, ax=ax[1, 0], shrink=0.8)
    ax[1, 1].imshow(cls, cmap=ListedColormap(CLASS_COLORS), vmin=0, vmax=5, interpolation="nearest")
    ax[1, 1].set_title("Classes")
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[k]) for k in CLASSES]
    ax[1, 1].legend(handles, CLASSES.values(), fontsize=8, loc="lower left", framealpha=0.9)
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"Latest date: {run['summary']['series'][-1]['label']}")
    fig.tight_layout()
    return _png(fig)


def series_chart(run: dict) -> str:
    s = run["summary"]["series"]
    labels = [x["label"] for x in s]
    x = np.arange(len(s))
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(x, [v["median_ndvi_land"] for v in s], "o-", color="#2e7d32", label="Median NDVI, land")
    ax.set_ylabel("NDVI")
    ax2 = ax.twinx()
    ax2.bar(x, [v["water_ha"] for v in s], color="#2b7bba", alpha=0.35, label="Open water (ha)")
    ax2.set_ylabel("Open water (ha)")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)
    ax.set_title("Vegetation and open water through time (note the acquisition month on each date)")
    fig.tight_layout()
    return _png(fig)


def small_multiples(run: dict) -> str:
    res = run["results"]
    n = len(res)
    cols = min(5, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(2.3 * cols, 2.5 * rows), squeeze=False)
    for a in axes.ravel():
        a.axis("off")
    for a, r, s in zip(axes.ravel(), res, run["summary"]["series"]):
        with rasterio.open(r.index_paths["NDVI"]) as src:
            f = max(1, src.width // 300)
            a.imshow(src.read(1, out_shape=(src.height // f, src.width // f)), cmap="RdYlGn", vmin=-0.2, vmax=0.8)
        a.set_title(s["label"], fontsize=9)
    fig.suptitle("NDVI on every date (same colour scale)")
    fig.tight_layout()
    return _png(fig)


def change_map(run: dict) -> str | None:
    c = run["change"]
    if not c:
        return None
    with rasterio.open(c["path"]) as src:
        d = src.read(1)
    fig, (ax, axh) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
    im = ax.imshow(d, cmap="BrBG", vmin=-0.4, vmax=0.4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"NDVI change {c['from']} to {c['to']} (green = gain)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    v = d[~np.isnan(d)]
    axh.hist(v, bins=80, range=(-0.8, 0.8), color="#5b7f95")
    for t in (-c["threshold"], c["threshold"]):
        axh.axvline(t, color="#c62828", ls="--", lw=1)
    axh.set_xlabel("NDVI change"); axh.set_ylabel("pixels"); axh.set_title("Distribution, land in both dates")
    fig.tight_layout()
    return _png(fig)


def _table(rows, cols):
    head = "".join(f"<th>{html.escape(h)}</th>" for _, h in cols)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(r.get(k, '')))}</td>" for k, _ in cols) + "</tr>" for r in rows)
    return f"<table><tr>{head}</tr>{body}</table>"


def write_report(run: dict, bands: dict, out_dir: str, title: str = "Multispectral pipeline report", notes: str = "") -> str:
    s = run["summary"]
    series_cols = [("label", "Date"), ("median_ndvi_land", "Median NDVI (land)"), ("p90_ndvi_land", "NDVI 90th pct"),
                   ("dense_veg_pct_of_land", "Dense vegetation % of land"), ("water_ha", "Open water (ha)"), ("valid_ha", "Area (ha)")]
    rows = [{**r, "dense_veg_pct_of_land": f"{r['dense_veg_pct_of_land']:.1f}"} for r in s["series"]]
    c = s["change"]
    change_html = ""
    img = change_map(run)
    if c:
        change_html = (f"<h2>Change between two dates</h2><p>Mean NDVI change on land: <b>{c['mean_change']:+.3f}</b>; "
                       f"gain above {c['threshold']}: <b>{c['gain_ha']} ha</b>; loss beyond {c['threshold']}: <b>{c['loss_ha']} ha</b>.</p>"
                       f"<img src='data:image/png;base64,{img}' alt='NDVI change'>")
    latest = s["series"][-1]
    area_rows = [{"class": k.replace(" (ha)", ""), "ha": v} for k, v in latest.items() if k.endswith("(ha)") and k != "water_ha"]
    g = s["grid"]
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:1120px;margin:24px auto;padding:0 16px;color:#1d2430}}
h1{{font-size:22px;margin-bottom:4px}} h2{{font-size:17px;margin-top:28px;border-bottom:1px solid #ddd;padding-bottom:4px}}
table{{border-collapse:collapse;width:100%;font-size:14px;margin:6px 0}} td,th{{border:1px solid #ddd;padding:6px 8px;text-align:left}}
th{{background:#f3f5f8}} .meta{{color:#555;font-size:13px}} img{{max-width:100%}}
</style></head><body>
<h1>{html.escape(title)}</h1>
<p class="meta">{len(s['images'])} images · bands {html.escape(', '.join(f'{k}={v}' for k, v in s['bands'].items()))} · water rule: {s.get('water_method', '')} · aligned to {html.escape(g['reference'])} ({g['crs']}, {g['pixel_m']:g} m pixels) · zones {s['zones']['size_m']:g} m</p>
{f"<p>{html.escape(notes)}</p>" if notes else ""}
<h2>Summary by date</h2>{_table(rows, series_cols)}
<img src="data:image/png;base64,{series_chart(run)}" alt="Time series">
<img src="data:image/png;base64,{small_multiples(run)}" alt="NDVI small multiples">
{change_html}
<h2>Latest date in detail</h2>
<img src="data:image/png;base64,{latest_maps(run, bands)}" alt="Latest date maps">
<h2>Class areas, latest date</h2>{_table(area_rows, [("class", "Class"), ("ha", "Hectares")])}
<p class="meta">Outputs: rasters/ (one Cloud Optimized GeoTIFF per index and date, class rasters, change raster), timeseries.csv,
zones.geojson ({s['zones']['count']} management zones with median NDVI, dense vegetation and water share), summary.json.</p>
</body></html>"""
    path = os.path.join(out_dir, "report.html")
    with open(path, "w") as f:
        f.write(page)
    return path
