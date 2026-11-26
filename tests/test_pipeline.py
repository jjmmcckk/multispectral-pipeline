"""Tests on synthetic 4 band images with known reflectance."""
import json

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from mspipe import classify, otsu_threshold, parse_bands, process_image, reference_grid, run
from mspipe.cli import main

CRS = "EPSG:26914"


def make_image(path, nir_field=200, red_field=40, water=True, shift_m=0.0, size=200, pixel=5.0):
    """A vegetated field (NIR high, red low) with a square lake (NIR very low) in one corner."""
    r = np.full((size, size), red_field, np.uint8)
    g = np.full((size, size), 60, np.uint8)
    b = np.full((size, size), 40, np.uint8)
    n = np.full((size, size), nir_field, np.uint8)
    if water:
        r[:50, :50], g[:50, :50], b[:50, :50], n[:50, :50] = 30, 45, 50, 10
    t = from_origin(500000 + shift_m, 5200000, pixel, pixel)
    with rasterio.open(path, "w", driver="GTiff", width=size, height=size, count=4, dtype="uint8",
                       crs=CRS, transform=t) as dst:
        for i, band in enumerate([r, g, b, n], start=1):
            dst.write(band, i)
    return path


def test_parse_bands():
    assert parse_bands("R=3,G=2,B=1,RE=4,NIR=5") == {"R": 3, "G": 2, "B": 1, "RE": 4, "NIR": 5}
    with pytest.raises(ValueError):
        parse_bands("G=2,B=1")


def test_ndvi_matches_formula(tmp_path):
    img = make_image(str(tmp_path / "x_20200101.tif"), water=False)
    res = process_image(img, parse_bands("R=1,G=2,B=3,NIR=4"), reference_grid(img), str(tmp_path / "out"), block=64)
    with rasterio.open(res.index_paths["NDVI"]) as src:
        ndvi = src.read(1)
    assert np.nanmean(ndvi) == pytest.approx((200 - 40) / (200 + 40), abs=1e-4)


def test_block_size_does_not_change_results(tmp_path):
    img = make_image(str(tmp_path / "x_20200101.tif"))
    bands, grid = parse_bands("R=1,G=2,B=3,NIR=4"), reference_grid(img)
    a = process_image(img, bands, grid, str(tmp_path / "a"), block=37)
    b = process_image(img, bands, grid, str(tmp_path / "b"), block=1024)
    with rasterio.open(a.index_paths["NDVI"]) as x, rasterio.open(b.index_paths["NDVI"]) as y:
        assert np.allclose(x.read(1), y.read(1), equal_nan=True)


def test_water_area_from_scene_threshold(tmp_path):
    img = make_image(str(tmp_path / "x_20200101.tif"))
    res = process_image(img, parse_bands("R=1,G=2,B=3,NIR=4"), reference_grid(img), str(tmp_path / "out"))
    assert res.class_area_ha["Open water"] == pytest.approx(50 * 50 * 25 / 10000)   # 6.25 ha lake


def test_otsu_splits_two_groups():
    values = np.r_[np.full(500, 0.05), np.full(500, 0.8)]
    t = otsu_threshold(values)
    assert 0.05 < t < 0.8


def test_classes_follow_ndvi_breaks():
    ndvi = np.array([[0.1, 0.3, 0.5, 0.7]])
    assert classify(ndvi, None).tolist() == [[2, 3, 4, 5]]
    assert classify(ndvi, np.array([[True, False, False, False]])).tolist() == [[1, 3, 4, 5]]


def test_images_are_aligned_to_the_first_date(tmp_path):
    a = make_image(str(tmp_path / "x_20190101.tif"), nir_field=150)
    b = make_image(str(tmp_path / "x_20200101.tif"), nir_field=200, shift_m=3.0)   # 3 m offset grid
    out = run([b, a], parse_bands("R=1,G=2,B=3,NIR=4"), str(tmp_path / "out"))
    s = out["summary"]
    assert [r["date"] for r in s["series"]] == ["20190101", "20200101"]
    assert s["grid"]["reference"] == "x_20190101.tif"
    assert s["change"]["mean_change"] > 0          # field greener in 2020


def test_cli_writes_outputs(tmp_path):
    a = make_image(str(tmp_path / "x_20190101.tif"), nir_field=150)
    b = make_image(str(tmp_path / "x_20200101.tif"))
    out = tmp_path / "res"
    assert main([a, b, "--out", str(out), "--zone", "250"]) == 0
    for name in ["report.html", "timeseries.csv", "zones.geojson", "summary.json"]:
        assert (out / name).exists(), name
    zones = json.loads((out / "zones.geojson").read_text())
    assert len(zones["features"]) == 16            # 1 km square in 250 m zones


def test_reflectance_image_with_red_edge(tmp_path):
    """Drone style input: 5 bands of float reflectance (0 to 1), red edge included, NDWI water rule."""
    size, path = 100, str(tmp_path / "flight_20240715.tif")
    vals = {"B": 0.04, "G": 0.08, "R": 0.05, "RE": 0.25, "NIR": 0.45}
    arrs = {k: np.full((size, size), v, np.float32) for k, v in vals.items()}
    for k, v in {"B": 0.06, "G": 0.07, "R": 0.04, "RE": 0.03, "NIR": 0.02}.items():
        arrs[k][:20, :20] = v                                   # pond: green above NIR
    with rasterio.open(path, "w", driver="GTiff", width=size, height=size, count=5, dtype="float32",
                       crs=CRS, transform=from_origin(500000, 5200000, 1.0, 1.0)) as dst:
        for i, k in enumerate(["B", "G", "R", "RE", "NIR"], start=1):
            dst.write(arrs[k], i)
    bands = parse_bands("B=1,G=2,R=3,RE=4,NIR=5")
    res = process_image(path, bands, reference_grid(path), str(tmp_path / "out"), scale=1.0, water="ndwi")
    assert "NDRE" in res.index_paths
    with rasterio.open(res.index_paths["NDRE"]) as src:
        assert np.nanmax(src.read(1)) == pytest.approx((0.45 - 0.25) / (0.45 + 0.25), abs=1e-4)
    assert res.class_area_ha["Open water"] == pytest.approx(20 * 20 / 10000)
