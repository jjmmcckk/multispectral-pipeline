"""Multispectral pipeline: aligned, block by block vegetation and water indices, classes, change and zones."""
from .core import CLASSES, INDICES, classify, otsu_threshold, water_mask, parse_bands, process_image, reference_grid, zone_stats
from .pipeline import ndvi_change, run
from .report import write_report

__all__ = ["CLASSES", "INDICES", "classify", "otsu_threshold", "water_mask", "parse_bands", "process_image", "reference_grid", "zone_stats",
           "ndvi_change", "run", "write_report"]
__version__ = "0.1.0"
