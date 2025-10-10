#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Minimal demo script for sewpy — Source Extractor Wrapper for Python.
This script can be run directly from the examples/ directory.
"""
from sewpy import SEW
import sys
import os
import logging

# Allow this demo to find sewpy even if not installed
sys.path.insert(0, os.path.abspath("../"))

# Configure basic logging (prints debug info to stdout)
logging.basicConfig(
    format="%(levelname)s: %(name)s(%(funcName)s): %(message)s",
    level=logging.DEBUG,
)

# Initialize SEW instance with example parameters
sew = SEW(
    params=["X_IMAGE", "Y_IMAGE", "FLUX_APER(3)", "FLAGS"],
    config={"DETECT_MINAREA": 10, "PHOT_APERTURES": "5, 10, 20"},
    sexpath="sex",  # assumes SExtractor is available as 'sex' in PATH
)

# If SExtractor is not available in PATH, specify full path:
# sexpath="/usr/local/bin/sextractor"

# Run SExtractor on a FITS image (replace with your actual file)
out = sew("image.fits")

# Print resulting catalog (astropy.table.Table)
print(out["table"])
