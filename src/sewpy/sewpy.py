#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sewpy: Thin Python wrapper for running SExtractor on FITS images.

Highlights
----------
- **Robust parameter & PSF data loading**:
  * Parameter names are read from the packaged `sextractor_params.txt` using
    `importlib.resources`. If not found, the module tries `sex -dp`, and finally
    a local fallback file.
  * Default PSF data come from the packaged `psf_default.npy`. If unavailable,
    it falls back to a neighbor `psf_array.py:arr` or local `.npy` file.
- **Clean design**: type hints, explicit behavior, minimal surprises.
- **Wheel- & editable-safe**: data files resolved via importlib rather than
  filesystem-relative paths.

Requirements
------------
Python ≥ 3.8, NumPy ≥ 1.23, Astropy ≥ 5.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union
import copy
import logging
import os
import re
import subprocess
import tempfile
import importlib.util
from importlib.resources import files as pkg_files

import numpy as np
import astropy.table as apt
from astropy.io import fits

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------

def _try_import_arr_from_neighbor(py_path: Path, module_name: str = "psf_array"):
    """Try to import a neighbor psf_array.py and return its `arr` if defined."""
    cand = py_path.with_name(f"{module_name}.py")
    if cand.is_file():
        spec = importlib.util.spec_from_file_location(module_name, cand)
        if spec and spec.loader:
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)  # type: ignore[attr-defined]
            if hasattr(mod, "arr"):
                return getattr(mod, "arr")
    return None


def _default_param_registry(sexpath: str) -> ParamRegistry:
    """
    Prefer the packaged 'sextractor_params.txt'; fallback to neighbor file or
    auto-discovery via `sex -dp`.
    """
    # 1) Packaged data
    try:
        p = pkg_files("sewpy").joinpath("sextractor_params.txt")
        return ParamRegistry(filepath=str(p))
    except Exception:
        pass

    # 2) Local development fallback
    here = Path(__file__).resolve()
    txt = here.with_name("sextractor_params.txt")
    if txt.is_file():
        return ParamRegistry(filepath=str(txt))

    # 3) Last resort
    return ParamRegistry(discover_from_sexpath=True)


# -------------------------------------------------------------------------
# Defaults
# -------------------------------------------------------------------------

DEFAULT_PARAMS: List[str] = [
    "XWIN_IMAGE",
    "YWIN_IMAGE",
    "AWIN_IMAGE",
    "BWIN_IMAGE",
    "THETAWIN_IMAGE",
    "BACKGROUND",
    "FLUX_AUTO",
]

DEFAULT_CONFIG: Dict[str, Union[str, float, int]] = {}


# -------------------------------------------------------------------------
# Param registry
# -------------------------------------------------------------------------

@dataclass
class ParamRegistry:
    """Holds a list of known SExtractor parameters for sanity checks."""

    names: Optional[List[str]] = None
    filepath: Optional[Union[str, Path]] = None
    discover_from_sexpath: bool = False

    def load(self, sexpath: str) -> List[str]:
        if self.names:
            return list(self.names)

        # Try from file
        if self.filepath:
            p = Path(self.filepath)
            if p.is_file():
                text = p.read_text(encoding="utf-8", errors="ignore")
                found = re.findall(r"#(\w+)", text)
                if not found:
                    found = [
                        ln.strip() for ln in text.splitlines()
                        if ln.strip() and not ln.strip().startswith("#")
                    ]
                return found

        # Try discovery via `sex -dp`
        if self.discover_from_sexpath:
            try:
                p = subprocess.Popen(
                    [sexpath, "-dp"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                out, err = p.communicate()
                text = (out + err).decode("utf-8", errors="ignore")
                found = re.findall(r"#(\w+)", text)
                if found:
                    return found
            except Exception:
                pass

        return []


# -------------------------------------------------------------------------
# Core class
# -------------------------------------------------------------------------

class SEW:
    """SExtractor runner with light conveniences."""

    def __init__(
        self,
        workdir: Optional[Union[str, Path]] = None,
        sexpath: str = "sex",
        params: Optional[List[str]] = None,
        config: Optional[Dict[str, Union[str, float, int]]] = None,
        configfilepath: Optional[Union[str, Path]] = None,
        nice: Optional[int] = None,
        loglevel: Optional[Union[int, str]] = None,
        param_registry: Optional[ParamRegistry] = None,
        psf_array: Optional[np.ndarray] = None,
        psf_filepath: Optional[Union[str, Path]] = None,
    ) -> None:
        if loglevel is not None:
            logger.setLevel(loglevel)

        self.sexpath = sexpath
        self.configfilepath = Path(configfilepath) if configfilepath else None
        self.nice = nice

        logger.info("SExtractor version: %s", self.get_version())

        # Working directory
        if workdir is not None:
            self.workdir = Path(workdir)
            self.workdir.mkdir(parents=True, exist_ok=True)
            self.tmp = False
        else:
            self.workdir = Path(tempfile.mkdtemp(prefix="sewpy_"))
            self.tmp = True

        # Parameters and config
        self.params = list(params) if params is not None else list(DEFAULT_PARAMS)
        self.config = dict(config) if config is not None else dict(DEFAULT_CONFIG)

        # Registry
        self._registry = param_registry or _default_param_registry(self.sexpath)
        self.fullparamlist = self._registry.load(self.sexpath)

        # PSF resources
        self._psf_array = psf_array
        self._psf_filepath = Path(psf_filepath) if psf_filepath else None

        # Load default PSF if nothing provided
        if "PSF_NAME" not in self.config and self._psf_array is None and self._psf_filepath is None:
            maybe_arr = None
            # 1) packaged npy
            try:
                npy = pkg_files("sewpy").joinpath("psf_default.npy")
                with npy.open("rb") as fh:
                    maybe_arr = np.load(fh)
            except Exception:
                pass
            # 2) neighbor psf_array.py
            if maybe_arr is None:
                maybe_arr = _try_import_arr_from_neighbor(Path(__file__).resolve())
            # 3) neighbor npy
            if maybe_arr is None:
                npy_fs = Path(__file__).with_name("psf_default.npy")
                if npy_fs.is_file():
                    try:
                        maybe_arr = np.load(str(npy_fs))
                    except Exception:
                        pass
            if maybe_arr is not None:
                self._psf_array = np.asarray(maybe_arr)
                self.config["PSF_NAME"] = str(self._psf_path())

        self._set_instance_config()
        self._check_config_stub()
        self._check_params_against_registry()

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    def get_version(self) -> str:
        """Return the SExtractor version string."""
        try:
            p = subprocess.Popen(
                [self.sexpath],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            out, err = p.communicate()
        except Exception as exc:
            raise RuntimeError(f"Cannot run SExtractor at '{self.sexpath}'") from exc
        text = (out or b"") + (err or b"")
        m = re.search(rb"[Vv]ersion (\d+(?:\.\d+)*)", text)
        if not m:
            raise RuntimeError(f"Cannot determine version from '{self.sexpath}'.")
        return m.group(1).decode("utf-8")

    def __str__(self):
        return f"SEW(workdir={self.workdir})"

    def _params_path(self) -> Path: return self.workdir / "params.txt"
    def _config_path(self) -> Path: return self.configfilepath or (self.workdir / "config.txt")
    def _conv_path(self) -> Path: return self.workdir / "conv.txt"
    def _psf_path(self) -> Path: return self.workdir / "default.psf"
    def _cat_path(self, imgname: str) -> Path: return self.workdir / f"{imgname}.cat.txt"
    def _assoc_path(self, imgname: str) -> Path: return self.workdir / f"{imgname}.assoc.txt"
    def _log_path(self, imgname: str) -> Path: return self.workdir / f"{imgname}.log.txt"

    # -----------------------------------------------------------------
    # Setup
    # -----------------------------------------------------------------

    def _set_instance_config(self) -> None:
        """Seed instance-level config keys."""
        if "PARAMETERS_NAME" not in self.config:
            self.config["PARAMETERS_NAME"] = str(self._params_path())

        if "FILTER_NAME" not in self.config:
            self.config["FILTER_NAME"] = str(self._conv_path())

        if "PSF_NAME" not in self.config:
            if self._psf_filepath and self._psf_filepath.exists():
                self.config["PSF_NAME"] = str(self._psf_filepath)
            elif self._psf_array is not None:
                self.config["PSF_NAME"] = str(self._psf_path())

        self.config.pop("CATALOG_NAME", None)  # managed per-image

    def _check_config_stub(self): return

    def _check_params_against_registry(self) -> None:
        if not self.fullparamlist:
            return
        bad = []
        for p in self.params:
            m = re.match(r"(\w+)\(\d+\)", p)
            clean = m.group(1) if m else p
            if clean not in self.fullparamlist:
                bad.append(p)
        if bad:
            logger.warning("Unknown SExtractor params: %s", ", ".join(bad))

    # -----------------------------------------------------------------
    # Writers
    # -----------------------------------------------------------------

    def _write_params(self) -> None:
        self._params_path().write_text("\n".join(self.params) + "\n", encoding="utf-8")

    def _write_default_config(self) -> None:
        if self.configfilepath:
            return
        path = self._config_path()
        if not path.exists():
            p = subprocess.Popen([self.sexpath, "-dd"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            path.write_bytes(out or b"")

    def _write_default_conv(self) -> None:
        path = self._conv_path()
        if not path.exists():
            path.write_text("CONV NORM\n1 2 1\n2 4 2\n1 2 1\n", encoding="utf-8")

    def _write_psf_if_needed(self) -> None:
        if self._psf_array is None:
            return
        path = Path(self.config.get("PSF_NAME", self._psf_path()))
        arr = np.asarray(self._psf_array, dtype=np.float32)
        if arr.ndim == 2:
            arr = arr[..., None]
        if arr.ndim == 3 and arr.shape[0] != arr.shape[1] and arr.shape[-1] == arr.shape[-2]:
            arr = np.moveaxis(arr, 0, -1)
        if arr.ndim != 3:
            raise ValueError("psf_array must be 2D or 3D.")
        ny, nx, nz = arr.shape
        flat = arr.reshape(ny * nx * nz)
        col = fits.Column(
            name="PSF_MASK",
            format=f"{ny*nx*nz}E",
            dim=f"({nx},{ny},{nz})",
            array=[flat]
        )
        tbhdu = fits.BinTableHDU.from_columns([col])
        fits.HDUList([fits.PrimaryHDU(), tbhdu]).writeto(path, overwrite=True)

    # -----------------------------------------------------------------
    # Core run
    # -----------------------------------------------------------------

    def __call__(
        self,
        imgfilepath: Union[str, Path],
        imgname: Optional[str] = None,
        assoc_cat: Optional[apt.Table] = None,
        assoc_xname: str = "x",
        assoc_yname: str = "y",
        returncat: bool = True,
        prefix: str = "",
        writelog: bool = True,
    ) -> Dict[str, Union[str, Path, apt.Table]]:
        start = datetime.now()
        imgpath = Path(imgfilepath)
        if not imgpath.exists():
            raise FileNotFoundError(f"Image not found: {imgpath}")
        if imgname is None:
            imgname = imgpath.stem

        imgconfig = copy.deepcopy(self.config)
        imgconfig["CATALOG_NAME"] = str(self._cat_path(imgname))

        self._write_default_config()
        self._write_params()
        self._write_default_conv()
        self._write_psf_if_needed()

        cmd = [self.sexpath, str(imgpath), "-c", str(self._config_path())]
        for k, v in imgconfig.items():
            cmd += [f"-{k}", str(v).replace(" ", "")]
        if self.nice is not None:
            cmd = ["nice", "-n", str(self.nice)] + cmd

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        if writelog:
            self._log_path(imgname).write_text(
                (out or b"").decode("utf-8", "ignore") +
                "\n####### stderr #######\n" +
                (err or b"").decode("utf-8", "ignore")
            )

        cat_path = self._cat_path(imgname)
        if not cat_path.exists():
            raise RuntimeError(f"SExtractor failed to produce catalog {cat_path}")

        result: Dict[str, Union[str, Path, apt.Table]] = {
            "catfilepath": str(cat_path),
            "workdir": str(self.workdir),
        }
        if writelog:
            result["logfilepath"] = str(self._log_path(imgname))
        if returncat:
            result["table"] = apt.Table.read(str(cat_path), format="ascii.sextractor")

        elapsed = (datetime.now() - start).total_seconds()
        logger.info("SExtractor completed in %.2f s", elapsed)
        return result