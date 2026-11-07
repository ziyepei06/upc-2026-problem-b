"""Locate external data roots used by the figure / audit scripts.

The paper's own outputs live in ``../results`` and are shipped with this repo.
A few scripts additionally read two datasets that are *not* redistributed here:

* the reference library (the tracking datasets used for the noise model and the
  gold-paper PDFs)                    -> ``paths.LIB`` / ``paths.DATASETS``
* the parallel GPT-pipeline results    -> ``paths.GPT``

Point ``UPC_ROOT`` at a directory containing both, for example::

    export UPC_ROOT=$HOME/upc-2026-data
    # <UPC_ROOT>/2026/UPC 2026 Problem B 资料库/08_数据集/...
    # <UPC_ROOT>/GPT/prediction_revision_v5/results/...

Scripts that need a missing file stop with a message naming the variable, so
the paper PDF and all shipped JSON results stay reproducible without them.
"""
from __future__ import annotations

import os

_ROOT = os.environ.get("UPC_ROOT", os.path.expanduser("~/UPC-2026"))
LIB = os.path.join(_ROOT, "2026", "UPC 2026 Problem B 资料库")
DATASETS = os.path.join(LIB, "08_数据集")
GPT = os.path.join(_ROOT, "GPT", "prediction_revision_v5", "results")


def dataset(name: str) -> str:
    """Path of one file in the reference 08_数据集 directory."""
    return os.path.join(DATASETS, name)


def npz(name: str) -> str:
    """Path of one result file produced by the parallel GPT pipeline."""
    return os.path.join(GPT, name)
