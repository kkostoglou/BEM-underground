"""Small readers for legacy list-directed Fortran validation artifacts."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
from numpy.typing import NDArray


_COMPLEX = re.compile(
    r"\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)\s*,\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)\s*\)"
)
_REAL = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?")


def read_complex_records(path: str | Path, record_width: int) -> NDArray[np.complex128]:
    """Read list-directed ``(real,imag)`` records into a two-dimensional array."""

    text = Path(path).read_text(encoding="ascii")
    values = [
        complex(float(real.replace("D", "E").replace("d", "e")),
                float(imag.replace("D", "E").replace("d", "e")))
        for real, imag in _COMPLEX.findall(text)
    ]
    if record_width <= 0 or len(values) % record_width:
        raise ValueError("complex record count is not divisible by record_width")
    return np.asarray(values, dtype=np.complex128).reshape(-1, record_width)


def read_pro_displacements(path: str | Path, vector_width: int) -> NDArray[np.complex128]:
    """Extract the printed complex displacement vectors from a BEMLAPL ``.PRO`` file."""

    text = Path(path).read_text(encoding="ascii")
    starts = list(re.finditer(r"DISPLACEMENTS AT STEP=\s*\d+", text))
    records: list[np.ndarray] = []
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.end():end]
        next_step = re.search(r"TIME STEP\s*:", block)
        if next_step:
            block = block[: next_step.start()]
        numbers = [float(value.replace("D", "E").replace("d", "e")) for value in _REAL.findall(block)]
        required = 2 * vector_width
        if len(numbers) < required:
            raise ValueError("incomplete displacement record in PRO file")
        pairs = np.asarray(numbers[:required]).reshape(vector_width, 2)
        records.append(pairs[:, 0] + 1j * pairs[:, 1])
    if not records:
        raise ValueError("no displacement records found in PRO file")
    return np.asarray(records, dtype=np.complex128)


def read_real_records(path: str | Path) -> list[NDArray[np.float64]]:
    """Read one list-directed real-array record per physical text line."""

    lines = Path(path).read_text(encoding="ascii").splitlines()
    records = [np.fromstring(line.replace("D", "E").replace("d", "e"), sep=" ") for line in lines]
    return [record for record in records if len(record)]


def read_pro_real_displacements(path: str | Path, vector_width: int) -> NDArray[np.float64]:
    """Extract real displacement vectors from a BEMTIM ``.PRO`` report."""

    text = Path(path).read_text(encoding="ascii")
    starts = list(re.finditer(r"DISPLACEMENTS AT STEP=\s*\d+", text))
    records = []
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.end():end]
        next_step = re.search(r"TIME STEP\s*:", block)
        if next_step: block = block[:next_step.start()]
        values = [float(value.replace("D", "E").replace("d", "e")) for value in _REAL.findall(block)]
        if len(values) < vector_width:
            raise ValueError("incomplete real displacement record in PRO file")
        records.append(values[:vector_width])
    return np.asarray(records, dtype=np.float64)
