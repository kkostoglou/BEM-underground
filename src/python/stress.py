"""Stress transformations translated from ``PRETIM4.FOR``."""

from __future__ import annotations

from typing import TypeVar


Scalar = TypeVar("Scalar", float, complex)


def stress_global(
    sxx_local: Scalar,
    syy_local: Scalar,
    sxy_local: Scalar,
    cosine: float,
    sine: float,
) -> tuple[Scalar, Scalar, Scalar]:
    sxx = sxx_local * cosine**2 - 2.0 * sxy_local * cosine * sine + syy_local * sine**2
    syy = sxx_local * sine**2 + 2.0 * sxy_local * cosine * sine + syy_local * cosine**2
    sxy = (sxx_local - syy_local) * cosine * sine + sxy_local * (cosine**2 - sine**2)
    return sxx, syy, sxy


def stress_local(
    sxx: Scalar, syy: Scalar, sxy: Scalar, cosine: float, sine: float
) -> tuple[Scalar, Scalar, Scalar]:
    sxx_local = sxx * cosine**2 + 2.0 * sxy * cosine * sine + syy * sine**2
    syy_local = sxx * sine**2 - 2.0 * sxy * cosine * sine + syy * cosine**2
    sxy_local = (syy - sxx) * cosine * sine + sxy * (cosine**2 - sine**2)
    return sxx_local, syy_local, sxy_local


def traction(
    sxx: Scalar, syy: Scalar, sxy: Scalar, normal_x: float, normal_y: float
) -> tuple[Scalar, Scalar]:
    return sxx * normal_x + sxy * normal_y, sxy * normal_x + syy * normal_y
