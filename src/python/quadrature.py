"""Small, dependency-light Gauss-Legendre quadrature rules.

The NumPy ``leggauss`` implementation calls a native eigensolver on recent
Windows builds.  Newton iteration is sufficient for the low orders used by
the translated Fortran routines and avoids that fragile native dependency.
"""
from __future__ import annotations

import math
import numpy as np
from numpy.typing import NDArray


def leggauss(order: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if order < 1:
        raise ValueError("order must be positive")
    points = np.zeros(order, dtype=float)
    weights = np.zeros(order, dtype=float)
    for i in range((order + 1) // 2):
        x = math.cos(math.pi * (i + 0.75) / (order + 0.5))
        for _ in range(50):
            p0, p1 = 1.0, x
            for degree in range(2, order + 1):
                p0, p1 = p1, ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree
            derivative = order * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 1.0e-15:
                break
        points[i] = -x
        points[-i - 1] = x
        weight = 2.0 / ((1.0 - x * x) * derivative * derivative)
        weights[i] = weight
        weights[-i - 1] = weight
    return points, weights
