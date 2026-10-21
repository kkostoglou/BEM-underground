"""Line-element geometry translated from ``CSHAPE`` and ``ELEMLEN``."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from .quadrature import leggauss


def shape_functions(
    xi: float, coordinates: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """Return shape values, derivatives, outward normal, and Jacobian.

    Node ordering follows the Fortran routines. Coordinates must contain two
    nodes for a linear element or three nodes for a quadratic element.
    """

    xy = np.asarray(coordinates, dtype=np.float64)
    if xy.shape == (2, 2):
        shape = np.array([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)])
        derivative = np.array([-0.5, 0.5])
    elif xy.shape == (3, 2):
        shape = np.array(
            [-0.5 * xi + 0.5 * xi**2, 1.0 - xi**2, 0.5 * xi + 0.5 * xi**2]
        )
        derivative = np.array([-0.5 + xi, -2.0 * xi, 0.5 + xi])
    else:
        raise ValueError("coordinates must have shape (2, 2) or (3, 2)")

    tangent = derivative @ xy
    jacobian = float(np.linalg.norm(tangent))
    if jacobian < 1.0e-11:
        raise ValueError("degenerate line element has a near-zero Jacobian")
    normal = np.array([tangent[1], -tangent[0]]) / jacobian
    return shape, derivative, normal, jacobian


def element_length(coordinates: ArrayLike, quadrature_points: int = 10) -> float:
    """Integrate line-element length using Gauss-Legendre quadrature."""

    if not 2 <= quadrature_points <= 10:
        raise ValueError("quadrature_points must be between 2 and 10")
    points, weights = leggauss(quadrature_points)
    return float(
        sum(
            shape_functions(float(point), coordinates)[3] * float(weight)
            for point, weight in zip(points, weights, strict=True)
        )
    )
