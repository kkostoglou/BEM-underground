"""Numerical boundary-element integration translated from ``INTEGR``."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .geometry import shape_functions
from .materials import Material
from .numerics import elastostatic_kernel, laplace_domain_kernel
from .time_domain import time_domain_kernel
from .quadrature import leggauss


def integrate_elastostatic_element(
    field_point: ArrayLike,
    coordinates: ArrayLike,
    *,
    poisson: float,
    shear_modulus: float,
    quadrature_points: int = 10,
    subdivisions: int = 1,
    singular_local_node: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integrate nonsingular displacement and traction kernels over an element.

    The returned arrays have shape ``(nodes, 2, 2)``. This implements the
    nonsingular branch first; collocation points lying on the element are
    omitted for ``singular_local_node``, matching the static subtraction used
    by the legacy ``INTEGR`` routine. Even-order quadrature avoids evaluating
    the logarithmic displacement kernel at the collocation point itself.
    """

    point = np.asarray(field_point, dtype=np.float64)
    xy = np.asarray(coordinates, dtype=np.float64)
    if point.shape != (2,):
        raise ValueError("field_point must contain two coordinates")
    if xy.shape not in ((2, 2), (3, 2)):
        raise ValueError("coordinates must have shape (2, 2) or (3, 2)")
    if not 2 <= quadrature_points <= 10:
        raise ValueError("quadrature_points must be between 2 and 10")
    if subdivisions < 1:
        raise ValueError("subdivisions must be positive")
    if singular_local_node is not None and not 0 <= singular_local_node < len(xy):
        raise ValueError("singular_local_node is outside the element")
    coincident = np.flatnonzero(np.linalg.norm(xy - point, axis=1) < 1.0e-11)
    if len(coincident) and singular_local_node is None:
        raise ValueError("singular collocation requires singular_local_node")

    displacement = np.zeros((len(xy), 2, 2), dtype=np.float64)
    traction_matrix = np.zeros_like(displacement)
    gauss_points, gauss_weights = leggauss(quadrature_points)
    increment = 2.0 / subdivisions

    for subdivision in range(subdivisions):
        start = -1.0 + increment * subdivision
        end = start + increment
        midpoint = 0.5 * (start + end)
        half_width = 0.5 * increment
        for gauss_point, weight in zip(gauss_points, gauss_weights, strict=True):
            xi = midpoint + half_width * gauss_point
            shape, _, normal, jacobian = shape_functions(float(xi), xy)
            source = shape @ xy
            if np.linalg.norm(source - point) < 1.0e-11:
                # Odd-order rules can place a point directly on the
                # collocation node.  The CSHAPE/INTEGR singular branch treats
                # this isolated quadrature contribution as zero; neighbouring
                # subelements carry the finite-part integral.
                continue
            kernel_u, kernel_t = elastostatic_kernel(
                point,
                source,
                normal,
                poisson_ratio=poisson,
                shear_modulus=shear_modulus,
            )
            factor = jacobian * half_width * weight
            displacement += shape[:, None, None] * kernel_u[None, :, :] * factor
            traction_matrix += shape[:, None, None] * kernel_t[None, :, :] * factor

    if singular_local_node is not None:
        traction_matrix[singular_local_node] = 0.0

    return displacement, traction_matrix


def integrate_laplace_element(
    field_point: ArrayLike,
    coordinates: ArrayLike,
    laplace_parameter: complex,
    material: Material,
    *,
    quadrature_points: int = 10,
    subdivisions: int = 1,
    singular_local_node: int | None = None,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Integrate ``KERNEL2``, including the legacy singular subtraction."""

    point = np.asarray(field_point, dtype=np.float64)
    xy = np.asarray(coordinates, dtype=np.float64)
    if point.shape != (2,) or xy.shape not in ((2, 2), (3, 2)):
        raise ValueError("invalid field point or element coordinates")
    if not 2 <= quadrature_points <= 10 or subdivisions < 1:
        raise ValueError("invalid quadrature_points or subdivisions")
    if singular_local_node is not None and not 0 <= singular_local_node < len(xy):
        raise ValueError("singular_local_node is outside the element")

    displacement = np.zeros((len(xy), 2, 2), dtype=np.complex128)
    traction_matrix = np.zeros_like(displacement)
    gauss_points, gauss_weights = leggauss(quadrature_points)
    increment = 2.0 / subdivisions
    poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))

    for subdivision in range(subdivisions):
        midpoint = -1.0 + increment * (subdivision + 0.5)
        half_width = 0.5 * increment
        for gauss_point, weight in zip(gauss_points, gauss_weights, strict=True):
            xi = midpoint + half_width * gauss_point
            shape, _, normal, jacobian = shape_functions(float(xi), xy)
            integration_point = shape @ xy
            dynamic_u, dynamic_t = laplace_domain_kernel(
                point, integration_point, normal, laplace_parameter, material
            )
            factor = jacobian * half_width * weight
            displacement += shape[:, None, None] * dynamic_u[None, :, :] * factor
            for local_node in range(len(xy)):
                traction_value = dynamic_t
                if local_node == singular_local_node:
                    _, static_t = elastostatic_kernel(
                        point,
                        integration_point,
                        normal,
                        poisson_ratio=poisson,
                        shear_modulus=material.shear_modulus,
                    )
                    traction_value = dynamic_t - static_t
                traction_matrix[local_node] += shape[local_node] * traction_value * factor
    return displacement, traction_matrix


def integrate_time_element(
    field_point: ArrayLike, coordinates: ArrayLike, step: int,
    total_time: float, step_count: int, material: Material, *,
    quadrature_points: int, subdivisions: int,
    singular_local_node: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integrate the direct time-domain kernel for one impulse step."""

    point = np.asarray(field_point, dtype=float); xy = np.asarray(coordinates, dtype=float)
    displacement = np.zeros((len(xy), 2, 2)); traction_matrix = np.zeros_like(displacement)
    points, weights = leggauss(quadrature_points)
    increment = 2.0 / subdivisions
    poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))
    for subdivision in range(subdivisions):
        midpoint = -1.0 + increment * (subdivision + 0.5); half_width = increment / 2.0
        for gauss_point, weight in zip(points, weights, strict=True):
            xi = midpoint + half_width * gauss_point
            shape, _, normal, jacobian = shape_functions(float(xi), xy)
            integration_point = shape @ xy
            kernel_u, kernel_t = time_domain_kernel(
                point, integration_point, normal, step, total_time, step_count, material
            )
            factor = jacobian * half_width * weight
            displacement += shape[:, None, None] * kernel_u[None] * factor
            for local_node in range(len(xy)):
                traction_value = kernel_t
                if step == 1 and local_node == singular_local_node:
                    _, static_t = elastostatic_kernel(
                        point, integration_point, normal,
                        poisson_ratio=poisson, shear_modulus=material.shear_modulus,
                    )
                    traction_value = kernel_t - static_t
                traction_matrix[local_node] += shape[local_node] * traction_value * factor
    return displacement, traction_matrix
