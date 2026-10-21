"""Evaluate zone kernels in supplied Fortran KZON element orientation."""
import numpy as np
from .assembly import _legacy_quadrature
from .integration import integrate_laplace_element, integrate_elastostatic_element


def evaluate_zone_laplace_kernels(coordinates, connectivity, parameter, materials):
    """Return unscaled traction, displacement, and static correction arrays.

    KZON uses one-based node numbers, zero triples for absent elements.
    Material order must match KZON zone order. Connectivity orientation must
    already be established by preprocessing. Uses the existing Python
    quadrature policy; custom BEMZ NGAUS/NSNG controls are not implemented here.
    """
    xy = np.asarray(coordinates, dtype=float)
    conn = np.asarray(connectivity)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
        raise ValueError("coordinates must be finite node pairs")
    if conn.ndim != 3 or conn.shape[2] != 3 or len(conn) != len(materials):
        raise ValueError("quadratic connectivity and material zone counts must agree")
    if not np.issubdtype(conn.dtype, np.integer) or np.any(conn < 0) or np.any(conn > len(xy)):
        raise ValueError("invalid Fortran node numbers")
    active = np.any(conn != 0, axis=2)
    if np.any(active & np.any(conn == 0, axis=2)):
        raise ValueError("active elements must have three node numbers")
    shape = (len(conn), len(xy), conn.shape[1], 3, 2, 2)
    h = np.zeros(shape, dtype=complex)
    g = np.zeros_like(h)
    diagonal = np.zeros((len(conn), len(xy), 2, 2))
    for z, material in enumerate(materials):
        poisson = material.lame_lambda / (2 * (material.lame_lambda + material.shear_modulus))
        for node in np.unique(conn[z]):
            if node == 0:
                continue
            point = xy[node - 1]
            for element in np.flatnonzero(active[z]):
                vertices = xy[conn[z, element] - 1]
                points, subdivisions, singular = _legacy_quadrature(point, vertices)
                options = dict(quadrature_points=points, subdivisions=subdivisions,
                               singular_local_node=singular)
                displacement, traction = integrate_laplace_element(
                    point, vertices, parameter, material, **options)
                _, static = integrate_elastostatic_element(point, vertices,
                    poisson=poisson, shear_modulus=material.shear_modulus, **options)
                g[z, node - 1, element] = displacement
                h[z, node - 1, element] = traction
                diagonal[z, node - 1] -= static.sum(axis=0)
    return h, g, diagonal
