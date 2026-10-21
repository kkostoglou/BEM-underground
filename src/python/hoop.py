"""Boundary stress recovery from LAPLB3/TIMB2 HOOP and POSTB2 nodal averaging."""
from __future__ import annotations

import numpy as np

from .geometry import element_length, shape_functions
from .loads import _closed_boundary_inward_material_normals
from .materials import viscoelastic_properties


def incident_local_stress(mesh, global_stress, *, normals=None):
    """PREBEM ZSTRLOC: normal, tangential, and shear stresses at each node."""
    normal = _closed_boundary_inward_material_normals(mesh) if normals is None else np.asarray(normals)
    tangent = np.column_stack((-normal[:, 1], normal[:, 0]))
    return np.stack([
        np.einsum("ni,snij,nj->sn", normal, global_stress, normal),
        np.einsum("ni,snij,nj->sn", tangent, global_stress, tangent),
        np.einsum("ni,snij,nj->sn", normal, global_stress, tangent),
    ], axis=-1)


def boundary_stress(mesh, material, displacement, traction, incident_stress, parameters=None):
    """Return total local stresses (sample, element, local node, rr/tt/rt).

    Preserve the original seven-equation system, six-point element length,
    and derivative normalization by half the element length (not the local
    Jacobian). Incident stresses use the PREBEM nodal frame, as in HOOP.
    """
    displacement, traction = np.asarray(displacement), np.asarray(traction)
    samples, nodes = len(displacement), len(mesh.coordinates)
    if displacement.shape != (samples, 2 * nodes) or traction.shape != displacement.shape:
        raise ValueError("displacement and traction must have shape (samples, 2 * nodes)")
    incident_stress = np.asarray(incident_stress)
    if incident_stress.shape != (samples, nodes, 3):
        raise ValueError("incident_stress must have shape (samples, nodes, 3)")
    if parameters is not None and len(parameters) != samples:
        raise ValueError("one Laplace parameter is required per sample")
    dtype = np.complex128 if parameters is not None else np.float64
    result = np.empty((samples, len(mesh.connectivity), 3, 3), dtype=dtype)
    u, t = displacement.reshape(samples, nodes, 2), traction.reshape(samples, nodes, 2)
    for element, connectivity in enumerate(mesh.connectivity):
        xy = mesh.coordinates[connectivity]
        half_length = 0.5 * element_length(xy, 6)
        for local_node, xi in enumerate((-1., 0., 1.)):
            shape, derivative, normal, _ = shape_functions(xi, xy)
            nx, ny = normal
            qx, qy = -ny, nx
            rhs = np.zeros((samples, 7), dtype=dtype)
            rhs[:, 1:3] = np.einsum("j,sjk->sk", shape, t[:, connectivity])
            du = np.einsum("j,sjk->sk", derivative, u[:, connectivity]) / half_length
            rhs[:, 4], rhs[:, 6] = du[:, 0], du[:, 1]
            for step in range(samples):
                lam, mu = (material.lame_lambda, material.shear_modulus)
                if parameters is not None:
                    lam, mu = viscoelastic_properties(parameters[step], material)[:2]
                matrix = np.zeros((7, 7), dtype=dtype)
                matrix[0, [0, 3, 6]] = [-1, lam + 2 * mu, lam]
                matrix[1, [0, 2]] = [nx, ny]
                matrix[2, [1, 2]] = [ny, nx]
                matrix[3, [2, 4, 5]] = [-1, mu, mu]
                matrix[4, [3, 5]] = [qx, qy]
                matrix[5, [1, 3, 6]] = [-1, lam, lam + 2 * mu]
                matrix[6, [4, 6]] = [qx, qy]
                sxx, syy, sxy = np.linalg.solve(matrix, rhs[step])[:3]
                result[step, element, local_node] = (
                    sxx * nx**2 + syy * ny**2 + 2 * sxy * nx * ny,
                    sxx * ny**2 + syy * nx**2 - 2 * sxy * nx * ny,
                    (syy - sxx) * nx * ny + sxy * (nx**2 - ny**2),
                )
            result[:, element, local_node] += incident_stress[:, connectivity[local_node]]
    return result


def nodal_stress(mesh, element_stress):
    """POSTB2 default: average values from elements sharing a node."""
    values = np.asarray(element_stress)
    result = np.zeros((len(values), len(mesh.coordinates), 3), dtype=values.dtype)
    counts = np.zeros(len(mesh.coordinates))
    for element, connectivity in enumerate(mesh.connectivity):
        for local_node, node in enumerate(connectivity):
            result[:, node] += values[:, element, local_node]
            counts[node] += 1
    if np.any(counts == 0):
        raise ValueError("all nodes must belong to an element")
    return result / counts[None, :, None]
