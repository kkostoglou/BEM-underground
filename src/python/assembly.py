"""Collocation matrix assembly for the elastostatic boundary equations."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .integration import integrate_elastostatic_element, integrate_laplace_element, integrate_time_element
from .materials import Material
from .mesh import BoundaryMesh
from .geometry import element_length


def _legacy_quadrature(point: NDArray[np.float64], coordinates: NDArray[np.float64]) -> tuple[int, int, int | None]:
    """Return NGP, NSUB, and IN using the selection in NMAIN/BMAIN."""

    distances = np.linalg.norm(coordinates - point, axis=1)
    local = np.flatnonzero(distances < 1.0e-2)
    if len(local):
        return 6, 4, int(local[-1])
    minimum = float(np.min(distances))
    length = element_length(coordinates, quadrature_points=6)
    ratio = minimum / length
    subdivisions = 3 if ratio < 1.0 else (2 if ratio <= 2.0 else 1)
    sub_length = length / subdivisions
    ratio_one = 1.0e-4 * sub_length
    ratio_two = sub_length / (4.0 * minimum)
    estimate = 0.5 * np.log(ratio_one) / np.log(ratio_two) + 1.0
    points = max(2, min(10, int(estimate)))
    return points, subdivisions, None


def assemble_elastostatic(
    mesh: BoundaryMesh, *, poisson: float, shear_modulus: float, progress=None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Assemble traction ``H`` and displacement ``G`` matrices.

    The free-term diagonal of ``H`` is completed with the same rigid-translation
    identity used by the Fortran code: each 2-by-2 block row sums to zero.
    """

    node_count = len(mesh.coordinates)
    size = 2 * node_count
    h_matrix = np.zeros((size, size), dtype=np.float64)
    g_matrix = np.zeros_like(h_matrix)

    for collocation in range(node_count):
        row = slice(2 * collocation, 2 * collocation + 2)
        point = mesh.coordinates[collocation]
        for element_index, connectivity in enumerate(mesh.connectivity):
            if progress is not None:
                progress(collocation + 1, node_count, element_index + 1, len(mesh.connectivity))
            coordinates = mesh.coordinates[connectivity]
            quadrature_points, subdivisions, singular_node = _legacy_quadrature(point, coordinates)
            displacement, traction = integrate_elastostatic_element(
                point,
                coordinates,
                poisson=poisson,
                shear_modulus=shear_modulus,
                quadrature_points=quadrature_points,
                subdivisions=subdivisions,
                singular_local_node=singular_node,
            )
            for local_node, global_node in enumerate(connectivity):
                column = slice(2 * global_node, 2 * global_node + 2)
                g_matrix[row, column] += displacement[local_node]
                h_matrix[row, column] += traction[local_node]

        diagonal = slice(2 * collocation, 2 * collocation + 2)
        block_sum = sum(
            (h_matrix[row, slice(2 * node, 2 * node + 2)] for node in range(node_count)),
            start=np.zeros((2, 2)),
        )
        h_matrix[row, diagonal] -= block_sum

    return h_matrix, g_matrix


def assemble_time_domain_step(
    mesh: BoundaryMesh, step: int, total_time: float, step_count: int,
    material: Material, *, jump: float = 1.0,
    node_types: NDArray[np.int64] | None = None,
    progress=None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Assemble BEMTIM impulse matrices for one 1-based lag step."""

    node_count = len(mesh.coordinates); size = 2 * node_count
    types = np.zeros(node_count, dtype=np.int64) if node_types is None else np.asarray(node_types)
    if types.shape != (node_count,):
        raise ValueError("node_types must have one entry per mesh node")
    h_matrix = np.zeros((size, size)); g_matrix = np.zeros_like(h_matrix)
    for collocation in range(node_count):
        row = slice(2 * collocation, 2 * collocation + 2); point = mesh.coordinates[collocation]
        for element_index, connectivity in enumerate(mesh.connectivity):
            if progress is not None:
                progress(collocation + 1, node_count, element_index + 1,
                         len(mesh.connectivity), "kernels")
            coordinates = mesh.coordinates[connectivity]
            points, subdivisions, singular = _legacy_quadrature(point, coordinates)
            displacement, traction = integrate_time_element(
                point, coordinates, step, total_time, step_count, material,
                quadrature_points=points, subdivisions=subdivisions,
                singular_local_node=singular,
            )
            for local_node, global_node in enumerate(connectivity):
                column = slice(2 * global_node, 2 * global_node + 2)
                g_matrix[row, column] += material.shear_modulus * displacement[local_node]
                h_matrix[row, column] += traction[local_node]
    if step == 1:
        poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))
        static_h, _ = assemble_elastostatic(
            mesh, poisson=poisson, shear_modulus=material.shear_modulus,
            progress=(lambda node, nodes, element, elements:
                      progress(node, nodes, element, elements, "static correction"))
            if progress is not None else None,
        )
        for node in range(node_count):
            block = slice(2 * node, 2 * node + 2)
            local_jump = 0.5 if types[node] == 2 else jump
            h_matrix[block, block] += static_h[block, block] + local_jump * np.eye(2)
    return h_matrix, g_matrix


def assemble_laplace_domain(
    mesh: BoundaryMesh,
    laplace_parameter: complex,
    material: Material,
    *,
    jump: float = 1.0,
    node_types: NDArray[np.int64] | None = None,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Assemble transformed ``H`` and shear-scaled ``G`` matrices.

    ``node_types`` uses PREBEM ``FREESURF`` codes.  Type 2 free-surface
    collocation has the legacy half jump term; ordinary nodes use ``jump``.
    Enclosing/interface rows require their separate legacy system and are not
    accepted by this homogeneous assembler.
    """

    node_count = len(mesh.coordinates)
    size = 2 * node_count
    h_matrix = np.zeros((size, size), dtype=np.complex128)
    g_matrix = np.zeros_like(h_matrix)
    if node_types is None:
        types = np.zeros(node_count, dtype=np.int64)
    else:
        types = np.asarray(node_types, dtype=np.int64)
        if types.shape != (node_count,):
            raise ValueError("node_types must have one entry per mesh node")
        if np.any(np.isin(types, (3, 4))):
            raise NotImplementedError("simulation and enclosing boundary matrix equations are not translated")
    poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))
    static_h, _ = assemble_elastostatic(
        mesh, poisson=poisson, shear_modulus=material.shear_modulus
    )

    for collocation in range(node_count):
        row = slice(2 * collocation, 2 * collocation + 2)
        point = mesh.coordinates[collocation]
        for connectivity in mesh.connectivity:
            coordinates = mesh.coordinates[connectivity]
            quadrature_points, subdivisions, singular_node = _legacy_quadrature(point, coordinates)
            displacement, traction = integrate_laplace_element(
                point,
                coordinates,
                laplace_parameter,
                material,
                quadrature_points=quadrature_points,
                subdivisions=subdivisions,
                singular_local_node=singular_node,
            )
            for local_node, global_node in enumerate(connectivity):
                column = slice(2 * global_node, 2 * global_node + 2)
                g_matrix[row, column] += material.shear_modulus * displacement[local_node]
                h_matrix[row, column] += traction[local_node]

        diagonal = slice(2 * collocation, 2 * collocation + 2)
        local_jump = 0.5 if types[collocation] == 2 else jump
        h_matrix[row, diagonal] += static_h[row, diagonal] + local_jump * np.eye(2)
    return h_matrix, g_matrix


def assemble_laplace_internal(
    mesh: BoundaryMesh, points: NDArray[np.float64], laplace_parameter: complex,
    material: Material,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Form the ``GENRINT`` boundary operators for interior points.

    The returned matrices are respectively the traction (``F``) and the
    shear-scaled displacement (``G``) operators.  Unlike boundary
    collocation, an interior point has no jump term.
    """

    targets = np.asarray(points, dtype=float)
    if targets.ndim != 2 or targets.shape[1] != 2:
        raise ValueError("points must have shape (count, 2)")
    size = 2 * len(mesh.coordinates)
    h_matrix = np.zeros((2 * len(targets), size), dtype=np.complex128)
    g_matrix = np.zeros_like(h_matrix)
    for index, point in enumerate(targets):
        row = slice(2 * index, 2 * index + 2)
        for connectivity in mesh.connectivity:
            coordinates = mesh.coordinates[connectivity]
            quadrature_points, subdivisions, singular_node = _legacy_quadrature(point, coordinates)
            if singular_node is not None:
                raise ValueError("an internal point must not lie on the boundary")
            displacement, traction = integrate_laplace_element(
                point, coordinates, laplace_parameter, material,
                quadrature_points=quadrature_points, subdivisions=subdivisions,
            )
            for local_node, global_node in enumerate(connectivity):
                column = slice(2 * global_node, 2 * global_node + 2)
                g_matrix[row, column] += material.shear_modulus * displacement[local_node]
                h_matrix[row, column] += traction[local_node]
    return h_matrix, g_matrix


def assemble_time_internal(
    mesh: BoundaryMesh, points: NDArray[np.float64], step: int,
    total_time: float, step_count: int, material: Material,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Assemble interior-point BEM impulse operators (no free-term jump)."""
    targets = np.asarray(points, dtype=float)
    if targets.ndim != 2 or targets.shape[1] != 2:
        raise ValueError("points must have shape (count, 2)")
    size = 2 * len(mesh.coordinates)
    h_matrix = np.zeros((2 * len(targets), size))
    g_matrix = np.zeros_like(h_matrix)
    for index, point in enumerate(targets):
        row = slice(2 * index, 2 * index + 2)
        for connectivity in mesh.connectivity:
            coordinates = mesh.coordinates[connectivity]
            quadrature_points, subdivisions, singular_node = _legacy_quadrature(point, coordinates)
            if singular_node is not None:
                raise ValueError("an internal point must not lie on the boundary")
            displacement, traction = integrate_time_element(
                point, coordinates, step, total_time, step_count, material,
                quadrature_points=quadrature_points, subdivisions=subdivisions,
            )
            for local_node, global_node in enumerate(connectivity):
                column = slice(2 * global_node, 2 * global_node + 2)
                g_matrix[row, column] += material.shear_modulus * displacement[local_node]
                h_matrix[row, column] += traction[local_node]
    return h_matrix, g_matrix
