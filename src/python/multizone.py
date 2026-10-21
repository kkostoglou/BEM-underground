"""MULTI-ZONE block-system scaffolding translated from BEMZ1-BEMZ4.

This module defines the global ordering used by the zone matrix assembly.  The
kernel and interface coefficient calculations are added independently so the
validated single-zone solver is not changed while the translation proceeds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .workflow import MultiZoneModel
from .mesh import BoundaryMesh, combine_meshes


@dataclass(frozen=True)
class MultiZoneBlocks:
    """Global block offsets for displacement/traction unknowns."""

    node_count: int
    boundary_count: int
    zone_count: int
    offsets: tuple[int, ...]
    size: int


def build_block_layout(model: MultiZoneModel) -> MultiZoneBlocks:
    """Build the deterministic global layout used by multi-zone assembly."""
    node_count = len(model.mesh.coordinates)
    boundary_count = len(model.boundary_zones)
    zone_count = len(model.materials)
    # Two displacement/traction components per boundary node, with one block
    # per zone interface.  This is the layout consumed by the upcoming GENR
    # and BOUND translations.
    block_size = 2 * node_count
    offsets = tuple(index * block_size for index in range(zone_count + 1))
    return MultiZoneBlocks(node_count, boundary_count, zone_count, offsets,
                           (zone_count + 1) * block_size)


def allocate_block_system(model: MultiZoneModel) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Allocate the complex global matrices for a MULTI-ZONE Laplace solve."""
    layout = build_block_layout(model)
    return (np.zeros((layout.size, layout.size), dtype=np.complex128),
            np.zeros((layout.size, layout.size), dtype=np.complex128))


def material_for_boundary(model: MultiZoneModel, boundary_number: int):
    """Return the material assigned to a one-based boundary number."""
    if boundary_number < 1 or boundary_number > len(model.boundary_zones):
        raise IndexError(f"boundary number {boundary_number} is outside the model")
    zone_id = int(model.boundary_zones[boundary_number - 1])
    return model.materials.get(zone_id)


def assemble_zone_local_operators(
    model: MultiZoneModel, laplace_parameter: complex, *, jump: float = 1.0
) -> dict[int, tuple[NDArray[np.complex128], NDArray[np.complex128]]]:
    """Assemble validated per-zone operators for the MULTI-ZONE translation.

    This is the first executable BEMZ stage: each declared zone gets the
    same collocation operators as the validated single-zone Laplace solver,
    using that zone's material.  Interface coupling and the zone-specific
    ``BOUND`` equations remain a separate stage and are deliberately not
    hidden behind a single-zone solve.
    """
    from .assembly import assemble_laplace_domain

    if not model.materials:
        raise ValueError("the deck contains no ZONES material records")
    operators: dict[int, tuple[NDArray[np.complex128], NDArray[np.complex128]]] = {}
    zone_meshes = extract_zone_meshes(model)
    for zone_id, material in sorted(model.materials.items()):
        zone_mesh = zone_meshes.get(zone_id)
        if zone_mesh is None:
            continue
        operators[zone_id] = assemble_laplace_domain(
            zone_mesh, laplace_parameter, material, jump=jump,
        )
    return operators


def zone_boundary_table(model: MultiZoneModel) -> tuple[tuple[int, int, str], ...]:
    """Return ``(boundary, zone, role)`` records in BEMZ ordering.

    Boundary numbers are one-based.  A boundary listed as the first boundary
    of a zone is its outer interface; subsequent listed boundaries are inner
    interfaces.  Negative boundary IDs from the DAT ``ZONES`` record retain
    their free-surface meaning and are reported with the sign removed.
    """
    rows: list[tuple[int, int, str]] = []
    for boundary, zone in enumerate(model.boundary_zones, start=1):
        if int(zone) == 0:
            continue
        rows.append((boundary, int(zone), "zone interface"))
    return tuple(rows)


def zone_interface_pairs(model: MultiZoneModel) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Return the nested boundary lists needed by the translated ``GENRF`` stage.

    The result is ``(zone_number, (boundary_numbers, ...))``.  It is kept
    separate from matrix assembly because BEMZ applies continuity and
    equilibrium across these interfaces rather than treating every boundary
    as an independent single-zone problem.
    """
    grouped: dict[int, list[int]] = {}
    for boundary, zone in enumerate(model.boundary_zones, start=1):
        if int(zone):
            grouped.setdefault(int(zone), []).append(boundary)
    return tuple((zone, tuple(boundaries)) for zone, boundaries in sorted(grouped.items()))


def boundary_slice(model: MultiZoneModel, boundary_number: int) -> tuple[slice, slice]:
    """Return node and element slices for a one-based boundary."""
    if not model.boundary_node_ranges or not model.boundary_element_ranges:
        raise ValueError("model does not contain boundary ranges")
    index = boundary_number - 1
    if index < 0 or index >= len(model.boundary_node_ranges):
        raise IndexError(f"boundary number {boundary_number} is outside the model")
    node_start, node_end = model.boundary_node_ranges[index]
    elem_start, elem_end = model.boundary_element_ranges[index]
    return slice(node_start, node_end), slice(elem_start, elem_end)


def extract_boundary_mesh(model: MultiZoneModel, boundary_number: int) -> BoundaryMesh:
    """Extract one boundary with locally renumbered quadratic connectivity."""
    node_slice, element_slice = boundary_slice(model, boundary_number)
    coordinates = model.mesh.coordinates[node_slice].copy()
    connectivity = model.mesh.connectivity[element_slice].copy()
    connectivity -= node_slice.start
    if np.any(connectivity < 0) or np.any(connectivity >= len(coordinates)):
        raise ValueError("boundary connectivity crosses its declared node range")
    return BoundaryMesh(coordinates, connectivity)


def extract_zone_meshes(model: MultiZoneModel) -> dict[int, BoundaryMesh]:
    """Build independent meshes containing the boundaries assigned to each zone."""
    result: dict[int, BoundaryMesh] = {}
    for zone, boundaries in zone_interface_pairs(model):
        result[zone] = combine_meshes([extract_boundary_mesh(model, b) for b in boundaries])
    return result


def assemble_zone_time_step_operators(
    model: MultiZoneModel, step: int, total_time: float, step_count: int
) -> dict[int, tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Assemble independent BEMTIM impulse operators for each zone."""
    from .assembly import assemble_time_domain_step

    zone_meshes = extract_zone_meshes(model)
    return {
        zone_id: assemble_time_domain_step(zone_meshes[zone_id], step, total_time, step_count, material)
        for zone_id, material in sorted(model.materials.items())
        if zone_id in zone_meshes
    }


def block_diagonal_zone_operators(
    operators: dict[int, tuple[NDArray, NDArray]],
) -> tuple[NDArray, NDArray]:
    """Combine independent zone operators into deterministic block diagonals."""
    if not operators:
        raise ValueError("at least one zone operator is required")
    ordered = [operators[key] for key in sorted(operators)]
    size = sum(pair[0].shape[0] for pair in ordered)
    h = np.zeros((size, size), dtype=ordered[0][0].dtype)
    g = np.zeros((size, size), dtype=ordered[0][1].dtype)
    offset = 0
    for local_h, local_g in ordered:
        if local_h.shape != local_g.shape or local_h.shape[0] != local_h.shape[1]:
            raise ValueError("zone operators must be square matrices of equal shape")
        end = offset + local_h.shape[0]
        h[offset:end, offset:end] = local_h
        g[offset:end, offset:end] = local_g
        offset = end
    return h, g


def interface_block_shape(model: MultiZoneModel) -> tuple[int, int]:
    """Return the scalar row/column dimensions for the zone interface block."""
    meshes = extract_zone_meshes(model)
    rows = sum(2 * len(mesh.coordinates) for mesh in meshes.values())
    return rows, rows
