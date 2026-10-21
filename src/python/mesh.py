"""Boundary mesh construction translated from the PREGEOM routines."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

from .input_deck import CompoundSegment, Section


@dataclass(frozen=True)
class BoundaryMesh:
    coordinates: NDArray[np.float64]
    connectivity: NDArray[np.int64]


def mesh_circle(section: Section) -> BoundaryMesh:
    """Create PREBEM-compatible quadratic elements for a circle.

    Connectivity is returned zero-based. Its columns retain the Fortran order
    ``first, middle, last``.
    """

    if section.shape != 1:
        raise ValueError("mesh_circle requires shape code 1")
    radius = section.geometry[0]
    ratio = section.geometry[1] if len(section.geometry) > 1 else 1.0
    center_x = section.geometry[2] if len(section.geometry) > 2 else 0.0
    center_y = section.geometry[3] if len(section.geometry) > 3 else 0.0
    rotation = section.geometry[4] if len(section.geometry) > 4 else 0.0
    if radius <= 0.0:
        raise ValueError("circle radius must be positive")
    if not math.isclose(ratio, 1.0):
        return mesh_ellipse(section)

    node_count = 2 * section.x_elements
    angles = math.radians(rotation) - np.arange(node_count) * (2.0 * math.pi / node_count)
    coordinates = np.column_stack(
        (center_x + radius * np.cos(angles), center_y + radius * np.sin(angles))
    )
    element = np.arange(section.x_elements, dtype=np.int64)
    middle = 2 * element
    connectivity = np.column_stack(((middle - 1) % node_count, middle, (middle + 1) % node_count))
    return BoundaryMesh(coordinates=coordinates, connectivity=connectivity)


def _closed_quadratic_mesh(coordinates: NDArray[np.float64]) -> BoundaryMesh:
    node_count = len(coordinates)
    if node_count % 2:
        raise ValueError("a closed quadratic mesh requires an even node count")
    middle = 2 * np.arange(node_count // 2, dtype=np.int64)
    connectivity = np.column_stack(((middle - 1) % node_count, middle, (middle + 1) % node_count))
    return BoundaryMesh(coordinates, connectivity)


def mesh_ellipse(section: Section) -> BoundaryMesh:
    """Mesh shape 1 with a non-unit second-semiaxis ratio."""

    if section.shape != 1 or len(section.geometry) < 2:
        raise ValueError("ellipse geometry requires semiaxis and ratio")
    major, ratio = section.geometry[:2]
    center_x = section.geometry[2] if len(section.geometry) > 2 else 0.0
    center_y = section.geometry[3] if len(section.geometry) > 3 else 0.0
    rotation = np.deg2rad(section.geometry[4] if len(section.geometry) > 4 else 0.0)
    if major <= 0.0 or ratio <= 0.0:
        raise ValueError("ellipse semiaxes must be positive")
    angles = -2.0 * np.pi * np.arange(2 * section.x_elements) / (2 * section.x_elements)
    local = np.column_stack((major * np.cos(angles), major * ratio * np.sin(angles)))
    transform = np.array([[np.cos(rotation), -np.sin(rotation)], [np.sin(rotation), np.cos(rotation)]])
    return _closed_quadratic_mesh(local @ transform.T + [center_x, center_y])


def mesh_rectangle(section: Section) -> BoundaryMesh:
    """Uniform quadratic square/rectangle boundary for shape codes 2 and 3."""

    if section.shape not in (2, 3):
        raise ValueError("mesh_rectangle requires shape code 2 or 3")
    if section.shape == 2:
        width = height = section.geometry[0]; nx = ny = section.x_elements
        offset = 1
    else:
        width, height = section.geometry[:2]; nx = section.x_elements
        ny = section.y_elements or section.x_elements; offset = 2
    center_x = section.geometry[offset] if len(section.geometry) > offset else 0.0
    center_y = section.geometry[offset + 1] if len(section.geometry) > offset + 1 else 0.0
    rotation = np.deg2rad(section.geometry[offset + 2] if len(section.geometry) > offset + 2 else 0.0)
    if min(width, height, nx, ny) <= 0:
        raise ValueError("rectangle dimensions and element counts must be positive")
    # MESH4 starts at the upper-right corner and traverses anticlockwise.
    px = section.geometry[offset + 3] if len(section.geometry) > offset + 3 else 0.0
    py = section.geometry[offset + 4] if len(section.geometry) > offset + 4 else 0.0
    # Sharp corners have two one-sided normals and stress values. Preserve
    # the original node numbers and append the other side of each corner.
    # An explicit legacy flag of 1 still requests shared corner nodes.
    corner_type = section.geometry[offset + 5] if len(section.geometry) > offset + 5 else 1
    discontinuous = corner_type == 2
    if px == py == 0:
        px = py = 1.0
    elif px == 0 or py == 0:
        # RATENUM matches the first half-element length of the other side.
        from scipy.optimize import brentq
        if px == 0:
            target = height * graded_line_nodes(ny, py)[1] / width
            px = _matching_grade(nx, target, py, brentq)
        else:
            target = width * graded_line_nodes(nx, px)[1] / height
            py = _matching_grade(ny, target, px, brentq)
    corners = np.array([[width/2, height/2], [-width/2, height/2],
                        [-width/2, -height/2], [width/2, -height/2]])
    counts = (nx, ny, nx, ny)
    pieces = []
    for i, (count, grade) in enumerate(zip(counts, (px, py, px, py), strict=True)):
        fractions = graded_line_nodes(count, grade)[:-1]
        pieces.append(corners[i] + fractions[:, None] * (corners[(i+1) % 4] - corners[i]))
    local = np.vstack(pieces)
    start = 4 * nx + 3 * ny
    sign = 1 if section.region in (2, 5) else -1
    ordering = (start + sign * np.arange(len(local))) % len(local)
    local = local[ordering]
    transform = np.array([[np.cos(rotation), -np.sin(rotation)], [np.sin(rotation), np.cos(rotation)]])
    coordinates = local @ transform.T + [center_x, center_y]
    if ny % 2:
        connectivity = _closed_quadratic_mesh(coordinates).connectivity
    else:
        first = 2 * np.arange(len(coordinates) // 2, dtype=np.int64)
        connectivity = np.column_stack((first, first + 1, (first + 2) % len(coordinates)))
    if discontinuous:
        original_corners = (2*nx, 2*(nx+ny), 4*nx+2*ny, 0)
        for original in original_corners:
            node = int(np.flatnonzero(ordering == original)[0])
            attached = np.argwhere(connectivity == node)
            if len(attached) != 2:
                raise ValueError("a discontinuous corner must join two element endpoints")
            element, slot = attached[1]
            connectivity[element, slot] = len(coordinates)
            coordinates = np.vstack((coordinates, coordinates[node]))
    return BoundaryMesh(coordinates, connectivity)


def graded_line_nodes(elements: int, ratio: float = 1.0) -> NDArray[np.float64]:
    """SYNTEL positions on [0, 1], including quadratic middles and both ends.

    Positive ratios grade from one end; negative ratios grade symmetrically
    from both ends. Each quadratic middle bisects its own element.
    """
    if elements <= 0 or ratio == 0 or not np.isfinite(ratio):
        raise ValueError("positive element count and finite nonzero grading ratio required")
    exponents = np.arange(elements, dtype=float)
    if ratio < 0:
        exponents = np.minimum(exponents, exponents[::-1])
    logarithms = exponents * np.log(abs(ratio))
    widths = np.exp(logarithms - np.max(logarithms))
    # Prevent pathological corner elements.  The Fortran generator limits
    # refinement through its integer search; an unconstrained geometric
    # progression here can otherwise produce elements many orders smaller
    # than the side and destabilize HOOP stress recovery.
    minimum_fraction = 1.0e-3
    widths = np.maximum(widths, minimum_fraction)
    widths /= widths.sum()
    endpoints = np.r_[0., np.cumsum(widths)]
    result = np.empty(2 * elements + 1)
    result[::2] = endpoints
    result[1::2] = (endpoints[:-1] + endpoints[1:]) / 2
    result[-1] = 1.0
    return result


def _matching_grade(elements, target, other_grade, root):
    """Numerically solve RATENUM's equal corner-spacing condition."""
    sign = -1 if other_grade < 0 else 1
    function = lambda logarithm: graded_line_nodes(elements, sign * math.exp(logarithm))[1] - target
    if not 0 < target < 0.5:
        raise ValueError("automatic grading cannot match the requested corner spacing")
    try:
        return sign * math.exp(root(function, -20., 20.))
    except ValueError as error:
        raise ValueError("automatic grading has no solution for these dimensions") from error



def mesh_straight_line(section: Section) -> BoundaryMesh:
    """Uniform open quadratic line mesh for shape code 4."""

    if section.shape != 4 or len(section.geometry) < 4:
        raise ValueError("straight-line geometry requires x0, y0, x1, y1")
    start = np.array(section.geometry[:2]); end = np.array(section.geometry[2:4])
    if np.linalg.norm(end - start) < 1.0e-11:
        raise ValueError("straight-line endpoints must differ")
    ratio = section.geometry[4] if len(section.geometry) > 4 else 1.0
    fractions = graded_line_nodes(section.x_elements, ratio or 1.0)
    coordinates = start + fractions[:, None] * (end - start)
    first = 2 * np.arange(section.x_elements, dtype=np.int64)
    connectivity = np.column_stack((first, first + 1, first + 2))
    return BoundaryMesh(coordinates, connectivity)


def _compound_segment_coordinates(segment: CompoundSegment) -> NDArray[np.float64]:
    """Create all nodes, including the terminal node, for ``GEOCOMB``."""
    if segment.shape == 4:
        if len(segment.geometry) < 4:
            raise ValueError("compound straight segment needs x0, y0, x1, y1")
        return np.linspace(segment.geometry[:2], segment.geometry[2:4], 2 * segment.elements + 1)
    if segment.shape == 1:
        if len(segment.geometry) < 7:
            raise ValueError("compound ellipse segment needs center, angles, semiaxis, ratio, rotation")
        x0, y0, start, finish, radius, ratio, rotation = segment.geometry[:7]
        if radius <= 0.0 or ratio <= 0.0:
            raise ValueError("compound ellipse semiaxes must be positive")
        angle = np.linspace(math.radians(start), math.radians(finish), 2 * segment.elements + 1)
        local = np.column_stack((radius * np.cos(angle), radius * ratio * np.sin(angle)))
        theta = math.radians(rotation)
        transform = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
        return local @ transform.T + [x0, y0]
    raise ValueError("compound segment shape must be 1 or 4")


def mesh_compound(section: Section) -> BoundaryMesh:
    """Translate the quadratic-node construction in PREBEM ``GEOCOMB``."""
    if section.shape >= 0 or not section.segments:
        raise ValueError("mesh_compound requires a negative shape with segments")
    pieces = [_compound_segment_coordinates(segment) for segment in section.segments]
    for previous, following in zip(pieces, pieces[1:], strict=False):
        if not np.allclose(previous[-1], following[0], atol=1.0e-9):
            raise ValueError("compound segments must join end-to-start")
    coordinates = np.vstack([piece[:-1] for piece in pieces])
    is_open = section.region in (3, 4, 5)
    if is_open:
        coordinates = np.vstack((coordinates, pieces[-1][-1]))
        first = 2 * np.arange(section.x_elements, dtype=np.int64)
        connectivity = np.column_stack((first, first + 1, first + 2))
    else:
        if not np.allclose(pieces[-1][-1], pieces[0][0], atol=1.0e-9):
            raise ValueError("closed compound boundaries must return to their first point")
        connectivity = _closed_quadratic_mesh(coordinates).connectivity
    return BoundaryMesh(coordinates, connectivity)


def mesh_section(section: Section) -> BoundaryMesh:
    """Dispatch a simple PREBEM section to its translated mesh generator."""

    if section.shape == 1:
        ratio = section.geometry[1] if len(section.geometry) > 1 else 1.0
        return mesh_circle(section) if np.isclose(ratio, 1.0) else mesh_ellipse(section)
    if section.shape in (2, 3): return mesh_rectangle(section)
    if section.shape == 4: return mesh_straight_line(section)
    if section.shape < 0:
        return mesh_compound(section)
    raise ValueError(f"unknown PREBEM shape code {section.shape}")


def combine_meshes(meshes: tuple[BoundaryMesh, ...] | list[BoundaryMesh]) -> BoundaryMesh:
    """Combine independent boundary meshes while preserving element ordering."""

    if not meshes:
        raise ValueError("at least one boundary mesh is required")
    coordinates: list[NDArray[np.float64]] = []
    connectivity: list[NDArray[np.int64]] = []
    offset = 0
    for mesh in meshes:
        if mesh.coordinates.ndim != 2 or mesh.coordinates.shape[1] != 2:
            raise ValueError("mesh coordinates must have shape (nodes, 2)")
        if mesh.connectivity.ndim != 2 or mesh.connectivity.shape[1] != 3:
            raise ValueError("quadratic connectivity must have shape (elements, 3)")
        coordinates.append(mesh.coordinates)
        connectivity.append(mesh.connectivity + offset)
        offset += len(mesh.coordinates)
    return BoundaryMesh(np.vstack(coordinates), np.vstack(connectivity))


def boundary_node_types(sections: tuple[Section, ...] | list[Section], meshes: tuple[BoundaryMesh, ...] | list[BoundaryMesh], *, simulation_boundaries: tuple[int, ...] = ()) -> NDArray[np.int64]:
    """Translate PREBEM ``FREESURF`` node classifications.

    The returned codes are 0 ordinary, 1 inclusion, 2 free surface, 3
    simulation surface, and 4 enclosing boundary.  Boundary numbers in
    ``simulation_boundaries`` are one-based, as in the original input.
    """
    if len(sections) != len(meshes):
        raise ValueError("sections and meshes must have the same length")
    types: list[NDArray[np.int64]] = []
    simulated = set(simulation_boundaries)
    for number, (section, mesh) in enumerate(zip(sections, meshes, strict=True), start=1):
        kind = 0
        if section.region == 2:
            kind = 1
        if section.region == 4:
            kind = 2
        if number in simulated:
            kind = 3
        if section.region in (5, 6):
            kind = 4
        types.append(np.full(len(mesh.coordinates), kind, dtype=np.int64))
    return np.concatenate(types) if types else np.empty(0, dtype=np.int64)


def validate_zone_boundaries(section_count: int, zones: object) -> None:
    """Validate MULTI-ZONE boundary references before numerical assembly."""
    for zone in zones:
        boundaries = getattr(zone, "boundaries", ()) or tuple(
            b for b in (getattr(zone, "outer_boundary", 0), getattr(zone, "inner_boundary", 0)) if b
        )
        for boundary in boundaries:
            if abs(int(boundary)) > section_count:
                raise ValueError(f"zone {getattr(zone, 'number', '?')} references boundary {boundary}, but only {section_count} exist")


def zone_ids_per_boundary(section_count: int, zones: object) -> NDArray[np.int64]:
    """Return the one-based zone assignment for each boundary (zero if none)."""
    validate_zone_boundaries(section_count, zones)
    result = np.zeros(section_count, dtype=np.int64)
    for zone in zones:
        for boundary in (getattr(zone, "boundaries", ()) or ()):
            index = abs(int(boundary)) - 1
            if result[index] not in (0, int(zone.number)):
                raise ValueError(f"boundary {abs(int(boundary))} is assigned to multiple zones")
            result[index] = int(zone.number)
    return result
