"""Translate PREBEM ZONES declarations into BEMZ KZON/NZTYPE arrays."""
import numpy as np
from .mesh import mesh_section, combine_meshes


def prepare_zone_connectivity(deck):
    """Return coordinates, one-based KZON, NZTYPE, NDTYPE and materials.

    PREBEM requires every input zone boundary clockwise. For a finite zone its
    first (outer) boundary is reversed to anti-clockwise and classified 5;
    positive inner boundaries are 6, and negative inner voids are 7. If BOU1
    is zero, the zone extends to infinity: positive inner surfaces are type 8
    and negative free surfaces type 9. Region-4 ground surfaces retain type 2
    and the half-space jump regardless of their zone-boundary sign.
    """
    if not deck.zones:
        raise ValueError("deck contains no ZONES module")
    meshes = [mesh_section(section) for section in deck.sections]
    combined = combine_meshes(meshes)
    zone_count = len(deck.zones)
    element_count = len(combined.connectivity)
    node_count = len(combined.coordinates)
    kzon = np.zeros((zone_count, element_count, 3), dtype=np.int64)
    nztype = np.zeros((zone_count, node_count), dtype=np.int64)
    ndtype = np.zeros(node_count, dtype=np.int64)
    element_offset = node_offset = 0
    ranges = []
    for section, mesh in zip(deck.sections, meshes, strict=True):
        element_slice = slice(element_offset, element_offset + len(mesh.connectivity))
        node_slice = slice(node_offset, node_offset + len(mesh.coordinates))
        ranges.append((element_slice, node_slice))
        if section.region == 4:
            ndtype[node_slice] = 2
        element_offset = element_slice.stop
        node_offset = node_slice.stop
    for zone_index, zone in enumerate(deck.zones):
        for position, signed_boundary in enumerate(zone.boundaries):
            boundary = abs(signed_boundary)
            if not 1 <= boundary <= len(ranges):
                raise ValueError(f"zone {zone.number} refers to missing boundary {boundary}")
            element_slice, node_slice = ranges[boundary - 1]
            connectivity = combined.connectivity[element_slice] + 1
            finite_outer = position == 0 and zone.outer_boundary != 0
            if finite_outer:
                connectivity = connectivity[:, ::-1]
                kind = 5
            elif zone.outer_boundary == 0:
                kind = 9 if signed_boundary < 0 else 8
            else:
                kind = 7 if signed_boundary < 0 else 6
            if deck.sections[boundary - 1].region == 4:
                kind = 2
            if np.any(kzon[zone_index, element_slice]):
                raise ValueError(f"boundary {boundary} is repeated within zone {zone.number}")
            kzon[zone_index, element_slice] = connectivity
            nztype[zone_index, node_slice] = kind
            if kind in (2, 7, 9):
                ndtype[node_slice] = kind
    materials = tuple(zone.material for zone in deck.zones)
    return combined.coordinates, kzon, nztype, ndtype, materials
