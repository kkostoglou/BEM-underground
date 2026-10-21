"""GENRF displacement block from BEMZ3.FOR, lines 530--734.

Inputs retain Fortran node numbering: zero connectivity means absent.
Kernels must already use KZON orientation, as in the original scratch file 8.
This routine assembles coefficients; it does not compute or solve the kernels.
"""
import numpy as np


def assemble_genrf(connectivity, traction_kernels, static_diagonal, node_types,
                   *, include_jump=True):
    """Return F and its (zero-based zone, one-based node) row-node ordering.

    Shapes: connectivity (zones,elements,3), kernels
    (zones,nodes,elements,3,2,2), diagonal (zones,nodes,2,2),
    types (zones,nodes). Columns are global node X,Y displacements.
    """
    connectivity = np.asarray(connectivity)
    types = np.asarray(node_types)
    kernels = np.asarray(traction_kernels)
    diagonal = np.asarray(static_diagonal)
    if types.ndim != 2 or connectivity.ndim != 3 or connectivity.shape[2] != 3:
        raise ValueError("expected zone/node types and quadratic zone connectivity")
    zones, nodes = types.shape
    if connectivity.shape[0] != zones or zones == 0:
        raise ValueError("zone counts must agree and be nonzero")
    if not np.issubdtype(connectivity.dtype, np.integer):
        raise ValueError("connectivity must contain integer Fortran node numbers")
    if np.any(connectivity < 0) or np.any(connectivity > nodes):
        raise ValueError("connectivity node outside global node range")
    elements = connectivity.shape[1]
    if kernels.shape != (zones, nodes, elements, 3, 2, 2):
        raise ValueError("incorrect traction kernel shape")
    if diagonal.shape != (zones, nodes, 2, 2):
        raise ValueError("incorrect static diagonal shape")
    ordering = tuple((z, int(n)) for z in range(zones)
                     for n in np.unique(connectivity[z]) if n > 0)
    result = np.zeros((2 * len(ordering), 2 * nodes),
                      dtype=np.result_type(kernels, diagonal, np.complex128))
    jumps = {5: 0., 7: 0., 6: 1., 8: 1., 9: 0., 2: .5}
    for index, (zone, node) in enumerate(ordering):
        row = slice(2 * index, 2 * index + 2)
        kind = int(types[zone, node - 1])
        # Preserve the override order in GENRF, including NZONE == 1.
        if kind in (8, 9, 2):
            jump = jumps[kind]
        elif zones == 1:
            jump = 1.
        elif kind in jumps:
            jump = jumps[kind]
        else:
            raise ValueError(f"undefined GENRF jump for zone {zone + 1}, node {node}, type {kind}")
        for element, connection in enumerate(connectivity[zone]):
            for local, target in enumerate(connection):
                if target:
                    column = slice(2 * (target - 1), 2 * target)
                    result[row, column] += kernels[zone, node - 1, element, local]
        column = slice(2 * (node - 1), 2 * node)
        result[row, column] += diagonal[zone, node - 1]
        if include_jump:
            result[row, column] += jump * np.eye(2)
    return result, ordering
