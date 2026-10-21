"""Interface traction block of GENRGL in BEMZ3.FOR (lines 742--1020)."""
import numpy as np


def assemble_genrgl(connectivity, displacement_kernels, node_types,
                    global_node_types, shear_moduli):
    """Return traction block, row ordering, and one-based column node IDs.

    Connectivity uses Fortran numbering (zero means absent). Kernels have
    shape (zones,nodes,elements,3,2,2), in KZON orientation, and must be
    unscaled displacement kernels as stored in scratch file 7. The returned
    block follows the 2*NNODE displacement columns of GENRF. Its unknowns
    are tractions divided by max(shear_moduli).
    """
    conn = np.asarray(connectivity)
    kernels = np.asarray(displacement_kernels)
    types = np.asarray(node_types)
    global_types = np.asarray(global_node_types)
    mu = np.asarray(shear_moduli, dtype=float)
    if types.ndim != 2 or conn.ndim != 3 or conn.shape[2] != 3:
        raise ValueError("expected zone/node types and quadratic connectivity")
    zones, nodes = types.shape
    if zones == 0 or conn.shape[0] != zones:
        raise ValueError("zone counts must agree and be nonzero")
    if not np.issubdtype(conn.dtype, np.integer) or np.any(conn < 0) or np.any(conn > nodes):
        raise ValueError("invalid Fortran node numbers")
    if global_types.shape != (nodes,):
        raise ValueError("global node types must have one entry per node")
    if mu.shape != (zones,) or not np.all(np.isfinite(mu)) or np.any(mu <= 0):
        raise ValueError("one positive finite shear modulus required per zone")
    if kernels.shape != (zones, nodes, conn.shape[1], 3, 2, 2):
        raise ValueError("incorrect displacement kernel shape")
    free = (2, 7, 9)
    column_nodes = tuple(int(n + 1) for n in np.flatnonzero(~np.isin(global_types, free)))
    columns = {node: i for i, node in enumerate(column_nodes)}
    rows = tuple((z, int(n)) for z in range(zones) for n in np.unique(conn[z]) if n)
    matrix = np.zeros((2 * len(rows), 2 * len(columns)), dtype=np.result_type(kernels, complex))
    scale = float(mu.max())
    for index, (zone, node) in enumerate(rows):
        row = slice(2 * index, 2 * index + 2)
        for element, connection in enumerate(conn[zone]):
            for local, target in enumerate(connection):
                if target == 0 or types[zone, target - 1] in free:
                    continue
                if target not in columns:
                    raise ValueError("zone interface type conflicts with global free-boundary type")
                column = 2 * columns[target]
                sign = -1 if types[zone, target - 1] in (6, 8) else 1
                matrix[row, column:column + 2] += sign * scale * kernels[zone, node - 1, element, local]
    return matrix, rows, column_nodes
