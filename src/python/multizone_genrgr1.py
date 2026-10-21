"""Free-boundary loading block from GENRGR1 in BEMZ3.FOR.

Column keys explicitly retain (zone index, Fortran node number): the legacy
NDF counter continues across zones. Kernels are unscaled scratch-file-7
displacement kernels, already evaluated in KZON orientation.
"""
import numpy as np


def assemble_genrgr1(connectivity, displacement_kernels, node_types,
                     global_node_types, shear_moduli):
    """Return loading matrix, collocation row keys, and loading column keys.

    Matrix multiplies prescribed tractions divided by max(shear_moduli).
    For a single zone the original code retains all global node columns.
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
    rows = tuple((z, int(n)) for z in range(zones) for n in np.unique(conn[z]) if n)
    free = (2, 7, 9)
    if zones == 1:
        keys = tuple((0, n) for n in range(1, nodes + 1))
    else:
        for z, n in rows:
            if (types[z, n - 1] in free) != (global_types[n - 1] in free):
                raise ValueError("zone and global free-boundary classifications conflict")
        keys = tuple((z, n) for z, n in rows if global_types[n - 1] in free)
    columns = {key: i for i, key in enumerate(keys)}
    matrix = np.zeros((2 * len(rows), 2 * len(keys)), dtype=np.result_type(kernels, complex))
    for index, (zone, node) in enumerate(rows):
        row = slice(2 * index, 2 * index + 2)
        for element, connection in enumerate(conn[zone]):
            for local, target in enumerate(connection):
                if target == 0 or (zones > 1 and types[zone, target - 1] not in free):
                    continue
                column = 2 * columns[zone, int(target)]
                matrix[row, column:column + 2] += mu.max() * kernels[zone, node - 1, element, local]
    return matrix, rows, keys
