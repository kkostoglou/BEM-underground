"""Direct time-domain zone kernels and coupled convolution solver."""
import numpy as np
from .assembly import _legacy_quadrature
from .integration import integrate_elastostatic_element, integrate_time_element
from .multizone_system import assemble_coupled_zone_system


def evaluate_zone_time_kernels(coordinates, connectivity, step, total_time,
                               step_count, materials, progress=None):
    xy = np.asarray(coordinates, dtype=float)
    conn = np.asarray(connectivity)
    shape = (len(conn), len(xy), conn.shape[1], 3, 2, 2)
    h = np.zeros(shape); g = np.zeros(shape)
    diagonal = np.zeros((len(conn), len(xy), 2, 2))
    active = np.any(conn != 0, axis=2)
    for zone, material in enumerate(materials):
        poisson = material.lame_lambda / (2 * (material.lame_lambda + material.shear_modulus))
        active_elements = np.flatnonzero(active[zone])
        active_nodes = np.unique(conn[zone])
        active_nodes = active_nodes[active_nodes != 0]
        for node in active_nodes:
            if node == 0:
                continue
            point = xy[node - 1]
            for element_position, element in enumerate(active_elements):
                if progress is not None:
                    progress(zone + 1, len(materials), int(node), len(xy),
                             element_position + 1, len(active_elements))
                vertices = xy[conn[zone, element] - 1]
                points, subdivisions, singular = _legacy_quadrature(point, vertices)
                options = dict(quadrature_points=points, subdivisions=subdivisions,
                               singular_local_node=singular)
                displacement, traction = integrate_time_element(point, vertices,
                    step, total_time, step_count, material, **options)
                g[zone, node - 1, element] = displacement
                h[zone, node - 1, element] = traction
                if step == 1:
                    _, static = integrate_elastostatic_element(point, vertices,
                        poisson=poisson, shear_modulus=material.shear_modulus,
                        **options)
                    diagonal[zone, node - 1] -= static.sum(axis=0)
    return h, g, diagonal


def assemble_deck_multizone_time_systems(deck, progress=None):
    from .materials import Material
    from .multizone_preprocess import prepare_zone_connectivity
    coordinates, connectivity, types, global_types, inputs = prepare_zone_connectivity(deck)
    materials = []
    for item in inputs:
        shear = item.young / (2 * (1 + item.poisson))
        lame = 2 * shear * item.poisson / (1 - 2 * item.poisson)
        materials.append(Material(lame, shear, item.density / 1000))
    systems = []
    for step in range(1, deck.steps + 1):
        if progress is not None:
            progress(f"TD kernels: lag {step}/{deck.steps}")
        def kernel_progress(zone, zones, node, nodes, element, elements):
            if progress is not None:
                progress(f"TD kernels: lag {step}/{deck.steps}, zone {zone}/{zones}, "
                         f"node {node}/{nodes}, element {element}/{elements}")
        h, g, diagonal = evaluate_zone_time_kernels(coordinates, connectivity,
            step, deck.total_time, deck.steps, materials,
            progress=kernel_progress if progress is not None else None)
        if progress is not None:
            progress(f"TD matrix: lag {step}/{deck.steps}")
        systems.append(assemble_coupled_zone_system(connectivity, h, g, diagonal,
            types, global_types, [m.shear_modulus for m in materials],
            include_jump=step == 1))
    return coordinates, connectivity, types, systems


def solve_coupled_history(systems, prescribed_history, progress=None):
    """BEMTIM trapezoidal convolution for coupled multi-zone operators."""
    if not systems:
        raise ValueError("at least one impulse system is required")
    known = np.asarray(prescribed_history, dtype=float)
    first = systems[0]
    if known.shape != (len(systems), len(first.loading_nodes), 2):
        raise ValueError("prescribed history does not match loading nodes")
    if any((s.row_nodes, s.interface_nodes, s.loading_nodes) !=
           (first.row_nodes, first.interface_nodes, first.loading_nodes)
           for s in systems[1:]):
        raise ValueError("impulse systems use inconsistent ordering")
    known = known.reshape(len(systems), -1) / first.traction_scale
    dtype = np.result_type(first.matrix, first.loading_matrix, known)
    unknown = np.zeros((len(systems), first.matrix.shape[1]), dtype=dtype)
    zero_known = np.zeros(known.shape[1]); zero_unknown = np.zeros(unknown.shape[1])
    for current in range(len(systems)):
        if progress is not None:
            progress(f"TD history solve: step {current + 1}/{len(systems)}")
        right = np.zeros(first.matrix.shape[0], dtype=dtype)
        for lag in range(current + 1):
            newer = known[current - lag]
            older = known[current - lag - 1] if current - lag else zero_known
            right += systems[lag].loading_matrix @ ((newer + older) / 2)
        for lag in range(1, current + 1):
            newer = unknown[current - lag]
            older = unknown[current - lag - 1] if current - lag else zero_unknown
            right -= systems[lag].matrix @ ((newer + older) / 2)
        unknown[current] = np.linalg.solve(first.matrix, right)
    split = 2 * first.node_count
    return (np.real_if_close(unknown[:, :split]),
            np.real_if_close(unknown[:, split:] * first.traction_scale))
