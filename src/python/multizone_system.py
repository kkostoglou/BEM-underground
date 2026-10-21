"""Coupled PREGENR matrix assembly from already evaluated BEMZ kernels.

This layer requires original-style oriented connectivity and node types;
it does not infer them from the earlier experimental zone mesh helpers.
"""
from dataclasses import dataclass
import numpy as np
from .multizone_genrf import assemble_genrf
from .multizone_genrgl import assemble_genrgl
from .multizone_genrgr1 import assemble_genrgr1


@dataclass(frozen=True)
class CoupledZoneSystem:
    matrix: np.ndarray
    loading_matrix: np.ndarray
    row_nodes: tuple
    interface_nodes: tuple
    loading_nodes: tuple
    node_count: int
    traction_scale: float

    def solve_tractions(self, prescribed):
        """Solve A x = B (t/GSCALE) for prescribed loading-column tractions.

        Returns global (nodes,2) displacements and (interfaces,2) physical
        tractions in interface_nodes order. No incident-field conversion or
        prescribed-displacement column exchange is performed here.
        """
        traction = np.asarray(prescribed)
        if traction.shape != (len(self.loading_nodes), 2):
            raise ValueError("prescribed tractions must follow loading_nodes with two components")
        if not np.all(np.isfinite(traction)):
            raise ValueError("prescribed tractions must be finite")
        rhs = self.loading_matrix @ (traction.ravel() / self.traction_scale)
        solution = np.linalg.solve(self.matrix, rhs)
        split = 2 * self.node_count
        return (solution[:split].reshape(-1, 2),
                solution[split:].reshape(-1, 2) * self.traction_scale)


def assemble_coupled_zone_system(connectivity, traction_kernels,
                                displacement_kernels, static_diagonal,
                                node_types, global_node_types, shear_moduli,
                                *, include_jump=True):
    """Assemble PREGENR's [F G_interface] and B loading matrix."""
    f, rows = assemble_genrf(connectivity, traction_kernels, static_diagonal,
                             node_types, include_jump=include_jump)
    b, load_rows, loading = assemble_genrgr1(connectivity, displacement_kernels,
        node_types, global_node_types, shear_moduli)
    if len(connectivity) == 1:  # PREGENR skips GENRGL for NZONE=1.
        matrix, interfaces = f, ()
    else:
        g, interface_rows, interfaces = assemble_genrgl(connectivity,
            displacement_kernels, node_types, global_node_types, shear_moduli)
        if interface_rows != rows:
            raise ValueError("interface row ordering disagrees with GENRF")
        matrix = np.hstack((f, g))
    if load_rows != rows:
        raise ValueError("loading row ordering disagrees with GENRF")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("zone membership does not produce a square coupled system")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(b)):
        raise ValueError("assembled operators must be finite")
    return CoupledZoneSystem(matrix, b, rows, interfaces, loading,
                             len(global_node_types), float(np.max(shear_moduli)))


def assemble_deck_multizone_system(deck, laplace_parameter):
    """Preprocess a ZONES deck and assemble its coupled Laplace system."""
    from .materials import Material
    from .multizone_kernels import evaluate_zone_laplace_kernels
    from .multizone_preprocess import prepare_zone_connectivity

    coordinates, connectivity, node_types, global_types, inputs = prepare_zone_connectivity(deck)
    materials = []
    for item in inputs:
        shear = item.young / (2 * (1 + item.poisson))
        lame = 2 * shear * item.poisson / (1 - 2 * item.poisson)
        materials.append(Material(lame, shear, item.density / 1000,
                                  item.viscosity, item.viscosity_model))
    h, g, diagonal = evaluate_zone_laplace_kernels(
        coordinates, connectivity, laplace_parameter, materials)
    system = assemble_coupled_zone_system(connectivity, h, g, diagonal,
        node_types, global_types, [item.shear_modulus for item in materials])
    return coordinates, connectivity, node_types, system
