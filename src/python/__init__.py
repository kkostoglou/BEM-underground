"""Python translation of the FinalNov93 Boundary Element package."""

from .assembly import assemble_elastostatic, assemble_laplace_domain, assemble_laplace_internal, assemble_time_domain_step
from .geometry import element_length, shape_functions
from .fortran_io import (
    read_complex_records, read_pro_displacements, read_pro_real_displacements,
    read_real_records,
)
from .input_deck import CompoundSegment, ExplosionPoint, InputDeck, MaterialInput, PlaneWave, Section, read_input_deck
from .integration import integrate_elastostatic_element, integrate_laplace_element, integrate_time_element
from .materials import Material, viscoelastic_properties
from .loads import (
    explosion_stress_spectrum, explosion_stress_time, explosion_tractions,
    explosion_tractions_time, boundary_material_normals, longitudinal_plane_wave_tractions,
    longitudinal_plane_wave_tractions_time,
    stress_tractions,
)
from .mesh import (
    BoundaryMesh, boundary_node_types, combine_meshes, mesh_circle, mesh_compound, mesh_ellipse, mesh_rectangle,
    mesh_section, mesh_straight_line,
)
from .numerics import bessel_k012_legacy, elastostatic_kernel, laplace_domain_kernel
from .solver import solve_mixed_boundary_system, solve_mixed_history, solve_traction_history
from .stress import stress_global, stress_local, traction
from .reflection import (
    ReflectedPStress, delayed_p_reflected_stress, free_surface_line, p_reflected_stress_components, p_reflection_path_delays, p_wave_amplitudes, p_wave_stress_coefficients, reflection_condition,
    sv_wave_amplitudes, sv_wave_stress_coefficients,
)
from .plane_reflection import plane_reflected_stress, reflected_ray_times, sv_reflected_stress_components
from .transforms import LaplaceGrid, inverse_laplace_spectrum, laplace_grid
from .time_domain import time_domain_kernel
from .workflow import MultiZoneModel, SolverResult, prepare_multizone_model, run_laplace_case, run_time_case, save_result
from .multizone_workflow import MultiZoneResult, run_multizone_laplace_case, run_multizone_time_case
from .multizone_cq import solve_coupled_convolution_quadrature
from .multizone import (
    MultiZoneBlocks, allocate_block_system, assemble_zone_local_operators,
    assemble_zone_time_step_operators, block_diagonal_zone_operators, build_block_layout, boundary_slice, extract_boundary_mesh, extract_zone_meshes, interface_block_shape, zone_boundary_table, zone_interface_pairs,
)

__all__ = [
    "Material",
    "MaterialInput",
    "CompoundSegment",
    "PlaneWave",
    "ExplosionPoint",
    "Section",
    "InputDeck",
    "BoundaryMesh",
    "boundary_node_types",
    "combine_meshes",
    "assemble_elastostatic",
    "assemble_laplace_domain",
    "assemble_laplace_internal",
    "assemble_time_domain_step",
    "bessel_k012_legacy",
    "elastostatic_kernel",
    "laplace_domain_kernel",
    "solve_mixed_boundary_system",
    "solve_mixed_history",
    "solve_traction_history",
    "element_length",
    "shape_functions",
    "read_input_deck",
    "mesh_circle",
    "mesh_compound",
    "mesh_ellipse",
    "mesh_rectangle",
    "mesh_straight_line",
    "mesh_section",
    "integrate_elastostatic_element",
    "integrate_laplace_element",
    "integrate_time_element",
    "stress_global",
    "stress_local",
    "traction",
    "reflection_condition",
    "ReflectedPStress",
    "delayed_p_reflected_stress",
    "free_surface_line",
    "p_reflection_path_delays",
    "p_reflected_stress_components",
    "p_wave_amplitudes",
    "sv_wave_amplitudes",
    "p_wave_stress_coefficients",
    "sv_wave_stress_coefficients",
    "plane_reflected_stress",
    "reflected_ray_times",
    "sv_reflected_stress_components",
    "viscoelastic_properties",
    "LaplaceGrid",
    "laplace_grid",
    "inverse_laplace_spectrum",
    "read_complex_records",
    "read_pro_displacements",
    "read_real_records",
    "read_pro_real_displacements",
    "longitudinal_plane_wave_tractions",
    "boundary_material_normals",
    "stress_tractions",
    "longitudinal_plane_wave_tractions_time",
    "explosion_stress_spectrum",
    "explosion_stress_time",
    "explosion_tractions",
    "explosion_tractions_time",
    "time_domain_kernel",
    "SolverResult",
    "run_laplace_case",
    "run_time_case",
    "save_result",
    "MultiZoneModel",
    "prepare_multizone_model",
    "MultiZoneResult",
    "run_multizone_laplace_case",
    "run_multizone_time_case",
    "solve_coupled_convolution_quadrature",
    "MultiZoneBlocks",
    "build_block_layout",
    "allocate_block_system",
    "assemble_zone_local_operators",
    "zone_boundary_table",
    "zone_interface_pairs",
    "boundary_slice",
    "extract_boundary_mesh",
    "extract_zone_meshes",
    "assemble_zone_time_step_operators",
    "block_diagonal_zone_operators",
    "interface_block_shape",
]
