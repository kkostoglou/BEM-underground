"""End-to-end workflows for the validated reference problem families."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .assembly import (
    assemble_laplace_domain, assemble_laplace_internal,
    assemble_time_domain_step, assemble_time_internal,
)
from .input_deck import ExplosionPoint, InputDeck, PlaneWave, SurfaceRamp, TractionNodes, read_input_deck
from .loads import (
    explosion_stress_spectrum, explosion_stress_time, explosion_tractions,
    explosion_tractions_time, longitudinal_plane_wave_tractions,
    longitudinal_plane_wave_tractions_time, plane_wave_stress_spectrum,
    plane_wave_stress_time, boundary_material_normals, stress_tractions,
)
from .hoop import boundary_stress, incident_local_stress, nodal_stress
from .materials import Material
from .plane_reflection import plane_reflected_stress
from .mesh import BoundaryMesh, boundary_node_types, combine_meshes, mesh_section, validate_zone_boundaries, zone_ids_per_boundary
from .solver import solve_mixed_boundary_system, solve_mixed_history, solve_traction_history
from .transforms import inverse_laplace_spectrum, laplace_grid


@dataclass(frozen=True)
class SolverResult:
    time: NDArray[np.float64]
    coordinates: NDArray[np.float64]
    connectivity: NDArray[np.int64]
    displacement: NDArray[np.float64]
    traction: NDArray[np.float64]
    displacement_spectrum: NDArray[np.complex128] | None = None
    traction_spectrum: NDArray[np.complex128] | None = None
    boundary_stress: NDArray[np.float64] | None = None
    element_stress: NDArray[np.float64] | None = None
    stress_spectrum: NDArray[np.complex128] | None = None
    internal_coordinates: NDArray[np.float64] | None = None
    internal_displacement: NDArray[np.float64] | None = None
    internal_displacement_spectrum: NDArray[np.complex128] | None = None


def _model(deck: InputDeck) -> tuple[BoundaryMesh, Material, NDArray[np.int64]]:
    validate_zone_boundaries(len(deck.sections), deck.zones)
    meshes = [mesh_section(section) for section in deck.sections]
    mesh = combine_meshes(meshes)
    # The thesis ramp models the material below an open line.  Its input line
    # is listed left-to-right, which makes CSHAPE point the element normal into
    # the material.  The exterior half-space equation requires the opposite
    # normal.  Reverse local connectivity (and therefore normals) while
    # retaining the original node numbers used by plots and boundary data.
    if (
        isinstance(deck.load, SurfaceRamp)
        and len(deck.sections) == 1
        and deck.sections[0].region in (3, 4)
    ):
        mesh = BoundaryMesh(mesh.coordinates, mesh.connectivity[:, ::-1].copy())
    shear = deck.material.young / (2.0 * (1.0 + deck.material.poisson))
    lame = 2.0 * shear * deck.material.poisson / (1.0 - 2.0 * deck.material.poisson)
    simulation = () if deck.simulation is None else (deck.simulation,)
    return mesh, Material(lame, shear, deck.material.density / 1000.0, deck.material.viscosity, deck.material.viscosity_model), boundary_node_types(deck.sections, meshes, simulation_boundaries=simulation)


def zone_materials(deck: InputDeck) -> dict[int, Material]:
    """Build solver material objects for each parsed MULTI-ZONE record."""
    materials: dict[int, Material] = {}
    for zone in deck.zones:
        mat = zone.material
        shear = mat.young / (2.0 * (1.0 + mat.poisson))
        lame = 2.0 * shear * mat.poisson / (1.0 - 2.0 * mat.poisson)
        materials[zone.number] = Material(lame, shear, mat.density / 1000.0, mat.viscosity, mat.viscosity_model)
    return materials


def zone_boundary_assignments(deck: InputDeck) -> NDArray[np.int64]:
    """Return zone IDs aligned with the deck's one-based boundary numbers."""
    return zone_ids_per_boundary(len(deck.sections), deck.zones)


@dataclass(frozen=True)
class MultiZoneModel:
    """Prepared geometry and material data for multi-zone matrix assembly."""

    mesh: BoundaryMesh
    boundary_zones: NDArray[np.int64]
    materials: dict[int, Material]
    node_types: NDArray[np.int64]
    boundary_node_ranges: tuple[tuple[int, int], ...] = ()
    boundary_element_ranges: tuple[tuple[int, int], ...] = ()


def prepare_multizone_model(deck: InputDeck) -> MultiZoneModel:
    """Prepare zone-aware inputs without changing the existing solver path."""
    meshes = [mesh_section(section) for section in deck.sections]
    node_ranges = []
    element_ranges = []
    node_offset = element_offset = 0
    for item in meshes:
        node_ranges.append((node_offset, node_offset + len(item.coordinates)))
        element_ranges.append((element_offset, element_offset + len(item.connectivity)))
        node_offset += len(item.coordinates)
        element_offset += len(item.connectivity)
    return MultiZoneModel(
        combine_meshes(meshes),
        zone_boundary_assignments(deck),
        zone_materials(deck),
        boundary_node_types(deck.sections, meshes, simulation_boundaries=() if deck.simulation is None else (deck.simulation,)),
        tuple(node_ranges), tuple(element_ranges),
    )


def _displacement_flags(deck: InputDeck, node_count: int) -> NDArray[np.bool_]:
    """Translate ``PRESCRIBED.CONDITIONS`` into the Fortran ``MBCU`` flags."""

    flags = np.zeros(2 * node_count, dtype=bool)
    for node, axis, condition in deck.prescribed_conditions:
        if node > node_count:
            raise ValueError(f"prescribed condition refers to node {node}, but mesh has {node_count} nodes")
        flags[2 * (node - 1) + (axis - 1)] = bool(condition)
    return flags


def _is_loaded_half_space(deck: InputDeck) -> bool:
    """Identify the open half-space surface used by the thesis ramp case.

    PREBEM's free-surface code cannot be used here because it also clears the
    applied surface traction.  The thesis model is an open loaded ground
    surface: it retains the prescribed ramp and uses the exterior half jump.
    """

    return (
        isinstance(deck.load, SurfaceRamp)
        and len(deck.sections) == 1
        and deck.sections[0].region in (3, 4)
    )


def run_laplace_case(deck_or_path: InputDeck | str | Path, *, progress=None) -> SolverResult:
    deck = deck_or_path if isinstance(deck_or_path, InputDeck) else read_input_deck(deck_or_path)
    if deck.domain not in (1, 4):
        raise ValueError("run_laplace_case requires DOMAIN 1 or 4")
    if deck.zones:
        from .multizone_workflow import run_multizone_laplace_case
        return run_multizone_laplace_case(deck, progress=progress)
    if progress is not None:
        progress("LD: preparing geometry and loading")
    mesh, material, node_types = _model(deck); grid = laplace_grid(deck.total_time, deck.steps)
    direct_traction = None
    if deck.inside_pressure:
        incident = np.zeros((len(grid.parameters), len(mesh.coordinates), 2, 2), dtype=complex)
        for i, parameter in enumerate(grid.parameters):
            incident[i, :, 0, 0] = deck.inside_pressure / parameter
            incident[i, :, 1, 1] = deck.inside_pressure / parameter
    elif isinstance(deck.load, TractionNodes):
        direct_traction = np.zeros((len(mesh.coordinates), 2), dtype=float)
        for node, tx, ty in deck.load.values:
            if 1 <= node <= len(mesh.coordinates): direct_traction[node - 1] = (tx, ty)
        incident = np.zeros((len(grid.parameters), len(mesh.coordinates), 2, 2), dtype=complex)
    elif isinstance(deck.load, SurfaceRamp):
        x = mesh.coordinates[:, 0] - deck.load.center_x
        half = deck.load.plateau_length / 2.0
        value = np.where(np.abs(x) <= half, deck.load.peak,
                         np.where(np.abs(x) <= half + deck.load.ramp_length,
                                  deck.load.peak * (half + deck.load.ramp_length - np.abs(x)) / deck.load.ramp_length, 0.0))
        direct_traction = np.column_stack((np.zeros_like(value), value))
        incident = np.zeros((len(grid.parameters), len(mesh.coordinates), 2, 2), dtype=complex)
    elif isinstance(deck.load, PlaneWave):
        incident = plane_wave_stress_spectrum(mesh, grid.parameters, material,
            longitudinal_stress=deck.load.longitudinal_stress,
            shear_stress=deck.load.transverse_shear, angle_degrees=deck.load.angle_degrees,
            transverse_stress=deck.load.vertical_stress, frequency=deck.load.frequency,
            wavefront_point=deck.load.wavefront_point)
    elif isinstance(deck.load, ExplosionPoint):
        incident = explosion_stress_spectrum(mesh, grid.parameters, material, source=deck.load)
    else:
        raise NotImplementedError(f"unsupported load type {type(deck.load).__name__}")
    meshes = [mesh_section(section) for section in deck.sections]
    normals = boundary_material_normals(meshes)
    if deck.free_field_stresses == 1 and np.any(node_types == 2):
        if not isinstance(deck.load, PlaneWave):
            raise NotImplementedError("reflected explosion loading is not translated")
        incident += plane_reflected_stress(mesh, deck.sections, meshes, deck.load, material,
                                          parameters=grid.parameters)
    if direct_traction is None:
        incident_traction_spectrum = stress_tractions(incident, normals)
    else:
        incident_traction_spectrum = np.repeat(
            direct_traction.reshape(1, -1), len(grid.parameters), axis=0
        ).astype(complex)
        incident_traction_spectrum /= grid.parameters[:, None]
    displacement_spectrum = np.empty_like(incident_traction_spectrum)
    traction_spectrum = np.empty_like(incident_traction_spectrum)
    flags = _displacement_flags(deck, len(mesh.coordinates))
    half_space = _is_loaded_half_space(deck)
    boundary_jump = 0.5 if half_space else 1.0
    for index, parameter in enumerate(grid.parameters):
        if progress is not None:
            progress(f"LD spectrum: sample {index + 1}/{len(grid.parameters)}")
        h_matrix, g_matrix = assemble_laplace_domain(
            mesh, parameter, material, jump=boundary_jump, node_types=node_types
        )
        prescribed = incident_traction_spectrum[index].reshape(-1).copy()
        prescribed[flags] = 0.0
        displacement_vector, traction_vector = solve_mixed_boundary_system(
            h_matrix, g_matrix, prescribed, flags,
            shear_modulus=material.shear_modulus,
        )
        displacement_spectrum[index] = displacement_vector
        traction_spectrum[index] = traction_vector
    displacement = inverse_laplace_spectrum(displacement_spectrum, deck.total_time).reshape(
        len(grid.time), len(mesh.coordinates), 2
    )
    traction = inverse_laplace_spectrum(traction_spectrum, deck.total_time).reshape(
        len(grid.time), len(mesh.coordinates), 2
    )
    # The inverse Laplace sum is periodic and Lanczos filtered. Applying it to
    # already-prescribed step tractions can introduce ringing at the end of the
    # record (for example, a constant -6.9 MPa ramp can become +3.6 MPa at the
    # last sample). Keep the exact input history on traction-prescribed DOFs;
    # only invert tractions that were solved from prescribed displacements.
    known_traction = None
    if direct_traction is not None:
        known_traction = direct_traction.reshape(-1)
    elif deck.inside_pressure:
        known_traction = (deck.inside_pressure * normals).reshape(-1)
    if known_traction is not None:
        traction_flat = traction.reshape(len(grid.time), -1)
        prescribed_traction = ~flags
        traction_flat[:, prescribed_traction] = known_traction[prescribed_traction]
    stress_spectrum = boundary_stress(mesh, material,
        displacement_spectrum.reshape(len(grid.parameters), -1),
        traction_spectrum.reshape(len(grid.parameters), -1),
        incident_local_stress(mesh, incident, normals=normals), grid.parameters)
    element_stress = inverse_laplace_spectrum(stress_spectrum, deck.total_time)
    if direct_traction is not None and not material.viscosity:
        # The known step traction was restored exactly above. Recover the
        # elastic time stresses from that same physical traction history;
        # inverting the stress spectrum separately reintroduces end-of-window
        # ringing and can even reverse the prescribed normal stress.
        element_stress = boundary_stress(
            mesh, material, displacement.reshape(len(grid.time), -1),
            traction.reshape(len(grid.time), -1),
            np.zeros((len(grid.time), len(mesh.coordinates), 3)),
        )
    if deck.inside_pressure:
        # PREBEM reports cavity internal-pressure hoop tension with the
        # opposite material-normal sign used by the free-field stress path.
        element_stress = -element_stress
    internal_coordinates = np.asarray(deck.internal_points, dtype=float).reshape(-1, 2)
    internal_spectrum = None
    internal_displacement = None
    if len(internal_coordinates):
        internal_spectrum = np.empty((len(grid.parameters), 2 * len(internal_coordinates)), dtype=np.complex128)
        for index, parameter in enumerate(grid.parameters):
            h_internal, g_internal = assemble_laplace_internal(mesh, internal_coordinates, parameter, material)
            internal_spectrum[index] = (
                g_internal @ (traction_spectrum[index].reshape(-1) / material.shear_modulus)
                - h_internal @ displacement_spectrum[index].reshape(-1)
            )
        internal_displacement = inverse_laplace_spectrum(internal_spectrum, deck.total_time)
    return SolverResult(
        grid.time, mesh.coordinates, mesh.connectivity, displacement, traction,
        displacement_spectrum, traction_spectrum, nodal_stress(mesh, element_stress),
        element_stress, stress_spectrum, internal_coordinates, internal_displacement,
        internal_spectrum,
    )


def run_time_case(deck_or_path: InputDeck | str | Path, *, progress=None) -> SolverResult:
    deck = deck_or_path if isinstance(deck_or_path, InputDeck) else read_input_deck(deck_or_path)
    if deck.domain not in (2, 5):
        raise ValueError("run_time_case requires DOMAIN 2 (direct TD) or DOMAIN 5 (CQ)")
    if deck.zones:
        from .multizone_workflow import run_multizone_time_case
        return run_multizone_time_case(deck, progress=progress)
    if deck.domain == 5:
        raise ValueError("DOMAIN 5 convolution quadrature requires a ZONES deck")
    if progress is not None:
        progress("TD: preparing geometry and loading")
    mesh, material, node_types = _model(deck)
    direct_traction = None
    # ExplosionPoint.frequency is the charge-wave scaling parameter, not a
    # harmonic-load frequency. Only PlaneWave frequency triggers the Laplace
    # harmonic path; direct-time explosion decks are supported below.
    if material.viscosity or (isinstance(deck.load, PlaneWave) and deck.load.frequency):
        raise ValueError("viscoelastic and harmonic loads require the Laplace workflow")
    if np.any(np.isin(node_types, (3, 4))):
        raise NotImplementedError("direct-time simulation and enclosing equations are not translated")
    if deck.inside_pressure:
        # Internal pressure is a prescribed radial boundary traction, not a
        # free-field stress wave. Keep the incident stress tensor zero and
        # construct the constant inward pressure history after normals exist.
        incident = np.zeros((deck.steps, len(mesh.coordinates), 2, 2), dtype=float)
    elif isinstance(deck.load, TractionNodes):
        direct_traction = np.zeros((len(mesh.coordinates), 2), dtype=float)
        for node, tx, ty in deck.load.values:
            if 1 <= node <= len(mesh.coordinates): direct_traction[node - 1] = (tx, ty)
        incident = np.zeros((deck.steps, len(mesh.coordinates), 2, 2), dtype=float)
    elif isinstance(deck.load, SurfaceRamp):
        x = mesh.coordinates[:, 0] - deck.load.center_x
        half = deck.load.plateau_length / 2.0
        value = np.where(np.abs(x) <= half, deck.load.peak,
                         np.where(np.abs(x) <= half + deck.load.ramp_length,
                                  deck.load.peak * (half + deck.load.ramp_length - np.abs(x)) / deck.load.ramp_length, 0.0))
        direct_traction = np.column_stack((np.zeros_like(value), value))
        incident = np.zeros((deck.steps, len(mesh.coordinates), 2, 2), dtype=float)
    elif isinstance(deck.load, PlaneWave):
        incident = plane_wave_stress_time(mesh, deck.total_time, deck.steps, material,
            longitudinal_stress=deck.load.longitudinal_stress,
            shear_stress=deck.load.transverse_shear, angle_degrees=deck.load.angle_degrees,
            transverse_stress=deck.load.vertical_stress,
            wavefront_point=deck.load.wavefront_point)
    elif isinstance(deck.load, ExplosionPoint):
        incident = explosion_stress_time(mesh, deck.total_time, deck.steps, material, source=deck.load)
    else:
        raise NotImplementedError(f"unsupported load type {type(deck.load).__name__}")
    meshes = [mesh_section(section) for section in deck.sections]
    normals = boundary_material_normals(meshes)
    if deck.free_field_stresses == 1 and np.any(node_types == 2):
        if not isinstance(deck.load, PlaneWave):
            raise NotImplementedError("reflected explosion loading is not translated")
        incident += plane_reflected_stress(mesh, deck.sections, meshes, deck.load, material,
                                          times=deck.total_time * np.arange(deck.steps) / deck.steps)
    if deck.inside_pressure:
        pressure_vector = (deck.inside_pressure * normals).reshape(-1)
        traction = np.tile(pressure_vector, (deck.steps, 1))
    elif direct_traction is not None:
        traction = np.tile(direct_traction.reshape(-1), (deck.steps, 1))
    else:
        traction = stress_tractions(incident, normals)
    # TIMB2 GENR explicitly clears actual free-surface tractions.
    traction[:, np.repeat(node_types == 2, 2)] = 0.0
    h_impulses = []; g_impulses = []
    half_space = _is_loaded_half_space(deck)
    boundary_jump = 0.5 if half_space else 1.0
    for step in range(1, deck.steps + 1):
        if progress is not None:
            progress(f"TD kernels: lag {step}/{deck.steps}")
        def kernel_progress(node, nodes, element, elements, stage):
            if progress is not None:
                progress(f"TD {stage}: lag {step}/{deck.steps}, "
                         f"node {node}/{nodes}, element {element}/{elements}")
        h_matrix, g_matrix = assemble_time_domain_step(
            mesh, step, deck.total_time, deck.steps, material,
            jump=boundary_jump, node_types=node_types,
            progress=kernel_progress if progress is not None else None,
        )
        h_impulses.append(h_matrix); g_impulses.append(g_matrix)
    if progress is not None:
        progress("TD: solving displacement history")
    flags = _displacement_flags(deck, len(mesh.coordinates))
    if flags.any():
        prescribed = traction.copy()
        prescribed[:, flags] = 0.0
        displacement, traction = solve_mixed_history(
            h_impulses, g_impulses, prescribed, flags, shear_modulus=material.shear_modulus,
            progress=progress,
        )
    else:
        displacement = solve_traction_history(
            h_impulses, g_impulses, traction, shear_modulus=material.shear_modulus,
            progress=progress,
        )
    if progress is not None:
        progress("TD: recovering boundary stresses")
    time = deck.total_time * np.arange(deck.steps) / deck.steps
    element_stress = boundary_stress(mesh, material, displacement, traction, incident_local_stress(mesh, incident, normals=normals))
    internal_coordinates = np.asarray(deck.internal_points, dtype=float).reshape(-1, 2)
    internal_displacement = None
    if len(internal_coordinates):
        internal_history = np.zeros((deck.steps, len(internal_coordinates), 2))
        boundary_nodes = np.full(len(internal_coordinates), -1, dtype=int)
        for point_index, point in enumerate(internal_coordinates):
            distances = np.linalg.norm(mesh.coordinates - point, axis=1)
            nearest = int(np.argmin(distances))
            if distances[nearest] <= 1.0e-9:
                boundary_nodes[point_index] = nearest
                internal_history[:, point_index] = displacement[:, 2 * nearest:2 * nearest + 2]
        interior_indices = np.flatnonzero(boundary_nodes < 0)
        if len(interior_indices):
            if progress is not None:
                progress(f"TD: internal points, lag 1/{deck.steps}")
            interior_coordinates = internal_coordinates[interior_indices]
            internal_values = np.zeros((deck.steps, 2 * len(interior_indices)))
            scaled_traction = traction / material.shear_modulus
            zero_boundary = np.zeros(displacement.shape[1])
            for lag in range(deck.steps):
                if progress is not None:
                    progress(f"TD: internal points, lag {lag + 1}/{deck.steps}")
                h_internal, g_internal = assemble_time_internal(
                    mesh, interior_coordinates, lag + 1,
                    deck.total_time, deck.steps, material,
                )
                for current in range(lag, deck.steps):
                    history = current - lag
                    previous = history - 1
                    average_traction = 0.5 * (
                        scaled_traction[history]
                        + (scaled_traction[previous] if previous >= 0 else zero_boundary)
                    )
                    average_displacement = 0.5 * (
                        displacement[history]
                        + (displacement[previous] if previous >= 0 else zero_boundary)
                    )
                    internal_values[current] += g_internal @ average_traction
                    internal_values[current] -= h_internal @ average_displacement
            internal_history[:, interior_indices] = internal_values.reshape(
                deck.steps, len(interior_indices), 2
            )
        internal_displacement = internal_history
    return SolverResult(time, mesh.coordinates, mesh.connectivity, displacement, traction,
        boundary_stress=nodal_stress(mesh, element_stress), element_stress=element_stress,
        internal_coordinates=internal_coordinates, internal_displacement=internal_displacement)


def save_result(result: SolverResult, path: str | Path) -> None:
    """Save a portable NumPy archive containing all available solution arrays."""

    values = {
        "time": result.time, "coordinates": result.coordinates,
        "connectivity": result.connectivity, "displacement": result.displacement,
        "traction": result.traction,
    }
    if result.displacement_spectrum is not None: values["displacement_spectrum"] = result.displacement_spectrum
    if result.traction_spectrum is not None: values["traction_spectrum"] = result.traction_spectrum
    if result.boundary_stress is not None: values["boundary_stress"] = result.boundary_stress
    if result.element_stress is not None: values["element_stress"] = result.element_stress
    if result.stress_spectrum is not None: values["stress_spectrum"] = result.stress_spectrum
    if result.internal_coordinates is not None: values["internal_coordinates"] = result.internal_coordinates
    if result.internal_displacement is not None: values["internal_displacement"] = result.internal_displacement
    if result.internal_displacement_spectrum is not None: values["internal_displacement_spectrum"] = result.internal_displacement_spectrum
    np.savez(path, **values)
