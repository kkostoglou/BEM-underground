"""End-to-end Laplace and time-domain workflows for coupled ZONES decks."""
from dataclasses import dataclass
import numpy as np
from .input_deck import ExplosionPoint, PlaneWave, read_input_deck
from .loads import (boundary_material_normals, explosion_stress_spectrum,
                    explosion_stress_time, plane_wave_stress_spectrum,
                    plane_wave_stress_time, stress_tractions)
from .materials import Material
from .mesh import BoundaryMesh, combine_meshes, mesh_section
from .multizone_system import assemble_deck_multizone_system
from .transforms import inverse_laplace_spectrum, laplace_grid


@dataclass(frozen=True)
class MultiZoneResult:
    time: np.ndarray
    coordinates: np.ndarray
    connectivity: np.ndarray
    displacement: np.ndarray
    traction: np.ndarray
    displacement_spectrum: np.ndarray
    traction_spectrum: np.ndarray
    boundary_stress: np.ndarray | None = None
    element_stress: np.ndarray | None = None
    stress_spectrum: np.ndarray | None = None
    internal_coordinates: np.ndarray | None = None
    internal_displacement: np.ndarray | None = None
    internal_displacement_spectrum: np.ndarray | None = None


def _material(item):
    shear = item.young / (2 * (1 + item.poisson))
    lame = 2 * shear * item.poisson / (1 - 2 * item.poisson)
    return Material(lame, shear, item.density / 1000,
                    item.viscosity, item.viscosity_model)


def recover_boundary_stresses(deck, meshes, displacement, traction, stress, parameters=None):
    """Recover one material-side stress per boundary in input node order.

    Shared boundaries use the finite zone's outer side (the liner side in
    the thesis plots). Never average stresses across a material interface:
    tangential stress can jump there. The cavity also uses liner material.
    Parameters select viscoelastic Laplace recovery; omit for elastic TD.
    """
    from .hoop import boundary_stress, incident_local_stress
    from .multizone_preprocess import prepare_zone_connectivity

    _, _, types, _, materials = prepare_zone_connectivity(deck)
    values = []
    offset = 0
    samples = len(displacement)
    u = displacement.reshape(samples, -1, 2)
    t = traction.reshape(samples, -1, 2)
    normals = boundary_material_normals(meshes)
    for mesh in meshes:
        nodes = slice(offset, offset + len(mesh.coordinates))
        members = np.flatnonzero(types[:, offset])
        # Types 6/8 are the inner side in the input connectivity; type 5
        # is a finite zone's outer boundary, whose normal is reversed.
        outer = [z for z in members if types[z, offset] == 5]
        zone = outer[0] if outer else members[0]
        reverse = types[zone, offset] == 5
        local_mesh = (BoundaryMesh(mesh.coordinates, mesh.connectivity[:, ::-1])
                      if reverse else mesh)
        incident = incident_local_stress(mesh, stress[:, nodes], normals=normals[nodes])
        if len(members) > 1:
            # BEMZ HOOP suppresses the incident contribution on interfaces.
            incident = np.zeros_like(incident)
        # Match BEMZ4.FOR HOOP: KZON reverses the outer connectivity, but
        # ZTR still receives the stored ZLOAD unchanged. Its proposed type-5
        # traction sign reversal is explicitly commented out in the source.
        # This is the legacy reporting convention used by the thesis plots.
        local = boundary_stress(
            local_mesh, _material(materials[zone]), u[:, nodes].reshape(samples, -1),
            t[:, nodes].reshape(samples, -1),
            incident, parameters,
        )
        values.append(local[:, :, ::-1] if reverse else local)
        offset = nodes.stop
    return np.concatenate(values, axis=1)


def run_multizone_laplace_case(deck_or_path, *, progress=None):
    """Solve a coupled DOMAIN 1/4 ZONES deck for prescribed incident traction."""
    deck = deck_or_path if hasattr(deck_or_path, 'zones') else read_input_deck(deck_or_path)
    if deck.domain not in (1, 4) or not deck.zones:
        raise ValueError("a Laplace-domain deck containing ZONES is required")
    if progress is not None:
        progress("LD: preparing zones and loading")
    meshes = [mesh_section(section) for section in deck.sections]
    mesh = combine_meshes(meshes)
    normals = boundary_material_normals(meshes)
    exterior = next((zone.material for zone in deck.zones
                     if zone.outer_boundary == 0), deck.material)
    material = _material(exterior)
    grid = laplace_grid(deck.total_time, deck.steps)
    if isinstance(deck.load, PlaneWave):
        stress = plane_wave_stress_spectrum(mesh, grid.parameters, material,
            longitudinal_stress=deck.load.longitudinal_stress,
            transverse_stress=deck.load.vertical_stress,
            shear_stress=deck.load.transverse_shear,
            angle_degrees=deck.load.angle_degrees,
            frequency=deck.load.frequency,
            wavefront_point=deck.load.wavefront_point)
    elif isinstance(deck.load, ExplosionPoint):
        stress = explosion_stress_spectrum(mesh, grid.parameters, material,
                                            source=deck.load)
    else:
        raise NotImplementedError(f"unsupported multi-zone load {type(deck.load).__name__}")
    incident = stress_tractions(stress, normals).reshape(deck.steps, -1, 2)
    displacement = np.empty((deck.steps, len(mesh.coordinates), 2), dtype=complex)
    traction = np.zeros_like(displacement)
    for sample, parameter in enumerate(grid.parameters):
        if progress is not None:
            progress(f"LD spectrum: sample {sample + 1}/{len(grid.parameters)}")
        coordinates, _, _, system = assemble_deck_multizone_system(deck, parameter)
        if not np.allclose(coordinates, mesh.coordinates):
            raise ValueError("preprocessor coordinate ordering changed during assembly")
        prescribed = np.array([incident[sample, node - 1]
                               for _, node in system.loading_nodes])
        solved_u, solved_t = system.solve_tractions(prescribed)
        displacement[sample] = solved_u
        for value, node in zip(solved_t, system.interface_nodes, strict=True):
            traction[sample, node - 1] = value
        for value, (_, node) in zip(prescribed, system.loading_nodes, strict=True):
            traction[sample, node - 1] = value
    displacement_flat = displacement.reshape(deck.steps, -1)
    traction_flat = traction.reshape(deck.steps, -1)
    from .hoop import nodal_stress
    if progress is not None:
        progress("LD: recovering boundary stresses")
    stress_values = recover_boundary_stresses(
        deck, meshes, displacement_flat, traction_flat, stress, grid.parameters)
    element_values = inverse_laplace_spectrum(stress_values, deck.total_time)
    boundary_values = nodal_stress(mesh, element_values)
    return MultiZoneResult(grid.time, mesh.coordinates, mesh.connectivity,
        inverse_laplace_spectrum(displacement_flat, deck.total_time),
        inverse_laplace_spectrum(traction_flat, deck.total_time),
        displacement_flat, traction_flat, boundary_values, element_values,
        stress_values)


def run_multizone_time_case(deck_or_path, *, method=None, progress=None):
    """March a coupled ZONES deck through the load history.

    DOMAIN 2 selects the original direct BEMTIM-style time kernels. DOMAIN 5
    selects convolution quadrature, which transforms BDF2 samples of the
    multizone Laplace operators into time-domain matrix weights.
    """
    deck = deck_or_path if hasattr(deck_or_path, 'zones') else read_input_deck(deck_or_path)
    if deck.domain not in (2, 5) or not deck.zones:
        raise ValueError("a DOMAIN 2 or DOMAIN 5 deck containing ZONES is required")
    domain_method = "direct" if deck.domain == 2 else "convolution_quadrature"
    if method is None:
        method = domain_method
    elif method not in ("convolution_quadrature", "direct"):
        raise ValueError("method must be 'convolution_quadrature' or 'direct'")
    elif method != domain_method:
        raise ValueError(
            f"DOMAIN {deck.domain} selects {domain_method!r}; use DOMAIN "
            f"{5 if method == 'convolution_quadrature' else 2} for {method!r}")
    if progress is not None:
        progress(f"TD ({'direct kernels' if method == 'direct' else 'BDF2 CQ'}): preparing zones and loading")
    if any(zone.material.viscosity for zone in deck.zones):
        raise ValueError("time-domain multi-zone materials cannot be viscoelastic")
    meshes = [mesh_section(section) for section in deck.sections]
    mesh = combine_meshes(meshes); normals = boundary_material_normals(meshes)
    exterior = next((zone.material for zone in deck.zones
                     if zone.outer_boundary == 0), deck.material)
    material = _material(exterior)
    if isinstance(deck.load, PlaneWave):
        stress = plane_wave_stress_time(mesh, deck.total_time, deck.steps, material,
            longitudinal_stress=deck.load.longitudinal_stress,
            transverse_stress=deck.load.vertical_stress,
            shear_stress=deck.load.transverse_shear,
            angle_degrees=deck.load.angle_degrees,
            wavefront_point=deck.load.wavefront_point)
    elif isinstance(deck.load, ExplosionPoint):
        stress = explosion_stress_time(mesh, deck.total_time, deck.steps,
                                       material, source=deck.load)
    else:
        raise NotImplementedError(f"unsupported multi-zone load {type(deck.load).__name__}")
    incident = stress_tractions(stress, normals).reshape(deck.steps, -1, 2)
    if method == "direct":
        from .multizone_time import assemble_deck_multizone_time_systems, solve_coupled_history
        coordinates, _, _, systems = assemble_deck_multizone_time_systems(deck, progress=progress)
        loading_nodes = systems[0].loading_nodes
        interface_nodes = systems[0].interface_nodes
    else:
        # GENRGR1's loading order depends only on zone membership. Read it
        # from one Laplace system before transforming the time history.
        from .multizone_system import assemble_deck_multizone_system
        coordinates, _, _, first = assemble_deck_multizone_system(
            deck, 1.5 * deck.steps / deck.total_time)
        loading_nodes = first.loading_nodes
        interface_nodes = first.interface_nodes
    prescribed = np.array([[incident[step, node - 1]
                            for _, node in loading_nodes]
                           for step in range(deck.steps)])
    if method == "direct":
        displacement, interface_traction = solve_coupled_history(
            systems, prescribed, progress=progress)
        interface_traction = interface_traction.reshape(
            deck.steps, len(interface_nodes), 2)
    else:
        from .multizone_cq import solve_coupled_convolution_quadrature
        displacement, interface_traction = solve_coupled_convolution_quadrature(
            deck, prescribed, progress=progress)
        displacement = displacement.reshape(deck.steps, -1)
    traction = np.zeros((deck.steps, len(coordinates), 2))
    for index, node in enumerate(interface_nodes):
        traction[:, node - 1, :] = interface_traction[:, index, :]
    for index, (_, node) in enumerate(loading_nodes):
        traction[:, node - 1, :] = prescribed[:, index, :]
    time = deck.total_time * np.arange(deck.steps) / deck.steps
    from .hoop import nodal_stress
    if progress is not None:
        progress("TD: recovering boundary stresses")
    element_stress = recover_boundary_stresses(
        deck, meshes, displacement.reshape(deck.steps, -1),
        traction.reshape(deck.steps, -1), stress,
    )
    recovered_stress = nodal_stress(mesh, element_stress)
    return MultiZoneResult(time, coordinates, mesh.connectivity, displacement,
        traction.reshape(deck.steps, -1), None, None,
        recovered_stress, element_stress, None)
