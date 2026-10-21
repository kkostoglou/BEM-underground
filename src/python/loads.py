"""Incident loading translated from PREBEM ``DATALAP/WINCID/WAVPLAN``."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .materials import Material, viscoelastic_properties
from .mesh import BoundaryMesh
from .input_deck import ExplosionPoint


def _closed_boundary_inward_material_normals(mesh: BoundaryMesh) -> NDArray[np.float64]:
    """Normals used for cavity loading, from previous/next boundary nodes."""

    return boundary_material_normals([mesh])


def boundary_material_normals(meshes: tuple[BoundaryMesh, ...] | list[BoundaryMesh]) -> NDArray[np.float64]:
    """Return PREBEM material-side normals for independent open/closed meshes.

    This is the normal construction used by ``REFLECT`` before converting a
    global free-field stress tensor to the boundary traction vector.
    """
    normals: list[NDArray[np.float64]] = []
    for mesh in meshes:
        coordinates = mesh.coordinates
        if len(coordinates) < 3:
            raise ValueError("a quadratic boundary needs at least three nodes")
        # Connectivity, not array order, defines neighbours. This also gives
        # each duplicate corner its own one-sided material normal.
        previous = np.full(len(coordinates), -1, dtype=int)
        following = np.full(len(coordinates), -1, dtype=int)
        for first, middle, last in mesh.connectivity:
            following[first] = middle
            previous[middle] = first
            following[middle] = last
            previous[last] = middle
        if np.any((previous < 0) & (following < 0)):
            raise ValueError("every boundary node must belong to an element")
        indices = np.arange(len(coordinates))
        tangent = coordinates[np.where(following < 0, indices, following)] - coordinates[np.where(previous < 0, indices, previous)]
        length = np.linalg.norm(tangent, axis=1)
        if np.any(length < 1.0e-11):
            raise ValueError("cannot determine a boundary normal at a duplicate node")
        normals.append(np.column_stack((-tangent[:, 1], tangent[:, 0])) / length[:, None])
    return np.vstack(normals)


def stress_tractions(stress: NDArray, normals: NDArray[np.float64]) -> NDArray:
    """Convert global stress histories to material-side traction histories."""
    values = np.asarray(stress)
    normal_array = np.asarray(normals, dtype=float)
    if values.ndim != 4 or values.shape[1:] != (len(normal_array), 2, 2):
        raise ValueError("stress must have shape (sample, node, 2, 2) matching normals")
    return np.einsum("snij,nj->sni", values, normal_array).reshape(len(values), -1)


def _legacy_wavefront_position(
    centered_coordinates: NDArray[np.float64], direction: NDArray[np.float64]
) -> float:
    """Reproduce WAVEFIR's flat-face roundoff advance to the adjacent node."""

    projections = centered_coordinates @ direction
    maximum = float(np.max(projections))
    tolerance = 1.0e-10 * max(1.0, float(np.max(np.abs(centered_coordinates))))
    support = centered_coordinates[np.abs(projections - maximum) <= tolerance]
    if len(support) >= 3:
        radii = np.sort(np.linalg.norm(support, axis=1))
        # On a flat face WAVEFIR's strict comparison sees the first tiny
        # trigonometric overshoot and advances from the midpoint to its neighbor.
        distinct = radii[np.r_[True, np.diff(radii) > tolerance]]
        if len(distinct) >= 2:
            return float(distinct[1])
    return maximum


def plane_wave_stress_spectrum(
    mesh: BoundaryMesh,
    laplace_parameters: ArrayLike,
    material: Material,
    *,
    longitudinal_stress: float,
    transverse_stress: float = 0.0,
    shear_stress: float = 0.0,
    angle_degrees: float = 0.0,
    frequency: float = 0.0,
    wavefront_point: tuple[float, float] | None = None,
) -> NDArray[np.complex128]:
    """WINCID incident global stress tensors, shaped (sample, node, 2, 2)."""

    parameters = np.asarray(laplace_parameters, dtype=np.complex128)
    if parameters.ndim != 1 or np.any(parameters == 0.0):
        raise ValueError("laplace_parameters must be a nonzero one-dimensional array")
    angle = np.deg2rad(angle_degrees)
    direction = np.array([np.cos(angle), np.sin(angle)])
    center = np.mean(mesh.coordinates, axis=0)
    centered = mesh.coordinates - center
    projections = centered @ direction
    wavefront = (_legacy_wavefront_position(centered, direction)
                 if wavefront_point is None else float((np.asarray(wavefront_point) - center) @ direction))
    distances = np.maximum(wavefront - projections, 0.0)
    result = np.zeros((len(parameters), len(mesh.coordinates), 2, 2), dtype=np.complex128)

    rotation = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    for step, parameter in enumerate(parameters):
        _, _, poisson, pressure_speed, shear_speed = viscoelastic_properties(parameter, material)
        additional_transverse = (
            longitudinal_stress * poisson / (1.0 - poisson)
            if longitudinal_stress < 0.0 else 0.0
        )
        additional_longitudinal = (
            transverse_stress * poisson / (1.0 - poisson)
            if transverse_stress < 0.0 else 0.0
        )
        p_stress = np.array(
            [
                [longitudinal_stress + additional_longitudinal, 0.0],
                [0.0, transverse_stress + additional_transverse],
            ],
            dtype=np.complex128,
        )
        sv_stress = np.array(
            [[0.0, shear_stress], [shear_stress, 0.0]], dtype=np.complex128
        )
        # DATALAP calls WINCID with the complete stress tensor for each active
        # wave family. Thus a combined P+SV deck intentionally applies that
        # tensor with both delays (a legacy cross-coupling preserved here).
        full_stress = p_stress + sv_stress
        global_stress = np.einsum("ik,kl,jl->ij", rotation, full_stress, rotation)
        harmonic = np.tanh(np.pi * parameter / (2 * frequency)) if frequency else 1.0
        p_delay = harmonic * np.exp(-parameter * distances / pressure_speed) / parameter
        sv_delay = harmonic * np.exp(-parameter * distances / shear_speed) / parameter
        if np.any(p_stress):
            result[step] += p_delay[:, None, None] * global_stress
        if np.any(sv_stress):
            result[step] += sv_delay[:, None, None] * global_stress
    return result


def plane_wave_stress_time(
    mesh: BoundaryMesh, total_time: float, step_count: int, material: Material, *,
    longitudinal_stress: float, transverse_stress: float = 0.0,
    shear_stress: float = 0.0, angle_degrees: float = 0.0,
    wavefront_point: tuple[float, float] | None = None,
) -> NDArray[np.float64]:
    """TWINCID incident global stress tensors, shaped (step, node, 2, 2)."""

    angle = np.deg2rad(angle_degrees)
    direction = np.array([np.cos(angle), np.sin(angle)])
    center = np.mean(mesh.coordinates, axis=0)
    centered = mesh.coordinates - center
    projections = centered @ direction
    wavefront = (_legacy_wavefront_position(centered, direction)
                 if wavefront_point is None else float((np.asarray(wavefront_point) - center) @ direction))
    distances = np.maximum(wavefront - projections, 0.0)
    poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))
    pressure_speed = np.sqrt(
        (material.lame_lambda + 2.0 * material.shear_modulus) / material.density
    )
    shear_speed = np.sqrt(material.shear_modulus / material.density)
    additional_transverse = longitudinal_stress * poisson / (1.0 - poisson) if longitudinal_stress < 0 else 0.0
    additional_longitudinal = transverse_stress * poisson / (1.0 - poisson) if transverse_stress < 0 else 0.0
    p_stress = np.array([
        [longitudinal_stress + additional_longitudinal, 0.0],
        [0.0, transverse_stress + additional_transverse],
    ])
    sv_stress = np.array([[0.0, shear_stress], [shear_stress, 0.0]])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    full_stress = np.einsum("ik,kl,jl->ij", rotation, p_stress + sv_stress, rotation)
    times = total_time * np.arange(step_count) / step_count
    p_arrived = times[:, None] >= (distances[None, :] / pressure_speed - 1.0e-9)
    sv_arrived = times[:, None] >= (distances[None, :] / shear_speed - 1.0e-9)
    stress = np.zeros((step_count, len(mesh.coordinates), 2, 2), dtype=np.float64)
    if np.any(p_stress):
        stress += p_arrived[:, :, None, None] * full_stress
    if np.any(sv_stress):
        stress += sv_arrived[:, :, None, None] * full_stress
    return stress


def _explosion_pressure(
    distance: NDArray[np.float64], speed: complex | float, source: ExplosionPoint,
    density: float,
) -> NDArray[np.complex128]:
    """ZEXPLOSQ/TEXPLOSQ amplitude before delay and coordinate rotation."""

    if np.any(distance < 1.0e-11):
        raise ValueError("the explosion point must not coincide with a boundary node")
    if source.formula == 1:
        scale = distance / source.wave_velocity ** (1.0 / 3.0)
        return source.coefficient * source.frequency * density * speed * scale ** (-source.exponent)
    if source.formula == 2:
        scale = distance / source.wave_velocity ** 0.5
        return source.coefficient * source.frequency * density * speed * scale ** (-source.exponent)
    return source.coefficient * source.frequency * source.wave_velocity ** source.exponent / distance ** source.secondary_exponent


def _explosion_global_stress(
    mesh: BoundaryMesh, axial: NDArray, shear: NDArray, poisson: complex | float,
    source: ExplosionPoint,
) -> NDArray:
    """EXPLMOR/TEXPLMOR radial-ray local stresses in global coordinates."""

    rays = mesh.coordinates - np.array([source.x, source.y])
    distance = np.linalg.norm(rays, axis=1)
    direction = rays / distance[:, None]
    cosine, sine = direction[:, 0], direction[:, 1]
    transverse = axial * poisson / (1.0 - poisson)
    sxx = axial * cosine**2 - 2.0 * shear * cosine * sine + transverse * sine**2
    syy = axial * sine**2 + 2.0 * shear * cosine * sine + transverse * cosine**2
    sxy = (axial - transverse) * cosine * sine + shear * (cosine**2 - sine**2)
    result = np.empty((len(distance), 2, 2), dtype=np.result_type(axial, shear, poisson))
    result[:, 0, 0] = sxx; result[:, 1, 1] = syy
    result[:, 0, 1] = result[:, 1, 0] = sxy
    return result


def explosion_stress_spectrum(
    mesh: BoundaryMesh, laplace_parameters: ArrayLike, material: Material, *, source: ExplosionPoint,
) -> NDArray[np.complex128]:
    """Full-space ``WAVEXPL`` incident stresses for a point explosion.

    This excludes the free-surface reflected rays handled by ``REFEXSV``.
    """

    parameters = np.asarray(laplace_parameters, dtype=np.complex128)
    if parameters.ndim != 1 or np.any(parameters == 0.0):
        raise ValueError("laplace_parameters must be a nonzero one-dimensional array")
    distance = np.linalg.norm(mesh.coordinates - np.array([source.x, source.y]), axis=1)
    result = np.zeros((len(parameters), len(distance), 2, 2), dtype=np.complex128)
    p_rate = source.p_amplitude
    for index, parameter in enumerate(parameters):
        _, _, poisson, pressure_speed, shear_speed = viscoelastic_properties(parameter, material)
        if p_rate:
            pressure = _explosion_pressure(distance, pressure_speed * p_rate, source, material.density)
            delay = np.exp(-parameter * distance / pressure_speed)
            if source.decay:
                axial = -pressure * (distance / pressure_speed) / (1.0 + parameter * distance / pressure_speed) * delay
            else:
                axial = -pressure * delay / parameter
            result[index] += _explosion_global_stress(mesh, axial, np.zeros_like(axial), poisson, source)
        if source.sv_amplitude:
            pressure = _explosion_pressure(distance, shear_speed * source.sv_amplitude, source, material.density)
            delay = np.exp(-parameter * distance / shear_speed)
            if source.decay:
                shear = -pressure * (distance / shear_speed) / (1.0 + parameter * distance / shear_speed) * delay
            else:
                shear = -pressure * delay / parameter
            result[index] += _explosion_global_stress(mesh, np.zeros_like(shear), shear, poisson, source)
    return result


def explosion_stress_time(
    mesh: BoundaryMesh, total_time: float, step_count: int, material: Material, *, source: ExplosionPoint,
) -> NDArray[np.float64]:
    """Full-space ``TWAVEXPL`` incident stresses for a point explosion."""

    if step_count <= 0 or total_time <= 0.0:
        raise ValueError("total_time and step_count must be positive")
    distance = np.linalg.norm(mesh.coordinates - np.array([source.x, source.y]), axis=1)
    if np.any(distance < 1.0e-11):
        raise ValueError("the explosion point must not coincide with a boundary node")
    poisson = material.lame_lambda / (2.0 * (material.lame_lambda + material.shear_modulus))
    pressure_speed = np.sqrt((material.lame_lambda + 2.0 * material.shear_modulus) / material.density)
    shear_speed = np.sqrt(material.shear_modulus / material.density)
    times = total_time * np.arange(step_count) / step_count
    result = np.zeros((step_count, len(distance), 2, 2), dtype=np.float64)
    def add_family(speed: float, amplitude: float, is_shear: bool) -> None:
        if amplitude == 0.0:
            return
        magnitude = np.real(_explosion_pressure(distance, speed * amplitude, source, material.density))
        arrival = distance / speed
        elapsed = times[:, None] - arrival[None, :]
        active = elapsed >= -1.0e-9
        values = -magnitude[None, :] * active
        if source.decay:
            values *= np.exp(-np.maximum(elapsed, 0.0) / arrival[None, :])
        for step in range(step_count):
            axial = np.zeros(len(distance)) if is_shear else values[step]
            shear = values[step] if is_shear else np.zeros(len(distance))
            result[step] += _explosion_global_stress(mesh, axial, shear, poisson, source)
    add_family(pressure_speed, source.p_amplitude, False)
    add_family(shear_speed, source.sv_amplitude, True)
    return result


def longitudinal_plane_wave_tractions(mesh, laplace_parameters, material, *,
        longitudinal_stress, transverse_stress=0.0, shear_stress=0.0, angle_degrees=0.0):
    """Opposite incident tractions applied to the scattered-field problem."""
    stress = plane_wave_stress_spectrum(mesh, laplace_parameters, material,
        longitudinal_stress=longitudinal_stress, transverse_stress=transverse_stress,
        shear_stress=shear_stress, angle_degrees=angle_degrees)
    normals = _closed_boundary_inward_material_normals(mesh)
    return np.einsum("snij,nj->sni", stress, normals).reshape(len(stress), -1)


def longitudinal_plane_wave_tractions_time(mesh, total_time, step_count, material, *,
        longitudinal_stress, transverse_stress=0.0, shear_stress=0.0, angle_degrees=0.0):
    """PRETIM's delayed Heaviside traction history."""
    stress = plane_wave_stress_time(mesh, total_time, step_count, material,
        longitudinal_stress=longitudinal_stress, transverse_stress=transverse_stress,
        shear_stress=shear_stress, angle_degrees=angle_degrees)
    normals = _closed_boundary_inward_material_normals(mesh)
    return np.einsum("snij,nj->sni", stress, normals).reshape(len(stress), -1)


def explosion_tractions(mesh, laplace_parameters, material, *, source: ExplosionPoint):
    stress = explosion_stress_spectrum(mesh, laplace_parameters, material, source=source)
    normals = _closed_boundary_inward_material_normals(mesh)
    return np.einsum("snij,nj->sni", stress, normals).reshape(len(stress), -1)


def explosion_tractions_time(mesh, total_time, step_count, material, *, source: ExplosionPoint):
    stress = explosion_stress_time(mesh, total_time, step_count, material, source=source)
    normals = _closed_boundary_inward_material_normals(mesh)
    return np.einsum("snij,nj->sni", stress, normals).reshape(len(stress), -1)
