"""Free-surface plane-wave reflection formulae from ``PRELAPL4.FOR``.

These routines provide the material- and angle-dependent coefficients used by
the legacy free-field loaders.  Ray construction and open-boundary assembly
remain separate concerns.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .materials import Material, viscoelastic_properties
from .stress import stress_global, stress_local


@dataclass(frozen=True)
class ReflectedPStress:
    """Local free-surface P/SV components returned by ``REFLD``."""

    p: tuple[complex, complex, complex]
    sv: tuple[complex, complex, complex]
    confined_p: tuple[complex, complex, complex]
    p_amplitude: complex
    sv_amplitude: complex


def free_surface_line(start: tuple[float, float], end: tuple[float, float]) -> tuple[tuple[float, float], float, tuple[float, float, float]]:
    """Translate ``SURFLINE`` for a straight reflecting free surface.

    Returns the outward normal, its global angle in ``[0, 2π)``, and the
    normalized line coefficients ``(CE, CF, CQ)`` for ``CE*x + CF*y - CQ``.
    """
    x1, y1 = start
    x2, y2 = end
    length = math.hypot(x2 - x1, y2 - y1)
    if length < 1.0e-11:
        raise ValueError("free-surface endpoints must differ")
    normal = ((y2 - y1) / length, (x1 - x2) / length)
    angle = math.atan2(normal[1], normal[0])
    if angle < 0.0:
        angle += 2.0 * math.pi
    return normal, angle, (normal[0], normal[1], (x1 * y2 - y1 * x2) / length)


def p_reflection_path_delays(
    points: ArrayLike, *, line: tuple[float, float, float], normal_angle: float,
    incident_angle: float, reflected_sv_angle: float, angle_case: int,
    wavefront: tuple[float, float, float], center: tuple[float, float],
    laplace_parameter: complex, material: Material,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.bool_], NDArray[np.bool_]]:
    """Translate the ray lengths and delays in ``PATP``.

    ``wavefront`` is the PREBEM line ``CA*x + CB*y - CP = 0`` at time zero.
    The returned P and SV delays include the complex wave speeds, followed by
    the two legacy behind-wavefront flags.
    """
    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("points must have shape (count, 2)")
    ce, cf, cq = line
    ca, cb, cp = wavefront
    normal = np.array((math.cos(normal_angle), math.sin(normal_angle)))
    distance_to_surface = np.abs(coordinates @ np.array((ce, cf)) - cq)
    theta_p = incident_angle if angle_case == 1 else 2.0 * math.pi - incident_angle
    theta_sv = reflected_sv_angle if angle_case == 1 else 2.0 * math.pi - reflected_sv_angle
    cosine_p, cosine_sv = math.cos(theta_p), math.cos(theta_sv)
    if abs(cosine_p) < 1.0e-11 or abs(cosine_sv) < 1.0e-11:
        raise ValueError("a reflected ray parallel to the surface has no finite path")
    p_leg = np.abs(distance_to_surface / cosine_p)
    sv_leg = np.abs(distance_to_surface / cosine_sv)
    tangential_p = p_leg * math.sin(theta_p)
    tangential_sv = sv_leg * math.sin(theta_sv)
    tangent = np.array((-normal[1], normal[0]))
    reflected_p_point = coordinates + distance_to_surface[:, None] * normal + tangential_p[:, None] * tangent
    reflected_sv_point = coordinates + distance_to_surface[:, None] * normal + tangential_sv[:, None] * tangent
    center_array = np.asarray(center, dtype=float)
    p_front = (reflected_p_point - center_array) @ np.array((ca, cb)) - cp
    sv_front = (reflected_sv_point - center_array) @ np.array((ca, cb)) - cp
    _lam, _mu, _nu, cp_speed, cs_speed = viscoelastic_properties(laplace_parameter, material)
    p_time = (np.abs(p_front) + p_leg) / cp_speed
    sv_time = np.abs(sv_front) / cp_speed + sv_leg / cs_speed
    parameter = complex(laplace_parameter)
    return np.exp(-parameter * p_time), np.exp(-parameter * sv_time), p_front > 0.0, sv_front > 0.0


def reflection_condition(wavefront_angle: float, surface_normal_angle: float) -> tuple[int, int | None]:
    """Translate ``REFLCOND``.

    Returns ``(kind, angle_case)`` where kind is 0 (none), 1 (reflection), or
    2 (wavefront parallel to the surface); ``angle_case`` is 1 or 2 when a
    reflection occurs.
    """
    pi = math.pi
    difference = wavefront_angle - surface_normal_angle
    if abs(abs(difference) - pi / 2.0) < 1.0e-6 or abs(abs(difference) - 3.0 * pi / 2.0) < 1.0e-6:
        return 2, None

    kind, angle_case = 0, None
    if 0.0 <= surface_normal_angle < pi / 2.0:
        if pi / 2.0 < difference < 3.0 * pi / 2.0:
            kind = 1
        if pi / 2.0 < difference < pi:
            angle_case = 1
        if pi < difference < 3.0 * pi / 2.0:
            angle_case = 2
    elif pi / 2.0 <= surface_normal_angle < pi:
        if pi / 2.0 < difference < 2.0 * pi - surface_normal_angle or -pi / 2.0 < difference < 0.0:
            kind = 1
        if pi / 2.0 < difference < pi:
            angle_case = 1
        if pi < difference < 2.0 * pi - surface_normal_angle or -pi / 2.0 < difference < 0.0:
            angle_case = 2
    elif pi <= surface_normal_angle < 3.0 * pi / 2.0:
        if pi / 2.0 < difference < 2.0 * pi - surface_normal_angle or -pi / 2.0 < difference < 0.0:
            kind = 1
        if pi / 2.0 < difference < 2.0 * pi - surface_normal_angle or -surface_normal_angle < difference < -pi:
            angle_case = 1
        if -pi < difference < 0.0:
            angle_case = 2
    elif 3.0 * pi / 2.0 <= surface_normal_angle < 2.0 * pi:
        if -3.0 * pi / 2.0 < difference < -pi / 2.0:
            kind = 1
        if -3.0 * pi / 2.0 < difference < -pi:
            angle_case = 1
        if -pi < difference < -pi / 2.0:
            angle_case = 2
    return kind, angle_case


def _properties(laplace_parameter: complex, material: Material) -> tuple[complex, complex, complex]:
    lam, mu, _nu, cp, cs = viscoelastic_properties(laplace_parameter, material)
    return lam, mu, cp / cs


def p_wave_amplitudes(
    incident_angle: float, reflected_p_angle: float, reflected_sv_angle: float,
    laplace_parameter: complex, material: Material,
) -> tuple[complex, complex]:
    """Return reflected P and SV amplitude ratios (``ZREFAMP``)."""
    _lam, _mu, ratio = _properties(laplace_parameter, material)
    sin_p = math.sin(2.0 * reflected_p_angle)
    sin_sv = math.sin(2.0 * reflected_sv_angle)
    cos_sv = math.cos(2.0 * reflected_sv_angle)
    denominator = sin_p * sin_sv + ratio * ratio * cos_sv**2
    if abs(denominator) < 1.0e-11:
        denominator = complex(1.0e-11)
    sine_incident = math.sin(2.0 * incident_angle)
    return ((sine_incident * sin_sv + ratio * ratio * cos_sv**2) / denominator,
            ratio * cos_sv * (sine_incident - sin_p) / denominator)


def p_reflected_stress_components(
    incident_angle: float, reflected_p_angle: float, reflected_sv_angle: float,
    incident_longitudinal: float, incident_transverse: float, surface_normal: tuple[float, float],
    laplace_parameter: complex, material: Material,
    *, global_incident_angle: float | None = None,
) -> ReflectedPStress:
    """Translate ``REFLD`` before its PATP ray delays are applied."""
    normal_x, normal_y = surface_normal
    _, _, poisson, _cp, _cs = viscoelastic_properties(laplace_parameter, material)
    p_ratio, sv_ratio = p_wave_amplitudes(
        incident_angle, reflected_p_angle, reflected_sv_angle, laplace_parameter, material
    )
    global_angle = incident_angle if global_incident_angle is None else global_incident_angle
    incident_global = stress_global(
        complex(incident_longitudinal), complex(incident_transverse), 0j,
        math.cos(global_angle), math.sin(global_angle),
    )
    s11, s22, s12 = stress_local(*incident_global, normal_x, normal_y)
    c11p, c22p, c12p, c11sv, c22sv, c12sv = p_wave_stress_coefficients(
        incident_angle, reflected_p_angle, reflected_sv_angle, laplace_parameter, material
    )
    p_stress = (c11p * s11, c22p * s22, c12p * s12)
    sv_stress = (c11sv * s11, c22sv * s22, c12sv * s12)
    reflected_axis = stress_local(
        *p_stress, math.cos(reflected_p_angle), math.sin(reflected_p_angle)
    )
    confined_axis = (
        reflected_axis[0] + reflected_axis[1] * poisson / (1.0 - poisson),
        reflected_axis[1] + reflected_axis[0] * poisson / (1.0 - poisson),
        reflected_axis[2],
    )
    confined_p = stress_global(
        *confined_axis, math.cos(reflected_p_angle), math.sin(reflected_p_angle)
    )
    return ReflectedPStress(p_stress, sv_stress, confined_p, p_ratio, sv_ratio)


def delayed_p_reflected_stress(
    points: ArrayLike, *, components: ReflectedPStress, line: tuple[float, float, float],
    normal_angle: float, incident_angle: float, reflected_sv_angle: float, angle_case: int,
    wavefront: tuple[float, float, float], center: tuple[float, float],
    laplace_parameter: complex, material: Material, use_confined_p: ArrayLike | bool = False,
) -> NDArray[np.complex128]:
    """Translate the delayed P+SV stress superposition at the end of ``REFP``.

    ``use_confined_p`` selects the confined reflected-P stress for non-free
    surface nodes, exactly as the Fortran branch on ``NREGI`` does.
    """
    coordinates = np.asarray(points, dtype=float)
    p_delay, sv_delay, p_behind, sv_behind = p_reflection_path_delays(
        coordinates, line=line, normal_angle=normal_angle, incident_angle=incident_angle,
        reflected_sv_angle=reflected_sv_angle, angle_case=angle_case,
        wavefront=wavefront, center=center, laplace_parameter=laplace_parameter, material=material,
    )
    confined = np.broadcast_to(np.asarray(use_confined_p, dtype=bool), (len(coordinates),))
    p_component = np.where(confined[:, None], np.asarray(components.confined_p), np.asarray(components.p))
    sv_component = np.broadcast_to(np.asarray(components.sv), (len(coordinates), 3)).copy()
    p_component[p_behind] = 0.0
    sv_component[sv_behind] = 0.0
    parameter = complex(laplace_parameter)
    local = (
        p_component * components.p_amplitude * p_delay[:, None] / parameter
        + sv_component * components.sv_amplitude * sv_delay[:, None] / parameter
    )
    cosine, sine = math.cos(normal_angle), math.sin(normal_angle)
    result = np.empty((len(coordinates), 2, 2), dtype=np.complex128)
    for index, (sxx, syy, sxy) in enumerate(local):
        gx, gy, gxy = stress_global(sxx, syy, sxy, cosine, sine)
        result[index] = ((gx, gxy), (gxy, gy))
    return result


def _sv_case(incident_angle: float, critical_angle: float, angle_case: int) -> int:
    # The Fortran comparison for KANGL=2 repeats the KANGL=1 predicates; this
    # intentionally preserves it, including its critical-angle tolerance.
    result = 1 if incident_angle > critical_angle else 3
    if incident_angle == critical_angle or abs(abs(critical_angle) - abs(incident_angle)) < 1.0e-3:
        result = 2
    return result


def sv_wave_amplitudes(
    incident_angle: float, reflected_p_angle: float, reflected_sv_angle: float,
    critical_angle: float, angle_case: int, laplace_parameter: complex, material: Material,
) -> tuple[complex, complex]:
    """Return reflected P and SV amplitude ratios (``ZREFAMSV``)."""
    _lam, _mu, ratio = _properties(laplace_parameter, material)
    case = _sv_case(incident_angle, critical_angle, angle_case)
    if case != 1:
        return 0j, 1 + 0j
    sin_p = math.sin(2.0 * reflected_p_angle)
    sin_sv = math.sin(2.0 * reflected_sv_angle)
    cos_sv = math.cos(2.0 * reflected_sv_angle)
    denominator = sin_p * sin_sv + ratio * ratio * cos_sv**2
    if abs(denominator) < 1.0e-11:
        denominator = complex(1.0e-11)
    return (ratio * math.sin(4.0 * reflected_sv_angle) / denominator,
            (sin_p * math.sin(2.0 * incident_angle) + ratio * ratio * cos_sv**2) / denominator)


def p_wave_stress_coefficients(
    incident_angle: float, _reflected_p_angle: float, reflected_sv_angle: float,
    laplace_parameter: complex, material: Material,
) -> tuple[complex, complex, complex, complex, complex, complex]:
    """Return P then SV (xx, yy, xy) stress coefficients from ``ZREFSTRP``."""
    lam, mu, ratio = _properties(laplace_parameter, material)
    cosine, sine = math.cos(incident_angle), math.sin(incident_angle)
    sin_sv, cos_sv = math.sin(2.0 * reflected_sv_angle), math.cos(2.0 * reflected_sv_angle)
    c11sv = ratio * mu * sin_sv / (lam + 2.0 * mu * cosine**2)
    c22sv = -ratio * mu * sin_sv / (lam + 2.0 * mu * sine**2)
    c12sv = 0j if abs(math.sin(2.0 * incident_angle)) < 1.0e-3 else -ratio * cos_sv / math.sin(2.0 * incident_angle)
    return -1 + 0j, -1 + 0j, 1 + 0j, c11sv, c22sv, c12sv


def sv_wave_stress_coefficients(
    incident_angle: float, reflected_p_angle: float, reflected_sv_angle: float,
    critical_angle: float, angle_case: int, laplace_parameter: complex, material: Material,
) -> tuple[complex, complex, complex, complex, complex, complex]:
    """Return P then SV (xx, yy, xy) stress coefficients from ``ZREFSTRSV``."""
    lam, mu, ratio = _properties(laplace_parameter, material)
    if _sv_case(incident_angle, critical_angle, angle_case) != 1:
        return 0j, 0j, 0j, -1 + 0j, -1 + 0j, -1 + 0j
    cosine_twice = math.cos(2.0 * incident_angle)
    if abs(cosine_twice) < 1.0e-11:
        return 0j, 0j, 0j, 1 + 0j, 1 + 0j, -1 + 0j
    sine_twice = math.sin(2.0 * incident_angle)
    return (
        ratio * math.cos(2.0 * reflected_sv_angle) / sine_twice,
        -(lam + 2.0 * mu * math.sin(reflected_p_angle)**2) / (ratio * mu * sine_twice),
        -math.sin(2.0 * reflected_p_angle) / (ratio * cosine_twice),
        1 + 0j, 1 + 0j, -1 + 0j,
    )
