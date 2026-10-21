"""Plane P/SV reflected fields from PRELAPL2 and PRETIM2.

Angles used by the reflection coefficients are local to the reflecting surface;
the rotation of the incident tensor uses its independent global angle.
"""
from __future__ import annotations

import math
import numpy as np

from .loads import _legacy_wavefront_position
from .materials import viscoelastic_properties
from .reflection import (
    ReflectedPStress, free_surface_line, reflection_condition,
    p_reflected_stress_components, sv_wave_amplitudes, sv_wave_stress_coefficients,
)
from .stress import stress_global, stress_local


def sv_reflected_stress_components(incident_angle, reflected_p_angle, reflected_sv_angle,
                                  critical_angle, angle_case, shear_stress,
                                  global_incident_angle, normal, parameter, material):
    """REFLDSV: local reflected tensors and confinement before ray delays."""
    ratios = sv_wave_amplitudes(incident_angle, reflected_p_angle, reflected_sv_angle,
                               critical_angle, angle_case, parameter, material)
    coefficients = sv_wave_stress_coefficients(
        incident_angle, reflected_p_angle, reflected_sv_angle,
        critical_angle, angle_case, parameter, material)
    incident = stress_global(0j, 0j, complex(shear_stress),
                             math.cos(global_incident_angle), math.sin(global_incident_angle))
    local = np.asarray(stress_local(*incident, *normal))
    p = np.asarray(coefficients[:3]) * local
    sv = np.asarray(coefficients[3:]) * local
    axis = stress_local(*p, math.cos(reflected_p_angle), math.sin(reflected_p_angle))
    poisson = viscoelastic_properties(parameter, material)[2]
    confined = stress_global(axis[0] + axis[1] * poisson / (1 - poisson),
                             axis[1] + axis[0] * poisson / (1 - poisson), axis[2],
                             math.cos(reflected_p_angle), math.sin(reflected_p_angle))
    return ReflectedPStress(tuple(p), tuple(sv), confined, *ratios)


def reflected_ray_times(points, line, normal_angle, p_angle, sv_angle, angle_case,
                        wavefront, center, parameter, material, incident_family):
    """PATP/PATSV travel times; the SV branch retains its cosine tolerance."""
    points = np.asarray(points, dtype=float)
    normal = np.array([math.cos(normal_angle), math.sin(normal_angle)])
    tangent = np.array([-normal[1], normal[0]])
    depth = np.abs(points @ np.asarray(line[:2]) - line[2])
    cp, cs = viscoelastic_properties(parameter, material)[3:]
    speed = cp if incident_family == "P" else cs
    times, behind = [], []
    for angle, reflected_speed in ((p_angle, cp), (sv_angle, cs)):
        theta = angle if angle_case == 1 else 2 * math.pi - angle
        cosine = math.cos(theta)
        if abs(cosine) < 1e-9:
            cosine = 1e-9
        leg = depth / cosine
        reflection_points = points + depth[:, None] * normal + (leg * math.sin(theta))[:, None] * tangent
        front = (reflection_points - center) @ np.asarray(wavefront[:2]) - wavefront[2]
        times.append(np.abs(front) / speed + np.abs(leg) / reflected_speed)
        # PATSV does not set the behind-front flags in the retained source.
        behind.append(front > 0 if incident_family == "P" else np.zeros(len(points), dtype=bool))
    return times, behind


def plane_reflected_stress(mesh, sections, meshes, load, material, *, parameters=None, times=None):
    """REFP/REFSV or REFPTIM/REFSVTIM global stress histories.

    Straight free-surface segments are treated as infinite reflecting lines,
    as in SURFLINE. Input supplies exactly one of parameters and times.
    """
    if (parameters is None) == (times is None):
        raise ValueError("supply either Laplace parameters or time samples")
    samples = np.asarray(parameters if parameters is not None else times)
    frequency_domain = parameters is not None
    result = np.zeros((len(samples), len(mesh.coordinates), 2, 2),
                      dtype=complex if frequency_domain else float)
    global_angle = math.radians(load.angle_degrees) % (2 * math.pi)
    direction = np.array([math.cos(global_angle), math.sin(global_angle)])
    center = np.mean(mesh.coordinates, axis=0)
    front = (*direction, _legacy_wavefront_position(mesh.coordinates - center, direction))
    confined = np.concatenate([np.full(len(m.coordinates), s.region != 4)
                              for s, m in zip(sections, meshes, strict=True)])
    ratio = math.sqrt((material.lame_lambda + 2 * material.shear_modulus) / material.shear_modulus)
    for section, surface in zip(sections, meshes, strict=True):
        if section.region != 4:
            continue
        start, finish = surface.coordinates[[0, -1]]
        normal, angle, line = free_surface_line(tuple(start), tuple(finish))
        if np.max(np.abs(surface.coordinates @ np.asarray(line[:2]) - line[2])) > 1e-8:
            raise ValueError("analytic plane reflection requires a straight free surface")
        kind, angle_case = reflection_condition(global_angle, angle)
        if kind != 1:
            continue
        # REFLCOND leaves KANGL unset at exact normal incidence. Either side
        # gives the same physical ray and tensor; choose side 1 deterministically.
        angle_case = angle_case or 1
        theta = abs(math.pi - abs(global_angle - angle))
        for family, amplitude in (("P", load.longitudinal_stress), ("SV", load.transverse_shear)):
            if amplitude == 0:
                continue
            p_angle = theta if family == "P" else math.asin(min(1., ratio * math.sin(theta)))
            sv_angle = math.asin(math.sin(theta) / ratio) if family == "P" else theta
            side = 1 if angle_case == 1 else -1
            incident = math.pi - side * theta
            reflected_p, reflected_sv = math.pi + side * p_angle, math.pi + side * sv_angle
            critical = math.pi - side * math.asin(1 / ratio)
            for index, sample in enumerate(samples):
                parameter = complex(sample) if frequency_domain else 1 + 0j
                if family == "P":
                    components = p_reflected_stress_components(
                        incident, reflected_p, reflected_sv, amplitude,
                        getattr(load, "vertical_stress", 0.), normal, parameter, material,
                        global_incident_angle=global_angle)
                else:
                    components = sv_reflected_stress_components(
                        incident, reflected_p, reflected_sv, critical, angle_case,
                        amplitude, global_angle, normal, parameter, material)
                delays, behind = reflected_ray_times(
                    mesh.coordinates, line, angle, p_angle, sv_angle, angle_case,
                    front, center, parameter, material, family)
                p = np.where(confined[:, None], components.confined_p, components.p)
                sv = np.broadcast_to(components.sv, p.shape)
                local = np.zeros_like(p)
                for tensor, scale, delay, masked in zip(
                    (p, sv), (components.p_amplitude, components.sv_amplitude), delays, behind, strict=True
                ):
                    weight = np.exp(-parameter * delay) / parameter if frequency_domain else (sample >= delay.real - 1e-9)
                    local += tensor * (scale * weight * ~masked)[:, None]
                for node, stress in enumerate(local):
                    xx, yy, xy = stress_global(*stress, *normal)
                    tensor = np.array([[xx, xy], [xy, yy]])
                    result[index, node] += tensor if frequency_domain else tensor.real
    return result
