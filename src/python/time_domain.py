"""Direct time-domain fundamental solution from BEMTIM ``KERNEL/KERNA``."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .materials import Material
from .quadrature import leggauss


def _integrated_3d_kernel(
    source: NDArray[np.float64],
    field: NDArray[np.float64],
    normal: NDArray[np.float64],
    caps: tuple[float, float, float],
    h: NDArray[np.float64],
    branch: int,
    *,
    pressure_speed: float,
    shear_speed: float,
    density: float,
    time_step: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    identity = np.eye(3)
    aa = np.zeros((3, 3)); ab = np.zeros((3, 3)); ac = np.zeros((3, 3))
    ba = np.zeros((3, 3, 3)); bb = np.zeros_like(ba); bc = np.zeros_like(ba)
    bd = np.zeros_like(ba); be = np.zeros_like(ba)
    points, weights = leggauss(8)
    c1, c2 = pressure_speed, shear_speed

    for point, weight in zip(points, weights, strict=True):
        for cap, group in zip(caps, range(3), strict=True):
            x3 = 0.5 * cap * (1.0 + point)
            delta = np.array([field[0] - source[0], field[1] - source[1], x3])
            radius = float(np.linalg.norm(delta)); radial = delta / radius
            multiplier = weight * cap
            for i in range(3):
                for j in range(3):
                    a2 = radial[i] * radial[j] / radius
                    a3 = identity[i, j] / radius
                    if group == 0:
                        ab[i, j] += a2 / c1**2 * multiplier
                    elif group == 1:
                        ac[i, j] += (-a2 + a3) / c2**2 * multiplier
                    else:
                        aa[i, j] += (3.0 * a2 - a3) * multiplier
                    for k in range(3):
                        b1 = radial[i] * radial[j] * radial[k] / radius**2
                        b2 = (
                            identity[i, j] * radial[k] + identity[i, k] * radial[j]
                            + identity[j, k] * radial[i]
                        ) / radius**2
                        if group == 0:
                            b3 = radial[k] * identity[i, j] * (1.0 - 2.0 * c2**2 / c1**2) / radius**2
                            b4 = 2.0 * b1 * radius / c2
                            bb[i, j, k] += (-2.0 * (6.0 * b1 - b2) * c2**2 / c1**2 - b3) * multiplier
                            bd[i, j, k] += (-b4 * c2**3 / c1**3 - b3 * radius / c1) * multiplier
                        elif group == 1:
                            b4 = 2.0 * b1 * radius / c2
                            b5 = (identity[i, k] * radial[j] + identity[j, k] * radial[i]) / radius**2
                            bc[i, j, k] += (2.0 * (6.0 * b1 - b2) - b5) * multiplier
                            be[i, j, k] += (b4 - b5 * radius / c2) * multiplier
                        else:
                            ba[i, j, k] -= 6.0 * c2**2 * (5.0 * b1 - b2) * multiplier

    scale = np.ones(8)
    scale[6:] = 1.0 / time_step
    if branch == 1:
        scale[6:] *= -1.0
    else:
        scale[:6] = 0.0
    cai, cb, cc = h[3] - h[2], h[0], h[1]
    displacement = (aa * scale[0] * cai + ab * scale[1] * cb + ac * scale[2] * cc) / (4.0 * np.pi * density)
    stress = (
        ba * scale[3] * cai + bb * scale[4] * cb + bc * scale[5] * cc
        + bd * scale[6] * cb + be * scale[7] * cc
    ) / (4.0 * np.pi)
    traction = np.einsum("ijk,j->ik", stress, normal)
    return displacement[:2, :2], traction[:2, :2]


def time_domain_kernel(
    source: ArrayLike,
    field: ArrayLike,
    normal: ArrayLike,
    step: int,
    total_time: float,
    step_count: int,
    material: Material,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return the BEMTIM impulse kernel for one 1-based time ``step``."""

    y = np.asarray(source, dtype=float); x = np.asarray(field, dtype=float)
    normal2 = np.asarray(normal, dtype=float)
    if y.shape != (2,) or x.shape != (2,) or normal2.shape != (2,):
        raise ValueError("source, field, and normal must have shape (2,)")
    if not 1 <= step <= step_count or total_time <= 0.0:
        raise ValueError("invalid time step or total time")
    c1 = float(np.sqrt((material.lame_lambda + 2.0 * material.shear_modulus) / material.density))
    c2 = float(np.sqrt(material.shear_modulus / material.density))
    dt = total_time / step_count
    radius = float(np.linalg.norm(x - y))
    displacement = np.zeros((2, 2)); traction = np.zeros((2, 2))
    if step * dt < radius / c1:
        return displacement, traction
    normal3 = np.array([normal2[0], normal2[1], 0.0])

    def contribution(time: float, sign: float, branch: int) -> None:
        nonlocal displacement, traction
        if time <= 1.0e-15 or time < radius / c1:
            return
        h = np.zeros(4); h[0] = sign
        if time >= radius / c2: h[1] = sign
        h[2] = sign / (2.0 * c1**2)
        h[3] = sign / (2.0 * c2**2)
        if time < radius / c2: h[3] = sign * 0.5 * (time / radius) ** 2
        mixed_speed = c2 if time >= radius / c2 else radius / time
        cap1 = np.sqrt((time * c1) ** 2 - radius**2)
        cap2 = np.sqrt(max(0.0, (time * c2) ** 2 - radius**2))
        cap12 = np.sqrt(max(0.0, (time * mixed_speed) ** 2 - radius**2))
        u, t = _integrated_3d_kernel(
            y, x, normal3, (cap1, cap2, cap12), h, branch,
            pressure_speed=c1, shear_speed=c2, density=material.density, time_step=dt,
        )
        displacement += u; traction += t

    for branch in (1, 2):
        contribution((step + 1 - branch) * dt, 1.0, branch)
        contribution((step - branch) * dt, -1.0, branch)
    return displacement, traction
