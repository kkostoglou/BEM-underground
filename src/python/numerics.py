"""Foundational kernels translated from ``BEMLAPL/LAPLB3.FOR``."""

from __future__ import annotations

import cmath

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .materials import Material, viscoelastic_properties


def bessel_k012_legacy(z: complex) -> tuple[complex, complex, complex]:
    """Return the legacy approximations to modified Bessel K0, K1, and K2.

    This is a direct structured translation of Fortran ``BESSEL``. The
    original solver never permits a zero argument; it is rejected explicitly
    here instead of relying on platform-specific divide-by-zero behavior.
    """

    z = complex(z)
    if z == 0:
        raise ValueError("the Bessel kernel is singular at z=0")

    tolerance = 1.0e-9
    gamma = 0.5772156649
    scale = 1.253314137
    zs = 8.0 * z

    if abs(z) <= tolerance:
        k0 = -cmath.log(z)
        k1 = 1.0 / z
    elif abs(z) <= 3.5:
        i0 = (
            1.0
            + z**2 / 4.0
            + z**4 / 64.0
            + z**6 / 2304.0
            + z**8 / 147456.0
            + z**10 / 14745600.0
            + z**12 / 2123366400.0
        )
        k0 = -(cmath.log(z / 2.0) + gamma) * i0 + (
            z**2 / 4.0
            + 1.5 * z**4 / 64.0
            + 1.8333333333 * z**6 / 2304.0
            + 2.0833333333 * z**8 / 147456.0
            + 2.2833333333 * z**10 / 14745600.0
            + 2.45 * z**12 / 2123366400.0
        )
        i1 = (z / 2.0) * (
            1.0
            + z**2 / 8.0
            + z**4 / 192.0
            + z**6 / 9216.0
            + z**8 / 737280.0
            + z**10 / 88473600.0
            + z**12 / 14863564800.0
        )
        k1 = 1.0 / z + cmath.log(z / 2.0) * i1 - (z / 4.0) * (
            -0.1544313298
            + 1.345568670 * z**2 / 8.0
            + 2.178902003 * z**4 / 192.0
            + 2.762235336 * z**6 / 9216.0
            + 3.212235336 * z**8 / 737280.0
            + 3.578902003 * z**10 / 88473600.0
            + 3.888425813 * z**12 / 14863564800.0
        )
    elif abs(z) > 100.0 or (z.imag > 0.0 and abs(z.real) < tolerance):
        k0 = 0.0j
        k1 = 0.0j
    else:
        denominator = cmath.sqrt(z) * cmath.exp(z)
        k0 = scale * (
            1.0
            - 1.0 / zs
            + 9.0 / (2.0 * zs**2)
            - 225.0 / (6.0 * zs**3)
            + 11025.0 / (24.0 * zs**4)
            - 893025.0 / (120.0 * zs**5)
            + 108056025.0 / (720.0 * zs**6)
        ) / denominator
        k1 = scale * (
            1.0
            + 3.0 / zs
            - 15.0 / (2.0 * zs**2)
            + 315.0 / (6.0 * zs**3)
            - 14175.0 / (24.0 * zs**4)
            + 1091475.0 / (120.0 * zs**5)
            - 127702575.0 / (720.0 * zs**6)
        ) / denominator

    k2 = (2.0 / z) * k1 + k0
    return k0, k1, k2


def elastostatic_kernel(
    source: ArrayLike,
    field: ArrayLike,
    normal: ArrayLike,
    *,
    shear_modulus: float,
    poisson_ratio: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Translate Fortran ``ZERNEL`` displacement and traction kernels."""

    y = np.asarray(source, dtype=np.float64)
    x = np.asarray(field, dtype=np.float64)
    n = np.asarray(normal, dtype=np.float64)
    if y.shape != (2,) or x.shape != (2,) or n.shape != (2,):
        raise ValueError("source, field, and normal must each have shape (2,)")
    if shear_modulus <= 0.0:
        raise ValueError("shear_modulus must be positive")
    if poisson_ratio == 1.0:
        raise ValueError("poisson_ratio=1 makes the kernel singular")

    delta = x - y
    radius = float(np.hypot(delta[0], delta[1]))
    if radius == 0.0:
        raise ValueError("source and field points must be distinct")
    radial = delta / radius
    radial_normal = float(radial @ n)
    identity = np.eye(2)

    displacement = (
        -(3.0 - 4.0 * poisson_ratio) * np.log(radius) * identity
        + np.outer(radial, radial)
    ) / (8.0 * np.pi * shear_modulus * (1.0 - poisson_ratio))

    traction_kernel = np.empty((2, 2), dtype=np.float64)
    for i in range(2):
        for j in range(2):
            traction_kernel[i, j] = -(
                radial_normal
                * (
                    (1.0 - 2.0 * poisson_ratio) * identity[i, j]
                    + 2.0 * radial[i] * radial[j]
                )
                - (1.0 - 2.0 * poisson_ratio)
                * (radial[i] * n[j] - radial[j] * n[i])
            ) / (4.0 * np.pi * (1.0 - poisson_ratio) * radius)
    return displacement, traction_kernel


def laplace_domain_kernel(
    source: ArrayLike,
    field: ArrayLike,
    normal: ArrayLike,
    laplace_parameter: complex,
    material: Material,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Translate the complex elastodynamic Fortran ``KERNEL2`` routine."""

    y = np.asarray(source, dtype=np.float64)
    x = np.asarray(field, dtype=np.float64)
    n = np.asarray(normal, dtype=np.float64)
    if y.shape != (2,) or x.shape != (2,) or n.shape != (2,):
        raise ValueError("source, field, and normal must each have shape (2,)")
    z = complex(laplace_parameter)
    if z == 0.0:
        raise ValueError("laplace_parameter must be nonzero")

    _, mu, _, pressure_speed, shear_speed = viscoelastic_properties(z, material)
    delta = x - y
    radius = float(np.hypot(delta[0], delta[1]))
    if radius == 0.0:
        raise ValueError("source and field points must be distinct")
    radial = delta / radius
    radial_normal = float(radial @ n)
    c1_argument = z * radius / pressure_speed
    c2_argument = z * radius / shear_speed
    k01, k11, k21 = bessel_k012_legacy(c1_argument)
    k02, k12, k22 = bessel_k012_legacy(c2_argument)

    speed_ratio = shear_speed / pressure_speed
    psi = k02 + (k12 - speed_ratio * k11) / c2_argument
    chi = k22 - speed_ratio**2 * k21
    dpsi = (
        -k02 / radius
        - (z / shear_speed + 2.0 / (c2_argument * radius)) * k12
        + speed_ratio**2 * k01 / radius
        + 2.0 * shear_speed**2 * k11 / (z * pressure_speed * radius**2)
    )
    dchi = (
        -2.0 * k02 / radius
        - (z / shear_speed + 4.0 / (c2_argument * radius)) * k12
        + 2.0 * speed_ratio**2 * k01 / radius
        + (
            shear_speed**2 * z / pressure_speed**3
            + 4.0 * shear_speed**2 / (z * pressure_speed * radius**2)
        )
        * k11
    )
    e1 = dpsi - chi / radius
    e2 = 2.0 * chi / radius
    e3 = 2.0 * dchi * radial_normal
    e4 = (pressure_speed**2 / shear_speed**2 - 2.0) * (
        dpsi - dchi - chi / radius
    )

    identity = np.eye(2, dtype=np.complex128)
    displacement = (psi * identity - chi * np.outer(radial, radial)) / (2.0 * np.pi * mu)
    traction_kernel = np.empty((2, 2), dtype=np.complex128)
    for i in range(2):
        for j in range(2):
            traction_kernel[i, j] = (
                e1 * (radial_normal * identity[i, j] + radial[j] * n[i])
                - e2 * (n[j] * radial[i] - 2.0 * radial[i] * radial[j] * radial_normal)
                - e3 * radial[i] * radial[j]
                + e4 * radial[i] * n[j]
            ) / (2.0 * np.pi)
    return displacement, traction_kernel
