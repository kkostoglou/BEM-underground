"""Material models translated from Fortran ``VISCOEL``."""

from __future__ import annotations

import cmath
from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    lame_lambda: float
    shear_modulus: float
    density: float
    viscosity: float = 0.0
    viscosity_model: int = 1

    def __post_init__(self) -> None:
        if self.shear_modulus <= 0.0 or self.density <= 0.0:
            raise ValueError("shear_modulus and density must be positive")
        if self.viscosity_model not in (1, 2):
            raise ValueError("viscosity_model must be 1 (Kelvin) or 2 (Maxwell)")


def viscoelastic_properties(
    laplace_parameter: complex, material: Material
) -> tuple[complex, complex, complex, complex, complex]:
    """Return complex lambda, mu, nu, P-wave speed, and S-wave speed."""

    z = complex(laplace_parameter)
    lam = complex(material.lame_lambda)
    mu = complex(material.shear_modulus)
    if abs(material.viscosity) >= 1.0e-9:
        if material.viscosity_model == 1:
            mu = material.shear_modulus + material.viscosity * z
            lam = material.lame_lambda - (2.0 / 3.0) * material.viscosity * z
        else:
            rate = material.shear_modulus / material.viscosity
            mu = material.shear_modulus * z / (z + rate)
            lam = (
                material.lame_lambda * z
                + rate * (material.lame_lambda + (2.0 / 3.0) * material.shear_modulus)
            ) / (z + rate)
    pressure_speed = cmath.sqrt((lam + 2.0 * mu) / material.density)
    shear_speed = cmath.sqrt(mu / material.density)
    poisson = lam / (2.0 * (lam + mu))
    return lam, mu, poisson, pressure_speed, shear_speed
