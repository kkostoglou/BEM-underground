"""Boundary-system rearrangement translated from Fortran ``BOUND``."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def solve_mixed_boundary_system(
    h_matrix: ArrayLike,
    scaled_g_matrix: ArrayLike,
    prescribed_values: ArrayLike,
    displacement_prescribed: ArrayLike,
    *,
    shear_modulus: float,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Return physical displacement and traction vectors.

    ``displacement_prescribed`` is the Python equivalent of ``MBCU``: false/0
    denotes prescribed traction and true/1 denotes prescribed displacement.
    ``scaled_g_matrix`` must use the legacy ``G * shear_modulus`` scaling.
    """

    h = np.asarray(h_matrix, dtype=np.complex128)
    g = np.asarray(scaled_g_matrix, dtype=np.complex128)
    prescribed = np.asarray(prescribed_values, dtype=np.complex128)
    flags = np.asarray(displacement_prescribed, dtype=bool)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or g.shape != h.shape:
        raise ValueError("H and G must be square matrices of equal shape")
    if prescribed.shape != (h.shape[0],) or flags.shape != prescribed.shape:
        raise ValueError("prescribed values and flags must match the matrix size")
    if shear_modulus <= 0.0:
        raise ValueError("shear_modulus must be positive")

    system = h.copy()
    right_operator = g.copy()
    working_values = prescribed.copy()
    working_values[~flags] /= shear_modulus
    for column in np.flatnonzero(flags):
        old_h = system[:, column].copy()
        system[:, column] = -right_operator[:, column]
        right_operator[:, column] = -old_h

    unknowns = np.linalg.solve(system, right_operator @ working_values)
    displacement = unknowns.copy()
    traction = prescribed.copy()
    displacement[flags] = prescribed[flags]
    traction[flags] = unknowns[flags] * shear_modulus
    return displacement, traction


def solve_traction_history(
    h_impulses: ArrayLike, scaled_g_impulses: ArrayLike,
    traction_history: ArrayLike, *, shear_modulus: float, progress=None,
) -> NDArray[np.float64]:
    """Translate BEMTIM ``SOLSCE`` for an all-traction boundary."""

    h = np.asarray(h_impulses, dtype=float); g = np.asarray(scaled_g_impulses, dtype=float)
    traction = np.asarray(traction_history, dtype=float) / shear_modulus
    if h.shape != g.shape or h.ndim != 3 or traction.shape != (len(h), h.shape[1]):
        raise ValueError("incompatible impulse matrices or traction history")
    displacement = np.zeros_like(traction); zero = np.zeros(h.shape[1])
    for current in range(len(h)):
        if progress is not None:
            progress(f"TD history solve: step {current + 1}/{len(h)}")
        right = np.zeros(h.shape[1])
        for lag in range(current + 1):
            newer = traction[current - lag]
            older = traction[current - lag - 1] if current - lag > 0 else zero
            right += g[lag] @ ((newer + older) / 2.0)
        for lag in range(1, current + 1):
            newer = displacement[current - lag]
            older = displacement[current - lag - 1] if current - lag > 0 else zero
            right -= h[lag] @ ((newer + older) / 2.0)
        displacement[current] = np.linalg.solve(h[0], right)
    return displacement


def solve_mixed_history(
    h_impulses: ArrayLike, scaled_g_impulses: ArrayLike,
    prescribed_values: ArrayLike, displacement_prescribed: ArrayLike, *, shear_modulus: float,
    progress=None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Translate BEMTIM ``BOUND`` plus ``SOLSCE`` for mixed boundaries.

    A true entry in ``displacement_prescribed`` denotes the legacy zero
    prescribed displacement condition.  False entries take their prescribed
    physical traction from ``prescribed_values``.  The column interchange is
    performed for every impulse matrix before evaluating the trapezoidal
    time-convolution recurrence.
    """
    h = np.asarray(h_impulses, dtype=float)
    g = np.asarray(scaled_g_impulses, dtype=float)
    prescribed = np.asarray(prescribed_values, dtype=float)
    flags = np.asarray(displacement_prescribed, dtype=bool)
    if h.shape != g.shape or h.ndim != 3 or prescribed.shape != (len(h), h.shape[1]):
        raise ValueError("incompatible impulse matrices or prescribed history")
    if flags.shape != (h.shape[1],) or shear_modulus <= 0.0:
        raise ValueError("invalid mixed-boundary flags or shear modulus")

    # BOUND interchanges H/G columns for displacement-prescribed degrees of
    # freedom.  The resulting A multiplies the unknown history and B the
    # known history.
    a = h.copy(); b = g.copy()
    for column in np.flatnonzero(flags):
        prior_h = a[:, :, column].copy()
        a[:, :, column] = -b[:, :, column]
        b[:, :, column] = -prior_h
    known = prescribed.copy()
    known[:, ~flags] /= shear_modulus
    known[:, flags] = prescribed[:, flags]
    unknown = np.zeros_like(known)
    zero = np.zeros(h.shape[1])
    inverse_first = np.linalg.inv(a[0])
    for current in range(len(h)):
        if progress is not None:
            progress(f"TD history solve: step {current + 1}/{len(h)}")
        right = np.zeros(h.shape[1])
        for lag in range(current + 1):
            newer = known[current - lag]
            older = known[current - lag - 1] if current - lag > 0 else zero
            right += b[lag] @ ((newer + older) / 2.0)
        for lag in range(1, current + 1):
            newer = unknown[current - lag]
            older = unknown[current - lag - 1] if current - lag > 0 else zero
            right -= a[lag] @ ((newer + older) / 2.0)
        unknown[current] = inverse_first @ right
    displacement = unknown.copy()
    traction = prescribed.copy()
    displacement[:, flags] = prescribed[:, flags]
    traction[:, flags] = unknown[:, flags] * shear_modulus
    return displacement, traction
