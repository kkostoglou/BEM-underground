"""Laplace-spectrum sampling and inversion from ``TRIGO``/``FASTINV``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LaplaceGrid:
    time: NDArray[np.float64]
    parameters: NDArray[np.complex128]
    convergence_factors: NDArray[np.float64]
    lanczos_factors: NDArray[np.float64]


def laplace_grid(total_time: float, sample_count: int) -> LaplaceGrid:
    """Construct the exact sampling grid used by Fortran ``LAPIN``/``TRIGO``."""

    if total_time <= 0.0 or sample_count < 2:
        raise ValueError("total_time must be positive and sample_count at least two")
    indices = np.arange(sample_count, dtype=np.float64)
    time = total_time * indices / sample_count
    alpha = 6.0 / total_time
    omega = 2.0 * np.pi / total_time
    parameters = alpha + 1j * omega * indices
    convergence = 2.0 * np.exp(alpha * time) / total_time
    arguments = np.pi * indices / sample_count
    lanczos = np.ones(sample_count)
    lanczos[1:] = np.sin(arguments[1:]) / arguments[1:]
    return LaplaceGrid(time, parameters.astype(np.complex128), convergence, lanczos)


def inverse_laplace_spectrum(
    spectrum: ArrayLike, total_time: float, *, axis: int = 0
) -> NDArray[np.float64]:
    """Apply the legacy Lanczos-filtered, unnormalised inverse DFT.

    ``spectrum`` is ordered at ``alpha + i*2*pi*k/T`` for ``k=0..N-1``.
    The returned axis contains values at ``t=k*T/N``.
    """

    values = np.asarray(spectrum, dtype=np.complex128)
    if values.ndim == 0:
        raise ValueError("spectrum must have at least one dimension")
    count = values.shape[axis]
    grid = laplace_grid(total_time, count)
    shape = [1] * values.ndim
    shape[axis] = count
    filtered = values * grid.lanczos_factors.reshape(shape)
    filtered = filtered.copy()
    dc_slice = [slice(None)] * values.ndim
    dc_slice[axis] = 0
    filtered[tuple(dc_slice)] = (
        0.5 * filtered[tuple(dc_slice)].real + 1j * filtered[tuple(dc_slice)].imag
    )
    transformed = (np.fft.ifft(filtered, axis=axis) * count).real
    # FASTINV writes AX(1)=2*AX(1) after FFCT. Although apparently intended to
    # restore the pre-transform DC value, it doubles the first output sample.
    first = [slice(None)] * values.ndim
    first[axis] = 0
    transformed[tuple(first)] *= 2.0
    return transformed * grid.convergence_factors.reshape(shape)
