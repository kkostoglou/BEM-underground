"""Time-marching BDF2 convolution quadrature for multizone elastodynamic BEM.

The coupled Laplace system is A(s) x(s) = B(s) p(s). Sampling A and B at
the BDF2 symbol and inverse-transforming the matrices gives causal boundary
influence weights A[j] and B[j]. At time n, the solver advances with

    A[0] x[n] = sum(j=0..n) B[j] p[n-j]
                  - sum(j=1..n) A[j] x[n-j].

The unknown vector x contains global displacements followed by interface
tractions. Thus the multizone compatibility and equilibrium equations are
enforced at every time step by the coupled PREGENR matrix assembly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .multizone_system import assemble_deck_multizone_system


@dataclass(frozen=True)
class CoupledTimeOperators:
    matrix: np.ndarray
    loading_matrix: np.ndarray
    node_count: int
    interface_nodes: tuple
    loading_nodes: tuple
    traction_scale: float


def _real_weights(spectrum, radius, count, label):
    weights = np.fft.ifft(spectrum, axis=0)[:count]
    weights *= radius ** -np.arange(count)[:, None, None]
    tolerance = 1e-7 * max(1.0, np.max(np.abs(weights.real)))
    if np.max(np.abs(weights.imag)) > tolerance:
        raise ValueError(f"{label} convolution weights lost real symmetry")
    return weights.real


def assemble_coupled_cq_weights(deck, *, progress=None):
    """Construct causal real time operators from coupled BEMZ matrices.

    Complex Laplace evaluations generate the weights only; no transformed
    load or response is solved here. The matrices retain physical node order.
    """
    if deck.steps < 2 or deck.total_time <= 0:
        raise ValueError("convolution quadrature requires at least two positive time steps")
    if not deck.zones:
        raise ValueError("convolution quadrature requires a ZONES deck")

    count = deck.steps
    samples = 1 << (2 * count - 1).bit_length()
    radius = 10.0 ** (-6.0 / samples)
    dt = deck.total_time / count
    matrix_spectrum = None
    loading_spectrum = None
    coordinates = None
    interface_nodes = None
    loading_nodes = None
    traction_scale = None
    node_count = None
    # Real operators have conjugate frequency pairs on the FFT contour.
    for sample in range(samples // 2 + 1):
        if progress is not None:
            progress(f"TD multizone weights: frequency {sample + 1}/{samples // 2 + 1}")
        z = radius * np.exp(-2j * np.pi * sample / samples)
        parameter = (1.5 - 2.0 * z + 0.5 * z * z) / dt
        xy, _, _, system = assemble_deck_multizone_system(deck, parameter)
        if sample == 0:
            coordinates = xy
            interface_nodes = system.interface_nodes
            loading_nodes = system.loading_nodes
            traction_scale = system.traction_scale
            node_count = system.node_count
            matrix_spectrum = np.empty((samples, *system.matrix.shape), dtype=complex)
            loading_spectrum = np.empty((samples, *system.loading_matrix.shape), dtype=complex)
        elif (system.interface_nodes != interface_nodes or
              system.loading_nodes != loading_nodes or
              system.node_count != node_count or
              system.traction_scale != traction_scale or
              not np.array_equal(xy, coordinates)):
            raise ValueError("multizone ordering changed between frequency samples")
        matrix_spectrum[sample] = system.matrix
        loading_spectrum[sample] = system.loading_matrix
        if 0 < sample < samples // 2:
            matrix_spectrum[-sample] = matrix_spectrum[sample].conj()
            loading_spectrum[-sample] = loading_spectrum[sample].conj()

    return CoupledTimeOperators(
        _real_weights(matrix_spectrum, radius, count, "system"),
        _real_weights(loading_spectrum, radius, count, "loading"),
        node_count, interface_nodes, loading_nodes, traction_scale)


def solve_coupled_convolution_quadrature(deck, prescribed_history, *, progress=None):
    """March coupled interface equations through a physical traction history.

    ``prescribed_history`` has shape ``(steps, loading nodes, 2)`` and uses
    physical stress units at times t = n * total_time / steps.
    """
    count = deck.steps
    known = np.asarray(prescribed_history, dtype=float)
    if known.ndim != 3 or known.shape[0] != count or known.shape[2] != 2:
        raise ValueError("prescribed history must have shape (steps, loading nodes, 2)")
    if not np.all(np.isfinite(known)):
        raise ValueError("prescribed history must be finite")
    operators = assemble_coupled_cq_weights(deck, progress=progress)
    if known.shape[1] != len(operators.loading_nodes):
        raise ValueError("prescribed history does not match loading nodes")

    # PREGENR scales traction columns by the largest zone shear modulus.
    scaled = known.reshape(count, -1) / operators.traction_scale
    solution = np.zeros((count, operators.matrix.shape[1]), dtype=float)
    for step in range(count):
        if progress is not None:
            progress(f"TD multizone march: step {step + 1}/{count}")
        rhs = np.zeros(solution.shape[1], dtype=float)
        for lag in range(step + 1):
            rhs += operators.loading_matrix[lag] @ scaled[step - lag]
        for lag in range(1, step + 1):
            rhs -= operators.matrix[lag] @ solution[step - lag]
        solution[step] = np.linalg.solve(operators.matrix[0], rhs)

    split = 2 * operators.node_count
    displacement = solution[:, :split].reshape(count, operators.node_count, 2)
    traction = (solution[:, split:].reshape(count, len(operators.interface_nodes), 2)
                * operators.traction_scale)
    return displacement, traction
