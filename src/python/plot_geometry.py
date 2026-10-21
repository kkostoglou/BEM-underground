"""Geometry plots follow element connectivity, including split corner nodes."""
from __future__ import annotations

import numpy as np

from .input_deck import ExplosionPoint, PlaneWave
from .loads import boundary_material_normals
from .mesh import BoundaryMesh


def boundary_plot_coordinates(coordinates, connectivity):
    """Return separate element polylines; never connect unrelated node IDs."""
    elements = np.asarray(coordinates)[np.asarray(connectivity, dtype=int)]
    separated = np.full((len(elements), 4, 2), np.nan)
    separated[:, :3] = elements
    return separated.reshape(-1, 2)


def draw_boundary(axes, coordinates, connectivity, *args, **kwargs):
    points = boundary_plot_coordinates(coordinates, connectivity)
    return axes.plot(points[:, 0], points[:, 1], *args, **kwargs)


def draw_explosion(axes, load, *, fontsize=7):
    if not isinstance(load, ExplosionPoint):
        return
    axes.plot([load.x], [load.y], marker='*', color='#b3261e', markersize=8,
              linestyle='none', label='Explosion', zorder=6)
    # Keep the source label inside compact inset axes.  A source near the
    # right side needs right alignment so the text does not run past the frame.
    xlim = axes.get_xlim()
    right_side = load.x >= (xlim[0] + xlim[1]) / 2 if xlim[1] > xlim[0] else False
    axes.annotate(f'Explosion ({load.x:g}, {load.y:g}) m', (load.x, load.y),
                  xytext=(-6 if right_side else 6, 10), textcoords='offset points',
                  ha='right' if right_side else 'left',
                  fontsize=fontsize, color='#b3261e', annotation_clip=True,
                  clip_on=True)


def draw_plane_wave(axes, coordinates, load, *, fontsize=6):
    """Show the incident P/SV wave front and direction at t=0."""
    if not isinstance(load, PlaneWave):
        return
    xy = np.asarray(coordinates)
    center = xy.mean(axis=0)
    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), 1e-12)
    angle = np.deg2rad(float(load.angle_degrees))
    normal = np.array([np.cos(angle), np.sin(angle)])
    tangent = np.array([-normal[1], normal[0]])
    if load.wavefront_point is None:
        front_center = center + normal * (.55 * span)
    else:
        offset = (np.asarray(load.wavefront_point) - center) @ normal
        front_center = center + normal * offset
    half = .22 * span
    front = np.vstack((front_center - half*tangent, front_center + half*tangent))
    axes.plot(front[:, 0], front[:, 1], color='#b3261e', lw=1.0, zorder=4)
    axes.annotate('', xy=center + normal*.18*span, xytext=center + normal*.48*span,
                  arrowprops=dict(arrowstyle='-|>', color='#b3261e', lw=.8))
    wave = 'P' if abs(load.longitudinal_stress) > 0 else 'SV'
    if abs(load.longitudinal_stress) > 0 and abs(load.transverse_shear) > 0:
        wave = 'P + SV'
    axes.text(front_center[0], front_center[1] + .27*span, f'{wave} wave, t=0',
              color='#b3261e', fontsize=fontsize, ha='center', va='bottom')


def draw_stress_profile(axes, coordinates, connectivity, stresses):
    """Offset each side by its signed stress, preserving corner discontinuity."""
    xy = np.asarray(coordinates)
    values = np.asarray(stresses)
    span = max(float(np.max(np.ptp(xy, axis=0))), 1e-12)
    peak = max(float(np.max(np.abs(values))), 1e-12)
    scale = .18 * span / peak
    normals = boundary_material_normals([BoundaryMesh(xy, np.asarray(connectivity))])
    profile = xy - scale * values[:, None] * normals
    draw_boundary(axes, xy, connectivity, color='#555555', linewidth=1, label='Boundary')
    draw_boundary(axes, profile, connectivity, color='#9b3c2d', linewidth=1.3,
                  label=f'Stress offset: {1/scale:.3g} MPa/m')
    # Label extrema on each straight side, including both values at a corner.
    indices = {int(np.argmin(values)), int(np.argmax(values))}
    _, inverse, counts = np.unique(xy.round(10), axis=0, return_inverse=True, return_counts=True)
    indices.update(np.flatnonzero(counts[inverse] > 1))
    for index in sorted(indices):
        normal = normals[index]
        axes.annotate(f'{values[index]:.1f}', profile[index],
                      xytext=tuple(9*normal), textcoords='offset points',
                      fontsize=6, ha='left' if normal[0] >= 0 else 'right',
                      va='bottom' if normal[1] >= 0 else 'top')
    # Keep the geometry circular/square without changing fixed plot limits.
    # ``datalim`` conflicts with time-history axes and produces repeated
    # Matplotlib warnings; ``box`` preserves the caller's limits.
    axes.set_aspect('equal', adjustable='box')
    axes.set_xlabel('X (m)')
    axes.set_ylabel('Y (m)')
    axes.margins(.2)


def draw_geometry_key(axes, coordinates, connectivity, nodes, load=None, *, inside_pressure=0.0, fontsize=6):
    """Label selected nodes on their own element side, with spaced leaders."""
    xy = np.asarray(coordinates)
    draw_boundary(axes, xy, connectivity, color='#555555', linewidth=1)
    axes.scatter(xy[:, 0], xy[:, 1], s=4 / 3, color='#333333', zorder=3)
    low, high = xy.min(axis=0), xy.max(axis=0)
    center = (low + high) / 2
    span = max(float(np.max(high - low)), 1e-12)
    _, inverse, counts = np.unique(xy.round(10), axis=0, return_inverse=True, return_counts=True)
    groups = {side: [] for side in range(4)}
    for node in nodes:
        index = int(node) - 1
        if not 0 <= index < len(xy):
            continue
        # An endpoint duplicated at a corner belongs to just one side.
        attached = [element for element in connectivity if index in element]
        tangent = xy[attached[0][-1]] - xy[attached[0][0]]
        if abs(tangent[0]) >= abs(tangent[1]):
            side = 0 if xy[index, 1] >= center[1] else 1
        else:
            side = 2 if xy[index, 0] >= center[0] else 3
        groups[side].append(index)
    for side, indices in groups.items():
        axis = 0 if side < 2 else 1
        indices.sort(key=lambda i: xy[i, axis])
        positions = np.array([xy[i, axis] for i in indices])
        for i in range(1, len(positions)):
            positions[i] = max(positions[i], positions[i-1] + .10 * span)
        if len(positions):
            positions -= max(0., positions[-1] - high[axis] - .08 * span)
        for index, position in zip(indices, positions):
            label = xy[index].copy()
            label[axis] = position
            label[1-axis] = (high[1-axis] + .20*span if side in (0, 2)
                             else low[1-axis] - .20*span)
            # Place the number just outside the node, without leader lines.
            radial = xy[index] - center
            norm = np.linalg.norm(radial) or 1.0
            label = xy[index] + radial / norm * (.045 * span)
            if counts[inverse[index]] > 1:
                # Coincident corner nodes belong to different sides; place
                # their labels outside those sides instead of on one diagonal.
                outward = np.array(((0., 1.), (0., -1.), (1., 0.), (-1., 0.))[side])
                label = xy[index] + outward * (.06 * span)
            axes.annotate(str(index+1), xy[index], xytext=label, ha='center', va='center',
                          fontsize=fontsize * 0.81,
                          bbox=dict(facecolor='white', edgecolor='none', pad=.2), zorder=5)
    draw_explosion(axes, load, fontsize=fontsize)
    draw_plane_wave(axes, xy, load, fontsize=fontsize)
    if inside_pressure:
        center = xy.mean(axis=0)
        axes.plot([center[0]], [center[1]], marker="*", color="#b3261e",
                  markersize=8, linestyle="none", zorder=7)
        axes.annotate(f"Central explosion\n{inside_pressure:g} MPa radial pressure",
                      center, xytext=(8, 8), textcoords="offset points",
                      fontsize=fontsize, color="#b3261e", ha="left", va="bottom",
                      bbox=dict(facecolor="white", edgecolor="none", alpha=.75, pad=.5),
                      annotation_clip=True)
    # Keep the inset self-explanatory when it is used beside result curves.
    # These are component directions, not additional loads on the mesh.
    legend = (r"$\mathbf{u}_{\mathrm{x}}$ →   $\mathbf{u}_{\mathrm{y}}$ ↑", r"$\mathbf{t}_{\mathrm{x}}$ →   $\mathbf{t}_{\mathrm{y}}$ ↑", r"$\sigma_{11}$ radial", r"$\sigma_{22}$ tangential", r"$\sigma_{12}$ shear")
    # Keep the component-direction key clear of the loading marker and its
    # label.  Explosion sources are often close to the square/circle, so place
    # the key in the lower-left corner of the inset for those cases.
    key_x, key_y = ((0.06, 0.16) if isinstance(load, ExplosionPoint) else (0.83, 0.56))
    axes.text(key_x, key_y, "\n".join(legend), transform=axes.transAxes,
              ha="left", va="top", fontsize=max(4.5, fontsize - 0.4),
              color="#334155", alpha=0.88, clip_on=False,
              bbox=dict(facecolor="white", edgecolor="#94a3b8", alpha=0.42, pad=1.5))
    right = max(high[0], load.x if isinstance(load, ExplosionPoint) else high[0])
    left = min(low[0], load.x if isinstance(load, ExplosionPoint) else low[0])
    bottom = min(low[1], load.y if isinstance(load, ExplosionPoint) else low[1])
    top = max(high[1], load.y if isinstance(load, ExplosionPoint) else high[1])
    axes.set_xlim(left-.35*span, right+.35*span)
    axes.set_ylim(bottom-.35*span, top+.35*span)
    axes.set_aspect('equal', adjustable='box')
    axes.axis('off')


def draw_liner_materials(axes, deck, *, fontsize=7):
    """Identify each bounded material zone inside its cavity on the input drawing."""
    for zone in deck.zones:
        boundary = abs(zone.outer_boundary)
        if not boundary or not zone.inner_boundary or boundary > len(deck.sections):
            continue
        section = deck.sections[boundary - 1]
        if section.shape != 1 or len(section.geometry) < 4:
            continue
        x, y = section.geometry[2:4]
        material = zone.material
        label = (f"Liner zone {zone.number}\n"
                 f"E = {material.young:,.0f} MPa\n"
                 f"ρ = {material.density:g} Mg/m³")
        inner_boundary = abs(zone.inner_boundary)
        if inner_boundary <= len(deck.sections):
            inner = deck.sections[inner_boundary - 1]
            if (inner.shape == 1 and len(inner.geometry) >= 4
                    and np.allclose(section.geometry[2:4], inner.geometry[2:4])
                    and np.isclose(section.geometry[1], 1)
                    and np.isclose(inner.geometry[1], 1)):
                label += (f"\n$R_i$ = {inner.geometry[0]:g} m"
                          f"\n$R_o$ = {section.geometry[0]:g} m")
        axes.text(x, y, label, ha="center", va="center", fontsize=fontsize,
                  color="#17324d", zorder=10,
                  bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                            edgecolor="#b8c4ce", alpha=0.92))
