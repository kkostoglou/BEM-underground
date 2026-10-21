"""Small desktop interface for running cases and reviewing portable results."""

from __future__ import annotations

import csv
import html
import webbrowser
import textwrap
import io
import multiprocessing
import os
import queue
import re
from datetime import datetime
import ctypes
import ctypes.wintypes
import time
import traceback
import uuid
from dataclasses import asdict, replace
from pathlib import Path
from tkinter import END, BOTH, LEFT, RIGHT, X, Y, BooleanVar, Canvas, StringVar, Tk, Toplevel, filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from .input_deck import ExplosionPoint, InputDeck, PlaneWave, SurfaceRamp, TractionNodes, read_input_deck
from .geometry_editor import GeometryEditor
from .plot_geometry import draw_boundary, draw_explosion, draw_geometry_key, draw_liner_materials, draw_stress_profile
from .case_catalog import ExampleCase, case_details, discover_examples, example_directory, identify_input_case
from .thesis_figures import chapter3_figures, thesis_graphs_for_case
from .gui_run_worker import run_gui_case
from .mesh import mesh_section, combine_meshes


REQUIRED_RESULT_KEYS = {"time", "coordinates", "connectivity", "displacement", "traction"}
RESULT_VIEWS = ("Time history", "Boundary profile", "Deformed geometry")
RESULT_QUANTITIES = ("displacement", "traction", "hoop stress", "radial stress", "shear stress")
STRESS_COMPONENTS = {"radial stress": 0, "hoop stress": 1, "shear stress": 2}
THESIS_CATEGORIES = (
    "2. LAPLACE DOMAIN", "3. TIME DOMAIN", "5. MULTI-ZONAL BEM", "6. WAVE PROPAGATION",
)
SELECT_EXAMPLE = "Select Example"
SELECT_SUBCATEGORY = "Select Sub-Category"
CUSTOM_CATEGORY = "CUSTOM INPUT FILE"
CUSTOM_SUBCATEGORY = "Custom DAT file"
THESIS_SUBCATEGORIES = {
    "2. LAPLACE DOMAIN": (
        "2.3.1 Plane Wave", "2.3.2 Explosion at a Point", "2.3.3 Explosion Inside a Cylindrical Cavity",
    ),
    "3. TIME DOMAIN": (
        "3.3.1 Plane Wave", "3.3.2 Explosion at a Point", "3.3.3 Explosion Inside a Cylindrical Cavity",
        "3.3.4 Solution at Internal Points of Half-plane",
    ),
    "5. MULTI-ZONAL BEM": ("5.3.1 Infinite Domain Solutions", "5.3.2 Semi-Infinite Domain Solutions"),
    "6. WAVE PROPAGATION": ("6.8.1 Half-space Solutions",),
}
INPUT_FIELDS = (
    ("steps", "Steps"), ("total_time", "Total time (s)"),
    ("young", "Young's modulus E (MPa)"), ("poisson", "Poisson's ratio ?"),
    ("density", "Density (Mg/m³)"), ("p", "Compressional P (MPa)"),
    ("sv", "Shear SV (MPa)"), ("angle", "Incidence angle (°)"),
    ("radius", "Cavity radius (m)"), ("height", "Rectangle Y-side length (m)"),
    ("load_x", "Explosion X coordinate (m)"), ("load_y", "Explosion Y coordinate (m)"),
    ("wave_velocity", "Explosion wave velocity (m/s)"), ("formula", "Explosion formula"),
    ("coefficient", "Explosion coefficient"), ("decay", "Explosion decay (0/1)"),
    ("coupling", "Explosion coupling F"),
)


def validate_input_numbers(deck: InputDeck, values: dict[str, str]) -> tuple[dict[str, float], dict[str, str]]:
    """Return parsed values and all field errors using the solver's input constraints."""
    parsed: dict[str, float] = {}
    errors: dict[str, str] = {}
    for key, label in INPUT_FIELDS:
        # Callers pass only fields applicable to the current DAT deck.  A
        # present-but-empty field is invalid; an absent field is disabled in
        # the form and must not prevent the case from running.
        if key not in values:
            continue
        raw = values[key].strip()
        if not raw:
            errors[key] = f"{label}: a value is required."
            continue
        try:
            value = float(raw)
        except ValueError:
            errors[key] = f"{label}: enter a number (for example, 0.015 or 1.5e-2)."
            continue
        if not np.isfinite(value):
            errors[key] = f"{label}: enter a finite number."
            continue
        parsed[key] = value
    if "steps" in parsed and (parsed["steps"] < 2 or not parsed["steps"].is_integer()):
        errors["steps"] = "Steps must be a whole number of at least 2."
    for key in ("total_time", "young", "density", "radius"):
        if key in parsed and parsed[key] <= 0:
            errors[key] = f"{dict(INPUT_FIELDS)[key]} must be positive."
    if "poisson" in parsed and not -1 < parsed["poisson"] < 0.5:
        errors["poisson"] = "Poisson's ratio must be greater than -1 and less than 0.5."
    return parsed, errors


def edited_input_deck(deck: InputDeck, values: dict[str, str]) -> InputDeck:
    """Validate editable values and build a solver input without changing the DAT file."""
    parsed, errors = validate_input_numbers(deck, values)
    if errors:
        raise ValueError("\n".join(errors.values()))
    material = replace(deck.material, young=parsed["young"], poisson=parsed["poisson"], density=parsed["density"])
    load = deck.load
    if isinstance(load, PlaneWave):
        load = replace(load, longitudinal_stress=parsed["p"],
                       transverse_shear=parsed["sv"], angle_degrees=parsed["angle"])
    elif isinstance(load, ExplosionPoint):
        load = replace(load, x=parsed["load_x"], y=parsed["load_y"], wave_velocity=parsed["wave_velocity"], p_amplitude=parsed["p"], sv_amplitude=parsed["sv"],
                       formula=int(parsed["formula"]), coefficient=parsed["coefficient"], decay=int(parsed["decay"]),
                       frequency=parsed["coupling"])
    sections = deck.sections
    if "radius" in parsed:
        geometry = list(sections[0].geometry); geometry[0] = parsed["radius"]
        if "height" in parsed: geometry[1] = parsed["height"]
        sections = (replace(sections[0], geometry=tuple(geometry)),)
    return replace(deck, steps=int(parsed["steps"]), total_time=parsed["total_time"],
                   material=material, load=load, sections=sections)


def load_result_archive(path: str | Path) -> dict[str, np.ndarray]:
    """Load and validate a solver result without leaving an open zip handle."""

    with np.load(path) as archive:
        missing = REQUIRED_RESULT_KEYS.difference(archive.files)
        if missing:
            raise ValueError(f"result is missing arrays: {', '.join(sorted(missing))}")
        result = {name: np.array(archive[name], copy=True) for name in archive.files}
    time = result["time"]
    coordinates = result["coordinates"]
    for name in ("displacement", "traction"):
        values = result[name]
        nodal_shape = (len(time), len(coordinates), 2)
        flat_shape = (len(time), 2 * len(coordinates))
        if values.shape == nodal_shape:
            # Laplace results retain the explicit (node, component) axes.
            # The graph and CSV code use the common interleaved layout
            # [x1, y1, x2, y2, ...], so normalize the archive on loading.
            result[name] = values.reshape(flat_shape)
        elif values.shape != flat_shape:
            raise ValueError(f"{name} has incompatible shape {values.shape}")
    if "boundary_stress" in result:
        stress = result["boundary_stress"]
        nodal_stress_shape = (len(time), len(coordinates), 3)
        flat_stress_shape = (len(time), 3 * len(coordinates))
        if stress.shape == flat_stress_shape:
            result["boundary_stress"] = stress.reshape(nodal_stress_shape)
        elif stress.shape != nodal_stress_shape:
            raise ValueError("boundary_stress must contain normal, hoop, and shear stress at every time and node")
    if "internal_displacement" in result and "internal_coordinates" in result:
        internal_coordinates = np.asarray(result["internal_coordinates"]).reshape(-1, 2)
        internal = np.asarray(result["internal_displacement"])
        internal_shape = (len(time), len(internal_coordinates), 2)
        flat_internal_shape = (len(time), 2 * len(internal_coordinates))
        if internal.shape == flat_internal_shape:
            result["internal_displacement"] = internal.reshape(internal_shape)
        elif internal.shape != internal_shape:
            raise ValueError(f"internal_displacement has incompatible shape {internal.shape}")
        result["internal_coordinates"] = internal_coordinates
    return result


def internal_displacement_histories(result: dict[str, np.ndarray], quantity: str) -> np.ndarray | None:
    """Return (time, internal point, component) histories for displacement plots."""

    if quantity != "displacement" or "internal_displacement" not in result:
        return None
    values = np.asarray(result["internal_displacement"])
    coordinates = np.asarray(result.get("internal_coordinates", ()))
    expected = (len(result["time"]), len(coordinates), 2)
    if values.shape == expected:
        return values
    if values.shape == (expected[0], 2 * expected[1]):
        return values.reshape(expected)
    raise ValueError(f"internal_displacement has incompatible shape {values.shape}")


def export_result_csv(result: dict[str, np.ndarray], path: str | Path) -> None:
    """Export time, displacement, and traction histories as readable columns."""

    node_count = len(result["coordinates"])
    coordinates = np.asarray(result["coordinates"])
    # Use two header rows: the first identifies each quantity and unit once;
    # the second identifies the node/component columns.  Blank cells preserve
    # the visual group span when opened in a spreadsheet.
    group_header = ["Time (s)"]
    node_header = ["Time"]
    for quantity in ("displacement", "traction"):
        unit = "m" if quantity == "displacement" else "MPa"
        label = "Displacement" if quantity == "displacement" else "Traction"
        group_header.extend([f"{label} ({unit})"] + [""] * (2 * node_count - 1))
        for node in range(1, node_count + 1):
            x, y = coordinates[node - 1]
            node_header.extend((f"Node {node} ({x:g}, {y:g}) X",
                                f"Node {node} ({x:g}, {y:g}) Y"))
    if "boundary_stress" in result:
        group_header.extend(["Boundary stress (MPa)"] + [""] * (3 * node_count - 1))
        for node in range(1, node_count + 1):
            x, y = coordinates[node - 1]
            location = f"Node {node} ({x:g}, {y:g})"
            node_header.extend(f"{name.title()} {location}" for name in ("radial", "hoop", "shear"))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(group_header)
        writer.writerow(node_header)
        for step, time_value in enumerate(result["time"]):
            writer.writerow(
                [float(time_value)]
                + result["displacement"][step].tolist()
                + result["traction"][step].tolist()
                + (result["boundary_stress"][step].ravel().tolist() if "boundary_stress" in result else [])
            )


def _deck_summary(deck: InputDeck, path: Path) -> str:
    if deck.domain == 1:
        domain = "Laplace domain"
    elif deck.domain == 2:
        domain = "Direct time domain (original TD kernels)"
    elif deck.domain == 5:
        domain = "Time domain (BDF2 convolution quadrature)"
    else:
        domain = str(deck.domain)
    load = deck.load
    if isinstance(load, PlaneWave):
        load_text = (
            f"Plane P/SV wave sweeping the boundary: P={load.longitudinal_stress:g} MPa, "
            f"shear={load.transverse_shear:g} MPa, angle={load.angle_degrees:g} degrees"
        )
    elif isinstance(load, SurfaceRamp):
        load_text = (f"Surface ramp: peak={load.peak:g} MPa, ramp length={load.ramp_length:g} m, "
                     f"plateau length={load.plateau_length:g} m, centre x={load.center_x:g} m")
    elif isinstance(load, TractionNodes):
        load_text = f"Direct nodal tractions: {len(load.values)} node assignments"
    else:
        load_text = (f"Explosive P/SV wave at ({load.x:g}, {load.y:g}) m; "
                     f"P multiplier={load.p_amplitude:g}, SV multiplier={load.sv_amplitude:g}, "
                     f"decay={load.decay:g}")
    reference_note = ""
    if path.stem.upper() == "IFEXRE44":
        reference_note = ("\nIFEXRE44 is the archived SV-only case. For PhD Figures 2.40-2.43, "
                          "select the rectangular explosion P-wave example THESIS_RECT_EXPLOSION.DAT.")
    shape_names = {1: "ellipse/circle", 2: "square", 3: "rectangle", 4: "straight line"}
    geometry_lines = []
    for index, section in enumerate(deck.sections, 1):
        name = shape_names.get(section.shape, f"combined shape ({-section.shape} segments)" if section.shape < 0 else str(section.shape))
        if section.shape == 2 and section.geometry:
            g = section.geometry
            grading = g[4] if len(g) > 4 else 1
            dis = g[5] if len(g) > 5 else 1
            grading_note = ("uniform elements" if grading == 1 else
                            f"symmetric grading toward both corners (ratio {abs(grading):g})" if grading < 0 else
                            f"progressive grading (ratio {grading:g})")
            dis_note = "separate corner nodes (discontinuous corner stresses)" if dis == 2 else "shared corner nodes"
            details = (f"side={g[0]:g} m, centre=({g[1]:g}, {g[2]:g}) m, rotation={g[3]:g} deg; "
                       f"grading={grading:g} ({grading_note}), DIS={dis:g} ({dis_note})")
        elif section.shape == 3 and section.geometry:
            g = section.geometry
            typx = g[5] if len(g) > 5 else 1.0
            abbreviated = len(g) == 7 and g[-1] in (1, 2)
            typy = 0.0 if abbreviated else (g[6] if len(g) > 6 else 1.0)
            dis = g[-1] if abbreviated else (g[7] if len(g) > 7 else 1.0)
            grading_text = f"X={typx:g}, Y={'automatic' if typy == 0 else f'{typy:g}'}"
            details = (f"X-side={g[0]:g} m, Y-side={g[1]:g} m, centre=({g[2]:g}, {g[3]:g}) m, "
                       f"rotation={g[4]:g} deg, grading ({grading_text}), DIS={dis:g}")
        elif section.shape == 1 and section.geometry:
            g = section.geometry
            details = f"semiaxis={g[0]:g} m, axis ratio={g[1] if len(g) > 1 else 1:g}, centre=({g[2] if len(g) > 2 else 0:g}, {g[3] if len(g) > 3 else 0:g}) m"
        elif section.shape == 4 and section.geometry:
            g = section.geometry
            details = f"from ({g[0]:g}, {g[1]:g}) m to ({g[2]:g}, {g[3]:g}) m, grading={g[4] if len(g) > 4 else 1:g}"
        else:
            details = "segment-defined"
        region_note = {1: "cavity", 2: "inclusion", 3: "open/free surface", 4: "open/interface", 5: "compound boundary"}.get(section.region)
        region_text = f"{section.region} ({region_note})" if region_note else str(section.region)
        element_note = "per side" if section.shape in (2, 3) else ""
        geometry_lines.append(f"Boundary {index}: {name}; region={region_text}; elements={section.x_elements} {element_note}"
                              + (f", y-elements={section.y_elements}" if section.y_elements else "")
                              + f"; {details}")
    return (
        f"Input: {path}\nTitle: {deck.title}\nAnalysis: {domain}\n"
        f"Boundaries: {len(deck.sections)}\nSteps: {deck.steps}\n"
        f"Total time: {deck.total_time:g}\nMaterial: E={deck.material.young:g} MPa, "
        f"Poisson={deck.material.poisson:g}, density={deck.material.density:g}\nGeometry:\n" +
        "\n".join(geometry_lines) + f"\nLoad: {load_text}{reference_note}"
    )


def _format_input_details(details: dict[str, object]) -> str:
    """Create a compact, readable view of every parsed DAT module."""
    lines: list[str] = []
    for key, value in details.items():
        if key == "title":
            continue
        text = str(value).replace("), ", "),\n    ") if key in {"sections", "internal_points", "prescribed_conditions"} else str(value)
        lines.append(f"{key}: {text}")
    return "\n".join(lines)


def _input_meaning(key: str, value: object) -> str:
    # Appendix F.6 module codes are shown using their full DAT path (for
    # example ``sections.0.shape`` or ``load.frequency``).
    leaf = key.rsplit(".", 1)[-1]
    # Accept both the canonical snake_case names and flattened DAT labels.
    normalized_key = key.replace(".", "_").lower()
    if key.startswith("plots["):
        plot_meanings = {
            1: "Boundary hoop stress: normal direction", 2: "Boundary hoop stress: tangential direction",
            3: "Boundary hoop stress: shear", 4: "Boundary hoop stress: third direction (plane strain)",
            5: "Deformed mesh shape", 6: "Plots at internal points", 11: "Cartesian hoop stress: normal",
            12: "Cartesian hoop stress: tangential", 13: "Cartesian hoop stress: shear",
            14: "Cartesian third-direction stress (plane strain)", 15: "Cartesian stress X (not implemented)",
            16: "Cartesian stress Y (not implemented)", 17: "Cartesian shear XY (not implemented)",
            18: "Third-direction boundary stress (plane strain)", 19: "Traction X (not implemented)",
            20: "Traction Y (not implemented)", 21: "Cartesian displacement X", 22: "Cartesian displacement Y",
            26: "Write boundary stress values only at NODE.PLOTS nodes", 27: "Element-by-element stress values (no averaging)",
            28: "Overlay NODE.PLOTS graphs on one page", 29: "NODE.PLOTS x-axis uses time-step numbers",
            30: "Preprocessing mesh plot",
        }
        return plot_meanings.get(int(value), "Plot selection code")
    if key.startswith("node_plots["):
        try:
            index = int(key.split("[", 1)[1].split("]", 1)[0]) + 1
            return (f"Selected boundary node #{index}; its displacement/stress history is used "
                    "for Cartesian plots (11-14, 21-22) and optional plot 26 output")
        except (ValueError, IndexError):
            return "Selected boundary node for Cartesian plots and optional plot 26 output"
    if leaf == "formula":
        return {1: "Cubic-root explosion formula (usually explosions)", 2: "Square-root explosion formula (usually blasting)", 3: "Other user-defined/SI explosion formula"}.get(int(value), "Explosion formula type")
    if leaf == "decay":
        return {0: "No stress-wave decay with time", 1: "Stress decays with time after wave arrival"}.get(int(value), "Explosion decay option")
    if leaf == "coefficient":
        return "Empirical explosion coefficient B (formula-dependent)"
    if leaf == "exponent":
        return "Explosion attenuation exponent N"
    if leaf == "secondary_exponent":
        return "Secondary exponent N1 (used by formula 3)"
    if leaf in {"p_amplitude", "sv_amplitude"}:
        return "Relative explosion-wave amplitude"
    if ".geometry[" in key:
        try:
            index = int(key.rsplit("[", 1)[1].split("]", 1)[0])
            shape = None
            # Section geometry ordering from Appendix F.6.3.
            if ".sections[" in key:
                shape_key = key.split(".geometry", 1)[0] + ".shape"
                # The caller does not pass the whole mapping, so retain the
                # generic positional explanation where shape is unavailable.
            meanings = {
                0: "Primary radius, side, or X-side length",
                1: "Axis ratio or Y-side length",
                2: "Centre X coordinate (or rotation/centre data)",
                3: "Centre Y coordinate (or rotation/centre data)",
                4: "Rotation angle in degrees",
                5: "Element grading/distribution type",
                6: "Second-side element grading/distribution type",
                7: "Corner discontinuity type",
            }
            return meanings.get(index, "Geometry parameter")
        except (ValueError, IndexError):
            pass
    maps = {
        "domain": {1: "Laplace transform domain (transient elastodynamic)", 2: "Time domain (transient elastodynamic)", 3: "Elastostatic", 4: "Laplace transform domain (steady-state / frequency)"},
        "analysis": {1: "Plane strain", 2: "Plane stress"},
        "units": {1: "SI International System: lengths m, time s, stress/modulus MPa, density Mg/m³", 2: "Imperial system: lengths ft/in, time s, stress/modulus in consistent imperial units", 3: "User-defined combination; all input quantities must be mutually consistent"},
        "free_field_stresses": {
            0: "Default method: apply the incident field to every boundary, including a discretized free surface; suitable when the model treats the free surface as another boundary of an infinite-region solution",
            1: "Harmonic-wave method: reflect the incident wave at the free surface using reflection laws, then calculate free field as incident plus reflected waves; use for steady harmonic loading with a free ground surface",
            2: "Subregion/equilibrium method: divide the ground into subregions and enforce equilibrium at the subsurface interface; use when the free field must be assembled from layered or separate regions",
        },
        "simulation": {
            1: "Simulation boundary: this boundary represents the free surface/interface used to simulate the continuation of the surrounding medium",
            2: "Enclosing boundary: this boundary surrounds an inclusion or internal zone and is treated as an internal simulation interface",
        },
    }
    if key == "units" and value is None:
        return "SI international system units (default when UNITS is omitted)"
    # Keep omitted optional modules visible with their complete option list.
    if normalized_key.endswith("free_field_stresses"):
        selected = maps["free_field_stresses"].get(int(value), "Unknown option") if value is not None else "0 (default)"
        return (f"Selected: {selected}. Options: 0 (default) incident field on all boundaries, including a discretized free surface; "
                "1 harmonic reflection at the free surface (incident + reflected field); "
                "2 subregions with equilibrium at subsurface interfaces. Omitted = 0.")
    if normalized_key.endswith("simulation"):
        selected = maps["simulation"].get(int(value), "No simulation/enclosing boundary") if value is not None else "omitted"
        return (f"Selected: {selected}. Options: 1 simulation/free-surface continuation boundary; "
                "2 enclosing boundary around an inclusion or internal zone. Omitted = no simulation boundary.")
    if leaf == "shape":
        return {1: "Ellipse or circle", 2: "Square", 3: "Rectangle", 4: "Straight line"}.get(int(value), "Combined shape / segment")
    if leaf == "region":
        return {1: "Closed external surface (cavity)", 2: "Closed internal surface (inclusion)", 3: "Open surface", 4: "Open free traction surface (reflecting)", 5: "Enclosing surface, anticlockwise", 6: "Enclosing surface, clockwise"}.get(int(value), "")
    if leaf == "viscosity_model" or leaf == "mode":
        return {1: "Kelvin solid", 2: "Maxwell fluid"}.get(int(value), "")
    if leaf == "young":
        return "Young's modulus E (MPa in SI output convention)"
    if leaf == "poisson":
        return "Poisson's ratio NU"
    if leaf == "density":
        return "Mass density RO (Mg/m³ in SI convention)"
    if leaf == "viscosity":
        return "Viscosity coefficient (same stress-time convention as E)"
    if leaf == "longitudinal_stress":
        return "Longitudinal P-wave stress/amplitude"
    if leaf == "transverse_shear":
        return "Transverse SV-wave stress/amplitude"
    if leaf == "angle_degrees":
        return "Incident wave angle (degrees)"
    if leaf in {"total_time", "time"}:
        return "Total or elapsed time (s)"
    if leaf == "wave_velocity":
        return "Explosion wave velocity (length/s in the selected unit system)"
    if leaf in {"x", "y"} and key.startswith("load."):
        return "Explosion source coordinate (selected length unit)"
    if leaf in {"x_elements", "y_elements"}:
        return "Number of boundary elements"
    if leaf == "geometry" or ".geometry[" in key:
        return "Shape dimensions, centre, rotation, and grading data"
    if leaf in {"internal_points", "prescribed_conditions", "element_nodes"}:
        return "Values read from the corresponding DAT module"
    try:
        meaning_map = maps.get(key) or maps.get(normalized_key) or maps.get(leaf)
        if meaning_map is not None:
            return meaning_map.get(int(value), f"{key} option {value}")
        return ""
    except (TypeError, ValueError):
        return ""


def _display_input_value(value: object) -> str:
    """Use compact, deterministic formatting for values shown in the DAT table."""
    if isinstance(value, float):
        if value != value:
            return "NaN"
        if value == float("inf") or value == float("-inf"):
            return str(value)
        return f"{value:.8g}"
    return str(value)


class BemApplication:
    def __init__(self, root: Tk) -> None:
        self.root = root
        self.root.title("Boundary Element Underground — Python Solver")
        self.root.geometry("1100x760")
        self.input_path = StringVar()
        self.output_path = StringVar()
        default_input_dir = example_directory()
        self.last_input_dir_file = Path(__file__).resolve().parents[2] / "examples" / "results" / ".last_input_folder"
        self.last_input_dir = default_input_dir
        if self.last_input_dir_file.is_file():
            try:
                remembered = Path(self.last_input_dir_file.read_text(encoding="utf-8").strip())
                if remembered.is_dir():
                    self.last_input_dir = remembered
            except OSError:
                pass
        self.quantity = StringVar(value="displacement")
        self.component = StringVar(value="X")
        self.node = StringVar(value="1")
        self.show_all = BooleanVar(value=False)
        self.view = StringVar(value="Time history")
        self.view_label = StringVar()
        self.quantity_label = StringVar()
        self.results_title = StringVar(value="No test selected")
        self.thesis_graph = StringVar()
        self.graph_zoom = 1.0
        self.graph_zoom_label = StringVar(value="100%")
        self.step = StringVar(value="1")
        self.case_filter = StringVar(value=SELECT_EXAMPLE)
        self.case_subcategory = StringVar(value=SELECT_SUBCATEGORY)
        self.case_choice = StringVar()
        self.examples = list(discover_examples())
        # Include every readable archived DAT in the selectors.  Recognized
        # thesis cases retain their precise application labels; other decks
        # use their TITLE as the selectable sub-category.
        known = {case.path.resolve() for case in self.examples}
        for path in sorted(example_directory().glob("*.DAT")):
            if path.resolve() in known:
                continue
            try:
                deck = read_input_deck(path)
            except (OSError, ValueError):
                continue
            family = "2. LAPLACE DOMAIN" if deck.domain == 1 else "3. TIME DOMAIN" if deck.domain == 2 else "5. MULTI-ZONAL BEM"
            subcategory = ("2.3.1 Plane Wave" if "PLANE" in deck.title.upper() else
                           "2.3.2 Explosion at a Point" if "EXP" in deck.title.upper() else
                           "Archived DAT cases")
            self.examples.append(ExampleCase(path, deck, family, path.stem.upper() + " — " + deck.title.strip(), "", subcategory))
        self.examples = tuple(self.examples)
        self.filtered_examples: tuple[ExampleCase, ...] = self.examples
        self.result: dict[str, np.ndarray] | None = None
        self.result_path: Path | None = None
        self.active_deck: InputDeck | None = None
        self.active_case: ExampleCase | None = None
        self.thesis_history_nodes: tuple[int, ...] | None = None
        self.active_thesis_section = ""
        self._running = False
        self._solver_process = None
        self._solver_queue = None
        self._solver_temporary_path: Path | None = None
        self._solver_error_path: Path | None = None
        self._solver_output_path: Path | None = None
        self._run_id: str | None = None
        self._solver_stage = ""
        self._progress_started = 0.0
        self.run_stage = StringVar(value="Ready")
        self.editor_deck: InputDeck | None = None
        self.editor_path: Path | None = None
        self.editor_case: ExampleCase | None = None
        self.input_values = {key: StringVar() for key, _ in INPUT_FIELDS}
        self.nodes_text = StringVar()
        self.input_reference = StringVar()
        self.input_boundaries = StringVar()
        self.input_description = StringVar()
        self.input_full_details = StringVar()
        self._loading_inputs = False
        self._inputs_valid = False
        self._modified_inputs = False
        self._plot_after_id: str | None = None
        self._build()
        # Start with a complete, editable geometry instead of disconnected
        # placeholder numbers that cannot be passed to the solver.
        default_case = next((case for case in self.examples if case.path.stem.upper() == "INFPL1"), None)
        if default_case is not None:
            self.case_filter.set(default_case.family)
            self.filter_examples()
            self.case_subcategory.set(default_case.subcategory)
            self.filter_subcategories()
            if self.editor_path != default_case.path.resolve():
                self.case_choice.set(default_case.label)
                self.choose_example()
        self.root.after_idle(self._center_on_primary_screen)

    def _center_on_primary_screen(self) -> None:
        """Place the window in the primary monitor's usable work area."""
        self.root.attributes("-fullscreen", False)
        self.root.state("normal")
        self.root.update_idletasks()
        work_left = 0
        work_top = 0
        work_right = self.root.winfo_screenwidth()
        work_bottom = self.root.winfo_screenheight()
        if hasattr(ctypes, "windll"):
            rectangle = ctypes.wintypes.RECT()
            # SPI_GETWORKAREA returns the primary monitor area excluding the
            # taskbar, which is the area users perceive as the desktop.
            if ctypes.windll.user32.SystemParametersInfoW(0x0030, 0, ctypes.byref(rectangle), 0):
                work_left, work_top = rectangle.left, rectangle.top
                work_right, work_bottom = rectangle.right, rectangle.bottom
        screen_width = work_right - work_left
        screen_height = work_bottom - work_top
        requested_width = max(self.root.winfo_width(), self.root.winfo_reqwidth())
        requested_height = max(self.root.winfo_height(), self.root.winfo_reqheight())
        width = min(requested_width, max(800, screen_width - 80))
        # Keep an additional visible strip above the work-area boundary so the
        # taskbar is never visually crowded or covered.
        height = min(requested_height, max(600, screen_height - 160))
        self.root.maxsize(max(800, screen_width - 40), max(600, screen_height - 80))
        left = work_left + max(0, (screen_width - width) // 2)
        top = work_top + max(0, (screen_height - height) // 2)
        self.root.geometry(f"{width}x{height}+{left}+{top}")

    def _build(self) -> None:
        # Clam honors custom colors on Windows, including pressed states.
        style = ttk.Style(self.root)
        style.theme_use("clam")
        style.configure(".", foreground="#17212b", background="#f0f0f0",
                        font=("Segoe UI", 10))
        style.map(".", foreground=[("disabled", "#535b63")])
        for control in ("TEntry", "TCombobox", "TSpinbox"):
            style.configure(control, foreground="#17212b", fieldbackground="white")
            style.map(control, foreground=[("disabled", "#535b63")],
                      fieldbackground=[("disabled", "#e5e7eb"), ("readonly", "white")])
        style.configure("TButton", foreground="#17212b", background="#dcecf7",
                        font=("Segoe UI", 10, "bold"), padding=(6, 0), width=0,
                        relief="raised", borderwidth=2)
        # First match wins: pressed must precede hover, disabled precedes both.
        style.map("TButton",
                  foreground=[("disabled", "#535b63"), ("pressed", "white"), ("active", "#17212b")],
                  background=[("disabled", "#e5e7eb"), ("pressed", "#245579"), ("active", "#c6e0f2")],
                  relief=[("disabled", "flat"), ("pressed", "sunken"), ("!pressed", "raised")])
        style.configure("Compact.TButton", padding=(8, 0))
        style.configure("Action.TButton", padding=(8, 1))
        style.configure("Run.Action.TButton", background="#f7b8b8", foreground="#5c1717")
        style.map("Run.Action.TButton",
                  foreground=[("disabled", "#715454"), ("pressed", "white"), ("active", "#5c1717")],
                  background=[("disabled", "#eed4d4"), ("pressed", "#a83232"), ("active", "#f29b9b")])

        outer = ttk.Frame(self.root, padding=10)
        outer.pack(fill=BOTH, expand=True)

        files = ttk.LabelFrame(outer, text="1. Choose an example case or input file", padding=(8, 3))
        files.pack(fill=X)
        families = (SELECT_EXAMPLE,) + THESIS_CATEGORIES + (CUSTOM_CATEGORY,)
        ttk.Label(files, text="Category:").grid(row=0, column=0, sticky="w")
        family_box = ttk.Combobox(files, textvariable=self.case_filter, values=families, state="readonly", width=52)
        self.filter_box = family_box
        family_box.grid(row=0, column=1, sticky="ew", padx=6)
        family_box.bind("<<ComboboxSelected>>", self.filter_examples)
        ttk.Label(files, text="Sub-Category:").grid(row=1, column=0, sticky="w", pady=(2, 0))
        self.subcategory_box = ttk.Combobox(files, textvariable=self.case_subcategory,
                                            values=(SELECT_SUBCATEGORY,), state="readonly", width=52)
        self.subcategory_box.grid(row=1, column=1, sticky="ew", padx=6, pady=(2, 0))
        self.subcategory_box.bind("<<ComboboxSelected>>", self.filter_subcategories)
        ttk.Label(files, text="Example:").grid(row=2, column=0, sticky="w", pady=(2, 0))
        self.case_box = ttk.Combobox(files, textvariable=self.case_choice, state="readonly")
        self.case_box.grid(row=2, column=1, sticky="ew", padx=6, pady=(2, 0))
        self.case_box.bind("<<ComboboxSelected>>", self.choose_example)
        ttk.Button(files, style="Compact.TButton", text="Browse DAT...", command=self.choose_input).grid(row=2, column=2, pady=(2, 0))
        ttk.Label(files, text="Input:").grid(row=3, column=0, sticky="w", pady=(2, 0))
        ttk.Entry(files, textvariable=self.input_path).grid(row=3, column=1, sticky="ew", padx=6, pady=(2, 0))
        ttk.Label(files, text="Result:").grid(row=4, column=0, sticky="w", pady=(2, 0))
        ttk.Entry(files, textvariable=self.output_path).grid(row=4, column=1, sticky="ew", padx=6, pady=(2, 0))
        ttk.Button(files, style="Compact.TButton", text="Save as...", command=self.choose_output).grid(row=4, column=2, pady=(2, 0))
        files.columnconfigure(1, weight=1)


        actions = ttk.Frame(outer, padding=(0, 3))
        actions.pack(fill=X)
        button_opts = dict(style="Action.TButton", cursor="hand2")
        self.run_button = ttk.Button(actions, text="RUN", command=self.run_case,
                                     style="Run.Action.TButton", cursor="hand2")
        self.run_button.pack(side=LEFT, padx=(0, 5))
        self.stop_button = ttk.Button(actions, text="STOP", command=self._cancel_run, state="disabled", **button_opts)
        self.stop_button.pack(side=LEFT, padx=5)
        self.report_button = ttk.Button(actions, text="Report", command=self.export_report, state="disabled", **button_opts)
        self.report_button.pack(side=LEFT, padx=5)
        self.csv_button = ttk.Button(actions, text="Export CSV…", command=self.export_csv, state="disabled", **button_opts)
        self.csv_button.pack(side=LEFT, padx=5)
        ttk.Button(actions, text="Open existing result…", command=self.open_result, **button_opts).pack(side=LEFT, padx=(5, 8))
        self.progress = ttk.Progressbar(actions, mode="indeterminate")
        self.progress.pack(side=RIGHT, fill=X, expand=True, padx=(20, 0))
        ttk.Label(outer, textvariable=self.run_stage, anchor="w").pack(fill=X, pady=(0, 3))

        body = ttk.Panedwindow(outer, orient="horizontal")
        body.pack(fill=BOTH, expand=True)
        left = ttk.Frame(body, padding=(0, 0, 8, 0))
        right = ttk.Frame(body)
        body.add(left, weight=1)
        body.add(right, weight=3)

        tabs = ttk.Notebook(left)
        tabs.pack(fill=BOTH, expand=True)
        inputs = ttk.Frame(tabs, padding=6)
        messages = ttk.Frame(tabs)
        tabs.add(inputs, text="Case inputs")
        tabs.add(messages, text="Messages")
        input_canvas = Canvas(inputs, highlightthickness=0, borderwidth=0, confine=True)
        self.input_viewport = input_canvas
        input_scroll = ttk.Scrollbar(inputs, orient="vertical", command=input_canvas.yview)
        input_canvas.configure(yscrollcommand=input_scroll.set)
        input_scroll.pack(side=RIGHT, fill=Y)
        input_canvas.pack(side=LEFT, fill=BOTH, expand=True)
        input_form = ttk.Frame(input_canvas, padding=6)
        input_window = input_canvas.create_window((0, 0), window=input_form, anchor="nw")
        input_form.bind("<Configure>", lambda e: input_canvas.configure(scrollregion=input_canvas.bbox("all")))
        input_canvas.after_idle(lambda: input_canvas.configure(scrollregion=input_canvas.bbox("all")))
        input_canvas.bind("<Configure>", lambda e: input_canvas.itemconfigure(input_window, width=e.width))
        def scroll_units(step):
            first, last = input_canvas.yview()
            if not ((step < 0 and first <= 0.0) or (step > 0 and last >= 1.0)):
                input_canvas.yview_scroll(step, "units")
            return "break"

        def scroll_inputs(event):
            if not event.delta:
                return "break"
            return scroll_units(-max(1, abs(event.delta) // 120) if event.delta > 0
                                else max(1, abs(event.delta) // 120))

        def bind_input_wheel(widget):
            widget.bind("<MouseWheel>", scroll_inputs)
            widget.bind("<Button-4>", lambda event: scroll_units(-1))
            widget.bind("<Button-5>", lambda event: scroll_units(1))
            for child in widget.winfo_children():
                bind_input_wheel(child)

        ttk.Label(input_form, textvariable=self.input_reference, wraplength=300,
                  anchor="w", justify="left").pack(fill=X)
        cells = ttk.Frame(input_form)
        cells.pack(fill=X, pady=6)
        self.geometry_button = ttk.Button(cells, text="Edit boundaries…", command=self._edit_geometry)
        self.geometry_button.grid(row=0, column=0, sticky="w")
        style.configure("Boundary.TEntry", foreground="#17212b")
        style.map("Boundary.TEntry", foreground=[("readonly", "#17212b")],
                  fieldbackground=[("readonly", "#f0f0f0")])
        self.boundaries_entry = ttk.Entry(cells, textvariable=self.input_boundaries,
                                          state="readonly", width=12, style="Boundary.TEntry")
        self.boundaries_entry.grid(row=0, column=1, sticky="ew", padx=(6, 0))
        self.input_entries = {}
        self.input_labels = {}
        style = ttk.Style(self.root)
        style.map("Input.TEntry", foreground=[("invalid", "#b00020")],
                  fieldbackground=[("invalid", "#fff0f0")])
        display_row = 1
        for key, label in INPUT_FIELDS:
            row = display_row
            label_widget = ttk.Label(cells, text=label)
            label_widget.grid(row=row, column=0, sticky="w", pady=2)
            self.input_labels[key] = label_widget
            entry = ttk.Entry(cells, textvariable=self.input_values[key], width=12, style="Input.TEntry")
            entry.grid(row=row, column=1, sticky="ew", padx=(6, 0), pady=2)
            self.input_entries[key] = entry
            self.input_values[key].trace_add("write", self._input_changed)
            display_row += 1
            if key == "formula":
                explanation = ttk.Label(cells, text="1: cubic root   2: square root   3: other type",
                                        foreground="#555555")
                explanation.grid(row=display_row, column=0, sticky="w", padx=(12, 0), pady=(0, 3))
                self.formula_explanation = explanation
                display_row += 1
        cells.columnconfigure(1, weight=1)
        self.reset_inputs_button = ttk.Button(input_form, text="Reset inputs", command=self._reset_inputs)
        self.reset_inputs_button.pack(anchor="w")
        description = ttk.Label(input_form, textvariable=self.input_description, wraplength=300, justify="left")
        self.input_description_label = description
        description.pack(fill=X, pady=(8, 0))
        ttk.Label(input_form, textvariable=self.input_full_details, wraplength=300, justify="left").pack(fill=X, pady=(8, 0))
        input_form.bind("<Configure>", lambda event: description.configure(wraplength=max(180, event.width - 12)), add="+")
        bind_input_wheel(input_canvas)
        self.log = ScrolledText(messages, width=38, wrap="word", state="disabled")
        self.log.pack(fill=BOTH, expand=True, pady=(4, 0))

        controls = ttk.LabelFrame(right, text="2. Results", padding=(6, 3))
        self.results_frame = controls
        controls.pack(fill=X)
        self.results_case_label = ttk.Label(controls, textvariable=self.results_title,
                                            font=("Segoe UI", 10, "bold"))
        self.results_case_label.pack(anchor="w", pady=(0, 2))
        # Retain the preset selection model for automatic example setup, but
        # leave it out of the layout: the visible controls edit the same view.
        self.thesis_graph_box = ttk.Combobox(controls, textvariable=self.thesis_graph,
                                             state="disabled", width=78)
        self.thesis_graph_box.bind("<<ComboboxSelected>>", self._select_thesis_graph)
        selectors = ttk.Frame(controls)
        selectors.pack(fill=X, pady=(0, 2))
        ttk.Label(selectors, text="View ( X ):").pack(side=LEFT)
        self.view_box = ttk.Combobox(selectors, textvariable=self.view_label, width=36, state="readonly")
        self.view_box.pack(side=LEFT, padx=(4, 12))
        self.view_box.bind("<<ComboboxSelected>>", self._select_view)
        ttk.Label(selectors, text="Quantity ( Y ):").pack(side=LEFT)
        self.quantity_box = ttk.Combobox(selectors, textvariable=self.quantity_label, width=34, state="readonly")
        self.quantity_box.pack(side=LEFT, padx=(4, 12))
        self.quantity_box.bind("<<ComboboxSelected>>", self._select_quantity)
        options = ttk.Frame(controls)
        options.pack(fill=X)
        controls = options
        ttk.Label(controls, text="Node:").pack(side=LEFT)
        self.node_box = ttk.Spinbox(controls, from_=1, to=1, textvariable=self.node, width=5, command=self._select_single_node)
        self.node_box.pack(side=LEFT, padx=(4, 12))
        self.node_box.bind("<KeyRelease>", self._single_node_changed)
        self.node_box.bind("<FocusOut>", self._request_plot)
        ttk.Label(controls, text="Nodes:").pack(side=LEFT)
        self.nodes_entry = ttk.Entry(controls, textvariable=self.nodes_text, width=18)
        self.nodes_entry.pack(side=LEFT, padx=(4, 12))
        self.nodes_entry.bind("<KeyRelease>", self._multiple_nodes_changed)
        self.nodes_entry.bind("<FocusIn>", self._multiple_nodes_changed)
        ttk.Button(controls, text="DAT nodes", command=self._restore_dat_nodes).pack(side=LEFT, padx=(0, 4))
        ttk.Button(controls, text="All nodes", command=self._restore_all_nodes).pack(side=LEFT, padx=(0, 12))
        ttk.Label(controls, text="Component:").pack(side=LEFT)
        component = ttk.Combobox(controls, textvariable=self.component, values=("X", "Y"), width=4, state="readonly")
        self.component_box = component
        component.pack(side=LEFT, padx=(4, 12))
        component.bind("<<ComboboxSelected>>", self._request_plot)
        ttk.Label(controls, text="Step:").pack(side=LEFT)
        self.step_box = ttk.Spinbox(controls, from_=1, to=1, textvariable=self.step, width=5, command=self._request_plot)
        self.step_box.pack(side=LEFT, padx=(4, 12))
        self.step_box.bind("<KeyRelease>", self._request_plot)
        self.step_box.bind("<FocusOut>", self._request_plot)
        ttk.Button(controls, text="Update graph", command=self._manual_update_graph).pack(side=LEFT)

        zoom_controls = ttk.Frame(right)
        zoom_controls.pack(fill=X, pady=(3, 0))
        ttk.Label(zoom_controls, text="Graph size:").pack(side=LEFT)
        self.zoom_out_button = ttk.Button(zoom_controls, text="Zoom out -",
                                          command=lambda: self._set_graph_zoom(self.graph_zoom - 0.25))
        self.zoom_out_button.pack(side=LEFT, padx=(4, 0))
        ttk.Label(zoom_controls, textvariable=self.graph_zoom_label, width=6, anchor="center").pack(side=LEFT)
        self.zoom_in_button = ttk.Button(zoom_controls, text="Zoom in +",
                                         command=lambda: self._set_graph_zoom(self.graph_zoom + 0.25))
        self.zoom_in_button.pack(side=LEFT)
        ttk.Button(zoom_controls, text="Reset size", command=lambda: self._set_graph_zoom(1.0)).pack(side=LEFT, padx=4)

        self.figure = Figure(figsize=(7, 5), dpi=100)
        self.axes = self.figure.add_subplot(111)
        graph_frame = ttk.Frame(right)
        graph_frame.pack(fill=BOTH, expand=True, pady=(3, 0))
        graph_frame.rowconfigure(0, weight=1)
        graph_frame.columnconfigure(0, weight=1)
        self.graph_viewport = Canvas(graph_frame, highlightthickness=0,
                                     width=700, height=500, xscrollincrement=24, yscrollincrement=24)
        self.graph_viewport.grid(row=0, column=0, sticky="nsew")
        vertical = ttk.Scrollbar(graph_frame, orient="vertical", command=self.graph_viewport.yview)
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal = ttk.Scrollbar(graph_frame, orient="horizontal", command=self.graph_viewport.xview)
        horizontal.grid(row=1, column=0, sticky="ew")
        self.graph_viewport.configure(xscrollcommand=horizontal.set, yscrollcommand=vertical.set)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.graph_viewport)
        graph_widget = self.canvas.get_tk_widget()
        self.graph_window = self.graph_viewport.create_window(0, 0, anchor="nw", window=graph_widget,
                                                              width=700, height=500)
        self.graph_viewport.configure(scrollregion=(0, 0, 700, 500))
        self.graph_viewport.bind("<Configure>", self._resize_graph_viewport)
        for widget in (self.graph_viewport, graph_widget):
            widget.bind("<MouseWheel>", self._scroll_graph)
            widget.bind("<Shift-MouseWheel>", self._scroll_graph)
            widget.bind("<Control-MouseWheel>", self._zoom_graph_wheel)
        self.axes.set_title("Select an input case or open an existing result")
        self.axes.grid(True, alpha=0.3)
        self.canvas.draw()
        self.filter_examples()

    def _set_graph_zoom(self, zoom: float) -> None:
        """Scale the whole figure, including text, without changing the data limits."""
        self.graph_zoom = round(max(0.5, min(3.0, zoom)), 2)
        self.graph_zoom_label.set(f"{self.graph_zoom:.0%}")
        self.zoom_out_button.configure(state="disabled" if self.graph_zoom <= 0.5 else "normal")
        self.zoom_in_button.configure(state="disabled" if self.graph_zoom >= 3.0 else "normal")
        self.figure.set_dpi(100 * self.graph_zoom)
        self._resize_graph_viewport()
        self.canvas.draw_idle()

    def _zoom_graph_wheel(self, event: object) -> str:
        if event.delta:
            self._set_graph_zoom(self.graph_zoom + (0.25 if event.delta > 0 else -0.25))
        return "break"

    def _resize_graph_viewport(self, event: object | None = None) -> None:
        """Keep the chosen scale when resizing the viewport; scroll larger graphs."""
        viewport_width = event.width if event is not None else self.graph_viewport.winfo_width()
        viewport_height = event.height if event is not None else self.graph_viewport.winfo_height()
        width = round(max(700, viewport_width) * self.graph_zoom)
        height = round(max(500, viewport_height) * self.graph_zoom)
        self.graph_viewport.itemconfigure(self.graph_window, width=width, height=height)
        self.graph_viewport.configure(scrollregion=(0, 0, width, height))
        if viewport_width >= width:
            self.graph_viewport.xview_moveto(0)
        if viewport_height >= height:
            self.graph_viewport.yview_moveto(0)

    def _scroll_graph(self, event: object) -> str:
        if event.delta:
            units = -int(event.delta / 120) if abs(event.delta) >= 120 else (-1 if event.delta > 0 else 1)
            if event.state & 0x0001:  # Shift + wheel scrolls horizontally.
                self.graph_viewport.xview_scroll(units, "units")
            else:
                self.graph_viewport.yview_scroll(units, "units")
        return "break"

    def _select_view(self, _event: object | None = None) -> None:
        self.view.set(RESULT_VIEWS[self.view_box.current()])
        self.thesis_history_nodes = None
        self.thesis_graph.set("")
        self._refresh_result_choices()
        if self.result is not None:
            self.plot_result()

    def _thesis_graph_options(self) -> tuple[tuple[object, ...], ...]:
        if self.active_case is None:
            return ()
        return tuple((f"Figure {graph.figure} - {graph.caption}", graph.view, graph.quantity,
                      graph.component, graph.step, graph.nodes)
                     for graph in thesis_graphs_for_case(self.active_case.path.stem, self.active_deck))

    def _refresh_thesis_graph_choices(self) -> None:
        options = self._thesis_graph_options()
        labels = tuple(option[0] for option in options)
        self.thesis_graph_box.configure(values=labels, state="readonly" if labels else "disabled")
        if self.thesis_graph.get() not in labels:
            if not options:
                self.thesis_graph.set("")
                self.thesis_history_nodes = None
                return
            current = (self.view.get(), self.quantity.get(), self.component.get())
            selected = next((option for option in options if option[1:4] == current), options[0])
            self.thesis_graph.set(selected[0])
            # Apply the preset without drawing a previous case's result while loading.
            self._select_thesis_graph(redraw=False)

    def _select_thesis_graph(self, _event: object | None = None, *, redraw: bool = True) -> None:
        selected = next((option for option in self._thesis_graph_options() if option[0] == self.thesis_graph.get()), None)
        if selected is None:
            return
        _, view, quantity, component, featured_step, nodes = selected
        self.view.set(view)
        self.quantity.set(quantity)
        self.component.set(component)
        if featured_step is not None:
            self.step.set(str(featured_step))
        self.thesis_history_nodes = nodes
        self._refresh_result_choices()
        if redraw and self.result is not None:
            self.plot_result()

    def _select_quantity(self, _event: object | None = None) -> None:
        self.quantity.set(RESULT_QUANTITIES[self.quantity_box.current()])
        self.thesis_history_nodes = None
        self.thesis_graph.set("")
        self._refresh_result_choices()
        if self.result is not None:
            self.plot_result()

    def _refresh_result_choices(self) -> None:
        """Use the graph caption's figure mapping in both dropdown menus."""
        component = 0 if self.component.get() == "X" else 1
        is_deformed = self.view.get() == "Deformed geometry"
        self.component_box.configure(state="disabled" if self.quantity.get() in STRESS_COMPONENTS or is_deformed else "readonly")
        self.quantity_box.configure(state="disabled" if is_deformed else "readonly")
        # A time-history plot uses all steps.  Thesis presets may also specify
        # a fixed set of history nodes, so the manual Node control is then not
        # an input to that graph.  Step selects the snapshot for the other
        # views only.
        node_locked = self.view.get() == "Time history" and bool(self.thesis_history_nodes)
        self.node_box.configure(state="disabled" if node_locked else "normal")
        self.step_box.configure(state="disabled" if self.view.get() == "Time history" else "normal")

        def label(name: str, view: str, quantity: str) -> str:
            caption = self._caption(view, quantity, component)
            if "unavailable" in caption or "no separate PhD figure" in caption:
                return name
            else:
                reference = "Fig. " + caption.split()[1]
            return f"{name} — {reference}"

        views = tuple(label(view, view, self.quantity.get()) for view in RESULT_VIEWS)
        quantities = tuple(label(quantity, self.view.get(), quantity) for quantity in RESULT_QUANTITIES)
        self.view_box.configure(values=views)
        self.quantity_box.configure(values=quantities)
        if self.case_filter.get() == SELECT_EXAMPLE:
            self.view_label.set("")
            self.quantity_label.set("")
            return
        self.view_label.set(views[RESULT_VIEWS.index(self.view.get())])
        self.quantity_label.set(quantities[RESULT_QUANTITIES.index(self.quantity.get())])

    def _request_plot(self, _event: object | None = None) -> None:
        """Debounce control changes and redraw as soon as values are valid."""
        self._refresh_result_choices()
        if self.result is None or not self.node.get().isdigit() or not self.step.get().isdigit():
            return
        if self._plot_after_id is not None:
            self.root.after_cancel(self._plot_after_id)
        self._plot_after_id = self.root.after(80, self._run_requested_plot)

    def _select_single_node(self) -> None:
        self.nodes_text.set("")
        self._refresh_result_choices()
        self._request_plot()

    def _single_node_changed(self, _event: object | None = None) -> str:
        """Typing a single node opts out of the multi-node selection."""
        if self.node.get().strip():
            self.thesis_history_nodes = None
            self.nodes_text.set("")
        self._request_plot()
        return "break" if _event is None else ""

    def _multiple_nodes_changed(self, _event: object | None = None) -> str:
        """Typing a node list clears the single-node selection."""
        if self.nodes_text.get().strip():
            self.thesis_history_nodes = None
            self.node.set("")
        self._request_plot()
        return ""

    def _restore_dat_nodes(self) -> None:
        if self.active_deck is None:
            return
        nodes = tuple(int(n) for n in self.active_deck.node_plots)
        if self.result is not None:
            nodes = tuple(n for n in nodes if 1 <= n <= len(self.result["coordinates"]))
        if nodes:
            self.node.set("")
            self.thesis_history_nodes = None
            self.nodes_text.set(", ".join(map(str, nodes)))
            self.nodes_entry.focus_set()
            self._request_plot()

    def _restore_all_nodes(self) -> None:
        if self.result is not None:
            node_count = len(self.result["coordinates"])
        elif self.active_deck is not None:
            try:
                meshes = tuple(mesh_section(section) for section in self.active_deck.sections)
                node_count = len(combine_meshes(meshes).coordinates)
            except (ValueError, TypeError):
                return
        else:
            return
        self.node.set("")
        self.thesis_history_nodes = None
        self.nodes_text.set(", ".join(str(i) for i in range(1, node_count + 1)))
        self.nodes_entry.focus_set()
        self._request_plot()

    def _run_requested_plot(self) -> None:
        self._plot_after_id = None
        self.plot_result()

    def _manual_update_graph(self) -> None:
        """Explicitly apply all graph controls and report the refresh."""
        # Combobox labels include the optional figure suffix; synchronize the
        # model from what is visibly selected before plotting.
        view_label = self.view_label.get()
        quantity_label = self.quantity_label.get()
        for index, name in enumerate(RESULT_VIEWS):
            if view_label == name or view_label.startswith(name + " "):
                self.view.set(name)
                break
        for index, name in enumerate(RESULT_QUANTITIES):
            if quantity_label == name or quantity_label.startswith(name + " "):
                self.quantity.set(name)
                break
        self._refresh_result_choices()
        self.plot_result()
        self.write_log(
            f"Graph updated: view={self.view.get()}, quantity={self.quantity.get()}, "
            f"nodes={self.nodes_text.get() or self.node.get()}, step={self.step.get()}"
        )

    def _update_run_button(self) -> None:
        disabled = (self._running or self.case_filter.get() == SELECT_EXAMPLE
                    or self.case_subcategory.get() == SELECT_SUBCATEGORY or not self._inputs_valid)
        self.run_button.configure(state="disabled" if disabled else "normal")

    def _load_inputs(self, deck: InputDeck, path: Path, case: ExampleCase | None = None) -> None:
        self.editor_deck, self.editor_path, self.editor_case = deck, path.resolve(), case
        self._source_editor_deck, self._source_editor_case = deck, case
        self._loading_inputs = True
        self.geometry_button.configure(state="normal")
        self.boundaries_entry.configure(state="readonly")
        self.reset_inputs_button.configure(state="normal")
        self.input_reference.set(f"Thesis test: {case.thesis_section}\nCategory: {case.family}" if case else "Custom input file")
        title = " ".join(deck.title.split()) or path.stem
        self.results_title.set(title)
        details = asdict(deck)
        self.input_full_details.set("All parsed DAT inputs:\n" + _format_input_details(details))
        self.input_boundaries.set(str(len(deck.sections)))
        values = dict(steps=deck.steps, total_time=deck.total_time, young=deck.material.young,
                      poisson=deck.material.poisson, density=deck.material.density)
        if isinstance(deck.load, PlaneWave):
            values.update(p=deck.load.longitudinal_stress, sv=deck.load.transverse_shear, angle=deck.load.angle_degrees)
        elif isinstance(deck.load, ExplosionPoint):
            # The common input panel uses P/SV fields for the relative
            # explosion amplitudes when an EXPLOSION.POINT deck is loaded.
            values.update(p=deck.load.p_amplitude, sv=deck.load.sv_amplitude,
                          load_x=deck.load.x, load_y=deck.load.y,
                          wave_velocity=deck.load.wave_velocity, formula=deck.load.formula,
                          coefficient=deck.load.coefficient, decay=deck.load.decay,
                          coupling=deck.load.frequency)
        if len(deck.sections) == 1 and deck.sections[0].shape in (1, 2, 3):
            # Geometry[0] is radius for a circle/ellipse and the primary
            # side (X dimension) for square/rectangle sections.
            values["radius"] = deck.sections[0].geometry[0]
            if deck.sections[0].shape == 3 and len(deck.sections[0].geometry) > 1:
                values["height"] = deck.sections[0].geometry[1]
        shape = deck.sections[0].shape if deck.sections else None
        visible_geometry = {"radius"} if shape in (1, 2) else ({"radius", "height"} if shape == 3 else set())
        visible_explosion = isinstance(deck.load, ExplosionPoint)
        if hasattr(self, "formula_explanation"):
            (self.formula_explanation.grid if visible_explosion else self.formula_explanation.grid_remove)()
        for key in ("height", "load_x", "load_y", "wave_velocity", "formula", "coefficient", "decay", "coupling"):
            label_widget = self.input_labels.get(key); entry_widget = self.input_entries.get(key)
            if label_widget is None or entry_widget is None:
                continue
            show = key in visible_geometry or (visible_explosion and key in {"load_x", "load_y", "wave_velocity", "formula", "coefficient", "decay", "coupling"})
            (label_widget.grid if show else label_widget.grid_remove)()
            (entry_widget.grid if show else entry_widget.grid_remove)()
        if "radius" in self.input_labels:
            self.input_labels["radius"].configure(
                text={1: "Cavity radius (m)", 2: "Square side length (m)", 3: "Rectangle X-side length (m)"}.get(shape, "Cavity size (m)")
            )
        for key, variable in self.input_values.items():
            variable.set(f"{values[key]:g}" if key in values else "")
            self.input_entries[key].configure(state="normal" if key in values else "disabled")
        # Load the DAT multi-node preset as soon as the case is selected.
        if deck.node_plots:
            self.nodes_text.set(", ".join(str(node) for node in deck.node_plots))
            self.node.set("")
        else:
            self.nodes_text.set("")
            self.node.set("1")
        self._loading_inputs = False
        self._input_changed()

    def _clear_input_cells(self) -> None:
        """Clear the panel while leaving fields available for direct entry."""
        self.editor_deck = self.editor_path = self.editor_case = None
        self.input_reference.set("No input file selected")
        self.results_title.set("No test selected")
        self.input_boundaries.set("—")
        self.input_description.set("")
        self.input_full_details.set("")
        self.input_description_label.configure(foreground="")
        self._inputs_valid = False
        self.boundaries_entry.configure(state="readonly")
        self.reset_inputs_button.configure(state="disabled")
        self.geometry_button.configure(state="disabled")
        defaults = {
            "steps": 20, "total_time": 0.015, "young": 20000,
            "poisson": 0.25, "density": 2.55, "p": -10, "sv": 0,
            "angle": 0, "radius": 5, "height": 4, "load_x": 8,
            "load_y": 0, "wave_velocity": 1850.50542, "formula": 1,
            "coefficient": 160, "decay": 0, "coupling": 1,
        }
        for key, variable in self.input_values.items():
            variable.set(f"{defaults[key]:g}" if key in defaults else "")
            self.input_entries[key].configure(state="normal")
            self.input_entries[key].state(["!invalid"])
        self._update_run_button()

    def _reset_inputs(self) -> None:
        if self.editor_deck is not None and not self._running and self.case_filter.get() != SELECT_EXAMPLE:
            self._load_inputs(self._source_editor_deck, self.editor_path, self._source_editor_case)
            self.active_deck, self.active_case = self.editor_deck, self.editor_case
            self.active_thesis_section = self.editor_case.thesis_section if self.editor_case else ""
            self._refresh_thesis_graph_choices()
            self._clear_result_graph()

    def _edit_geometry(self) -> None:
        if self._running or self.editor_deck is None:
            return
        try:
            # Preserve material/loading edits when the geometry dialog applies.
            deck = edited_input_deck(self.editor_deck, {
                key: value.get() for key, value in self.input_values.items()
                if not self.input_entries[key].instate(["disabled"])
            })
        except ValueError as error:
            messagebox.showerror("Input error", str(error), parent=self.root)
            return
        GeometryEditor(self.root, deck, self._apply_geometry)

    def _apply_geometry(self, deck: InputDeck) -> None:
        source, source_case = self._source_editor_deck, self._source_editor_case
        self._load_inputs(deck, self.editor_path)
        self._source_editor_deck, self._source_editor_case = source, source_case
        self.input_reference.set("Custom geometry (edited inputs)")
        self.active_deck, self.active_case = deck, None
        self.active_thesis_section = ""
        self.thesis_history_nodes = None
        self.nodes_text.set(", ".join(str(node) for node in deck.node_plots))
        self.node.set("" if deck.node_plots else "1")
        self._modified_inputs = True
        self._refresh_thesis_graph_choices()
        self._clear_result_graph()
        self.write_log(f"Geometry updated: {len(deck.sections)} boundaries. Press RUN to calculate.")

    def _input_changed(self, *_args: object) -> None:
        if self._loading_inputs or self.editor_deck is None or self.case_filter.get() == SELECT_EXAMPLE:
            return
        values = {key: value.get() for key, value in self.input_values.items()
                  if not self.input_entries[key].instate(["disabled"])}
        _, errors = validate_input_numbers(self.editor_deck, values)
        self._inputs_valid = not errors
        for key, entry in self.input_entries.items():
            entry.state(["invalid" if key in errors else "!invalid"])
        self._update_run_button()
        self.input_description_label.configure(foreground="#b00020" if errors else "")
        if errors:
            self.input_description.set("Please correct these inputs:\n" + "\n".join(errors.values()))
            return
        self.input_description.set("")

    def filter_examples(self, _event: object | None = None) -> None:
        selected = self.case_filter.get()
        self._update_run_button()
        subcategories = (CUSTOM_SUBCATEGORY,) if selected == CUSTOM_CATEGORY else tuple(dict.fromkeys(
            case.subcategory for case in self.examples if case.family == selected
        )) or THESIS_SUBCATEGORIES.get(selected, ())
        self.subcategory_box.configure(values=(SELECT_SUBCATEGORY,) + subcategories)
        self.case_subcategory.set(SELECT_SUBCATEGORY)
        self.filtered_examples = ()
        self._clear_case_selection()

    def filter_subcategories(self, _event: object | None = None) -> None:
        selected = self.case_subcategory.get()
        self.filtered_examples = tuple(case for case in self.examples if case.family == self.case_filter.get()
                                       and case.subcategory == selected)
        self.case_box.configure(values=tuple(case.label for case in self.filtered_examples))
        if self.filtered_examples:
            self.case_choice.set(self.filtered_examples[0].label)
            self.choose_example()
        else:
            self._clear_case_selection()
        self._update_run_button()

    def _clear_case_selection(self) -> None:
        self.case_box.configure(values=())
        self.case_choice.set("")
        self.input_path.set("")
        self.output_path.set("")
        self.thesis_graph.set("")
        self.thesis_graph_box.configure(values=(), state="disabled")
        self.active_case = None
        self.thesis_history_nodes = None
        self._clear_input_cells()
        self._clear_result_graph()

    def choose_example(self, _event: object | None = None) -> None:
        label = self.case_choice.get()
        case = next((item for item in self.filtered_examples if item.label == label), None)
        if case is None:
            return
        previous_path = self.editor_path.resolve() if self.editor_path is not None else None
        # Preserve controls only for a true re-selection of an already loaded
        # case with a displayed result; initial selector population must still
        # apply the case's thesis preset.
        same_case = self.active_case is not None and self.result is not None and previous_path == case.path.resolve()
        previous_controls = (self.view.get(), self.quantity.get(), self.component.get(), self.node.get(), self.step.get())
        self.input_path.set(str(case.path))
        self.active_deck = case.deck
        self.active_case = case
        self.active_thesis_section = case.thesis_section
        self._load_inputs(case.deck, case.path, case)
        self._refresh_thesis_graph_choices()
        self._refresh_result_choices()
        self.output_path.set(str(Path("examples/results/npz") / f"{case.path.stem}_result.npz"))
        self.write_log(case_details(case), clear=True)
        if same_case:
            self.view.set(previous_controls[0]); self.quantity.set(previous_controls[1]); self.component.set(previous_controls[2])
            self.node.set(previous_controls[3]); self.step.set(previous_controls[4])
            self.thesis_graph.set(""); self._refresh_result_choices()
            if self.result is not None:
                self.plot_result()
        else:
            self._clear_result_graph()

    def _clear_result_graph(self) -> None:
        """Avoid showing a graph from a different case before this case is run."""
        self.result = None
        self.result_path = None
        self.csv_button.configure(state="disabled")
        if getattr(self, "node_inset", None) is not None:
            self.node_inset.remove()
            self.node_inset = None
        self.axes.clear()
        self.axes.set_title("Press Run selected case to display this graph")
        self.axes.grid(True, alpha=0.3)
        self.figure.texts.clear()
        self.canvas.draw()

    def write_log(self, text: str, clear: bool = False) -> None:
        self.log.configure(state="normal")
        if clear:
            self.log.delete("1.0", END)
        self.log.insert(END, text.rstrip() + "\n")
        self.log.see(END)
        self.log.configure(state="disabled")

    def choose_input(self) -> None:
        default_dir = self.last_input_dir if self.last_input_dir.is_dir() else example_directory()
        path = filedialog.askopenfilename(
            title="Select PREBEM input",
            initialdir=str(default_dir),
            filetypes=(("PREBEM data", "*.DAT"), ("All files", "*.*")),
        )
        if not path:
            return
        self.last_input_dir = Path(path).resolve().parent
        try:
            self.last_input_dir_file.parent.mkdir(parents=True, exist_ok=True)
            self.last_input_dir_file.write_text(str(self.last_input_dir), encoding="utf-8")
        except OSError:
            pass
        try:
            deck = read_input_deck(path)
            case = identify_input_case(Path(path), deck, self.examples)
            self.examples = tuple(item for item in self.examples if item.path.resolve() != case.path.resolve()) + (case,)
            self.case_filter.set(case.family)
            self.filter_box.configure(values=(SELECT_EXAMPLE,) + tuple(dict.fromkeys(
                THESIS_CATEGORIES + tuple(item.family for item in self.examples) + (CUSTOM_CATEGORY,))))
            self.subcategory_box.configure(values=(SELECT_SUBCATEGORY,) + tuple(dict.fromkeys(
                item.subcategory for item in self.examples if item.family == case.family)))
            self.case_subcategory.set(case.subcategory)
            self.filtered_examples = tuple(item for item in self.examples
                                           if item.family == case.family and item.subcategory == case.subcategory)
            self.case_box.configure(values=tuple(item.label for item in self.filtered_examples))
            self.case_choice.set(case.label)
            self.input_path.set(path)
            self.output_path.set(str(Path("examples/results/npz") / f"{Path(path).stem}_result.npz"))
            self.active_case = case
            self.active_deck = deck
            self.active_thesis_section = case.thesis_section if case else ""
            self._load_inputs(deck, Path(path), case)
            self._refresh_thesis_graph_choices()
            self._refresh_result_choices()
            self.write_log(_deck_summary(deck, Path(path)), clear=True)
            self._clear_result_graph()
            self._show_dat_details(Path(path), deck, case)
        except Exception as error:
            self.write_log(f"Input cannot be read:\n{error}", clear=True)

    def _show_dat_details(self, path: Path, deck: InputDeck, case: ExampleCase | None) -> None:
        """Show DAT-derived metadata, parsed inputs, and messages in a separate form."""
        old = getattr(self, "dat_details_window", None)
        if old is not None and old.winfo_exists():
            old.destroy()
        window = Toplevel(self.root)
        self.dat_details_window = window
        window.title(f"DAT details — {path.name}")
        window.geometry("720x620")
        window.transient(self.root)
        # Reserve the footer before packing expanding content.  Otherwise the
        # notebook's requested height can squeeze Continue out at high DPI.
        actions = ttk.Frame(window)
        actions.pack(side="bottom", fill=X, padx=10, pady=(0, 10))
        ttk.Button(actions, text="Continue", command=window.destroy).pack(side="right")
        summary = ttk.LabelFrame(window, text="Case identification", padding=8)
        summary.pack(fill=X, padx=10, pady=10)
        category = case.family if case else ("1. ELASTOSTATIC DOMAIN" if deck.domain == 0 else "2. LAPLACE DOMAIN" if deck.domain == 1 else "3. TIME DOMAIN" if deck.domain == 2 else "5. MULTI-ZONAL BEM")
        subcategory = case.subcategory if case else "Custom DAT input"
        title = deck.title.strip() or path.stem
        for row, (label, value) in enumerate((("Category", category), ("Sub-category", subcategory), ("Example", title), ("Input", str(path)))):
            ttk.Label(summary, text=f"{label}:", font=("TkDefaultFont", 9, "bold")).grid(row=row, column=0, sticky="nw", padx=(0, 8), pady=2)
            ttk.Label(summary, text=value, wraplength=600, justify="left").grid(row=row, column=1, sticky="nw", pady=2)
        tabs = ttk.Notebook(window)
        tabs.pack(fill=BOTH, expand=True, padx=10, pady=(0, 10))
        inputs_frame = ttk.Frame(tabs)
        tabs.add(inputs_frame, text="Case inputs")
        inputs = ttk.Treeview(inputs_frame, columns=("field", "value", "meaning"), show="headings")
        inputs.heading("field", text="DAT field")
        inputs.heading("value", text="Value")
        inputs.heading("meaning", text="Meaning")
        inputs.tag_configure("group", background="#e8eef7", font=("TkDefaultFont", 9, "bold"))
        inputs.column("field", width=210, anchor="w", stretch=False)
        inputs.column("value", width=150, anchor="w", stretch=False)
        inputs.column("meaning", width=520, anchor="w", stretch=True)
        def rows(mapping: dict[str, object], prefix: str = ""):
            for key, value in mapping.items():
                name = f"{prefix}.{key}" if prefix else key
                if key == "title":
                    continue
                if isinstance(value, dict):
                    yield from rows(value, name)
                elif isinstance(value, (list, tuple)):
                    for index, item in enumerate(value):
                        item_name = f"{name}[{index}]"
                        if isinstance(item, dict):
                            yield from rows(item, item_name)
                        elif isinstance(item, (list, tuple)):
                            for sub_index, sub_item in enumerate(item):
                                yield f"{item_name}[{sub_index}]", sub_item
                        else:
                            yield item_name, item
                else:
                    yield name, value
        parsed_rows = [("load.module", "EXPLOSION.POINT" if isinstance(deck.load, ExplosionPoint) else "PLANE.WAVE")]
        parsed_rows.extend(rows(asdict(deck)))
        # These optional modules may be omitted from a DAT file.  Keep them
        # visible so the user can see the applicable default and meaning.
        parsed_keys = {key for key, _ in parsed_rows}
        if "free_field_stresses" not in parsed_keys:
            parsed_rows.append(("free_field_stresses", None))
        if "simulation" not in parsed_keys:
            parsed_rows.append(("simulation", None))
        groups = {
            "general": "General problem",
            "geometry": "Geometry",
            "material": "Material",
            "loading": "Loading",
            "output": "Plots and output",
        }
        current_group = None
        for key, value in parsed_rows:
            if key.startswith("sections["):
                group = "geometry"
            elif key.startswith("material"):
                group = "material"
            elif key.startswith("load"):
                group = "loading"
            elif key.startswith("plots") or key.startswith("node_plots") or key in {"element_nodes", "internal_points"}:
                group = "output"
            else:
                group = "general"
            if group != current_group:
                inputs.insert("", "end", values=(f"--- {groups[group]} ---", "", ""), tags=("group",))
                current_group = group
            if key == "units" and value is None:
                shown_value = "1 (default)"
            elif key == "free_field_stresses" and value is None:
                shown_value = "0 (default)"
            elif key == "simulation" and value is None:
                shown_value = "omitted"
            else:
                shown_value = _display_input_value(value)
            meaning = _input_meaning(key, value)
            if ".geometry[" in key:
                section_index = key.split("sections[", 1)[1].split("]", 1)[0]
                shape = asdict(deck)["sections"][int(section_index)]["shape"]
                geometry_meanings = {
                    1: ("Radius / semiaxis (selected length unit)", "Axis ratio b/a (dimensionless)", "Centre X (length)", "Centre Y (length)", "Rotation (degrees)"),
                    2: ("Square side length (selected length unit)", "Centre X (length)", "Centre Y (length)", "Rotation (degrees)", "X-side grading ratio", "Y-side grading ratio", "Corner discontinuity: 1 shared node; 2 duplicate corner nodes"),
                    3: ("Rectangle X-side length (selected length unit)", "Rectangle Y-side length (selected length unit)", "Centre X (length)", "Centre Y (length)", "Rotation (degrees)", "X-side grading ratio: 1 uniform; >1 increasing; 0<p<1 decreasing; p<0 symmetric grading", "Y-side grading ratio: 1 uniform; 0 means calculate automatically; negative values grade toward both corners", "Corner discontinuity: 1 shared corner node; 2 duplicate corner nodes for accurate corner stress"),
                }
                try:
                    meaning = geometry_meanings.get(shape, ())[int(key.rsplit("[", 1)[1].split("]", 1)[0])]
                except (IndexError, ValueError):
                    pass
                length_unit = {1: "m", 2: "imperial length unit", 3: "user-defined length unit"}.get(deck.units or 1, "length unit")
                meaning = meaning.replace("(length)", f"({length_unit})").replace("(selected length unit)", f"({length_unit})")
            if key == "load.module":
                meaning = "Point explosion source (EXPLOSION.POINT)" if value == "EXPLOSION.POINT" else "Plane-wave source (PLANE.WAVE)"
            field_leaf = key.rsplit(".", 1)[-1]
            if key.startswith("load.") and field_leaf in {"x", "y", "wave_velocity"}:
                load_unit = {1: "m", 2: "imperial length unit", 3: "user-defined length unit"}.get(deck.units or 1, "length unit")
                meaning = (f"Explosion source coordinate ({load_unit})" if field_leaf in {"x", "y"}
                           else f"Explosion wave velocity ({load_unit}/s)")
            if key == "load.frequency" and isinstance(deck.load, ExplosionPoint):
                meaning = ("Explosion coupling factor F (dimensionless; F=1.0 means a fully coupled charge "
                           "fills the borehole; decoupled charges have an air gap or smaller charge diameter)")
            inputs.insert("", "end", values=(key, shown_value, meaning))
        inputs_scroll = ttk.Scrollbar(inputs_frame, orient="vertical", command=inputs.yview)
        inputs.configure(yscrollcommand=inputs_scroll.set)
        inputs.place(relx=0, rely=0, relwidth=0.975, relheight=1)
        inputs_scroll.place(relx=0.975, rely=0, relwidth=0.025, relheight=1)
        # Make the complete DAT-derived list (including explosion loading
        # fields) accessible with the mouse wheel as well as the scrollbar.
        inputs.bind("<MouseWheel>", lambda event: inputs.yview_scroll(-int(event.delta / 120) or -1, "units"))
        inputs.bind("<Button-4>", lambda event: inputs.yview_scroll(-1, "units"))
        inputs.bind("<Button-5>", lambda event: inputs.yview_scroll(1, "units"))
        tabs.select(inputs_frame)
        messages = ScrolledText(tabs, wrap="word", state="normal")
        messages.insert("1.0", _deck_summary(deck, path))
        messages.configure(state="disabled")
        tabs.add(messages, text="Messages")
        window.update_idletasks()
        tabs.select(0)

    def choose_output(self) -> None:
        path = filedialog.asksaveasfilename(title="Save numerical result", defaultextension=".npz", filetypes=(("NumPy result", "*.npz"),))
        if path:
            output = Path(path)
            if output.suffix.lower() != ".npz":
                output = output.with_suffix(".npz")
            # When a result is already displayed, Save as... should persist
            # that result immediately.  Before the first run it remains the
            # output destination for the next calculation.
            if self.result is not None:
                output.parent.mkdir(parents=True, exist_ok=True)
                np.savez(output, **self.result)
                self.result_path = output
                self.output_path.set(str(output))
                self.write_log(f"Saved result:\n{output.resolve()}")
                messagebox.showinfo("Result saved", f"Result saved to:\n{output.resolve()}", parent=self.root)
            else:
                self.output_path.set(str(output))

    def run_case(self) -> None:
        if self._running or self.case_filter.get() == SELECT_EXAMPLE or self.case_subcategory.get() == SELECT_SUBCATEGORY:
            return
        try:
            input_path = Path(self.input_path.get())
            if not input_path.is_file():
                raise ValueError("Select an existing .DAT input file first.")
            output_path = Path(self.output_path.get())
            if output_path.suffix.lower() != ".npz":
                output_path = output_path.with_suffix(".npz")
                self.output_path.set(str(output_path))
            deck = read_input_deck(input_path)
            if self.editor_path != input_path.resolve():
                case = next((case for case in self.examples if case.path.resolve() == input_path.resolve()), None)
                self._load_inputs(deck, input_path, case)
            original_deck = deck
            if self.editor_deck is not None:
                deck = self.editor_deck
            deck = edited_input_deck(
                deck,
                {key: value.get() for key, value in self.input_values.items()
                 if not self.input_entries[key].instate(["disabled"])},
            )
            self._modified_inputs = deck != original_deck
            self.active_deck = deck
            self.active_thesis_section = self.editor_case.thesis_section if self.editor_case else ""
        except Exception as error:
            messagebox.showerror("Input error", str(error), parent=self.root)
            return
        self._running = True
        self._solver_stage = "Preparing solver"
        self._progress_started = time.monotonic()
        self.run_stage.set("Preparing solver — 0 s")
        self.stop_button.configure(state="normal")
        self._update_run_button()
        self.progress.start(10)
        if self.editor_case:
            self.write_log(case_details(replace(self.editor_case, deck=deck)), clear=True)
        self.write_log("\nRunning solver. Please wait…")
        run_id = uuid.uuid4().hex
        self._run_id = run_id
        self._solver_output_path = output_path
        self._solver_temporary_path = output_path.with_name(
            f".{output_path.stem}.{run_id}.running.npz")
        self._solver_error_path = output_path.with_name(
            f".{output_path.stem}.{run_id}.error.txt")
        try:
            context = multiprocessing.get_context("spawn")
            self._solver_queue = context.Queue()
            self._solver_process = context.Process(
                target=run_gui_case,
                args=(deck, self._solver_temporary_path, self._solver_error_path,
                      self._solver_queue),
                daemon=True,
            )
            self._solver_process.start()
        except Exception:
            self._run_failed(traceback.format_exc())
            return
        self.root.after(100, lambda: self._poll_solver_stage(run_id))

    def _clear_solver_run(self) -> None:
        if self._solver_queue is not None:
            self._solver_queue.close()
            self._solver_queue = None
        for path in (self._solver_temporary_path, self._solver_error_path):
            if path is not None:
                path.unlink(missing_ok=True)
        self._solver_process = None
        self._solver_temporary_path = None
        self._solver_error_path = None
        self._solver_output_path = None
        self._run_id = None

    def _poll_solver_stage(self, run_id: str) -> None:
        if not self._running or self._run_id != run_id:
            return
        if self._solver_queue is not None:
            while True:
                try:
                    self._solver_stage = self._solver_queue.get_nowait()
                except queue.Empty:
                    break
        elapsed = int(time.monotonic() - self._progress_started)
        self.run_stage.set(f"{self._solver_stage} — {elapsed} s")
        process = self._solver_process
        if process is not None and process.is_alive():
            self.root.after(100, lambda: self._poll_solver_stage(run_id))
            return
        if process is not None:
            process.join(timeout=0)
        try:
            if process is None or process.exitcode != 0:
                error = self._solver_error_path
                details = (error.read_text(encoding="utf-8") if error and error.exists()
                           else f"Solver process exited with code {process.exitcode if process else 'unknown'}.")
                self._run_failed(details)
                return
            temporary = self._solver_temporary_path
            output = self._solver_output_path
            if temporary is None or output is None:
                raise RuntimeError("Solver output location was lost")
            loaded = load_result_archive(temporary)
            os.replace(temporary, output)
            self._clear_solver_run()
            self._run_finished(output, loaded)
        except Exception:
            self._run_failed(traceback.format_exc())

    def _cancel_run(self) -> None:
        if self._running:
            process = self._solver_process
            if process is not None and process.is_alive():
                process.terminate()
                process.join(timeout=0.25)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=0.25)
            self._run_cancelled()

    def _run_cancelled(self) -> None:
        self.progress.stop(); self._running = False; self._update_run_button()
        self.stop_button.configure(state="disabled")
        self._clear_solver_run()
        self.run_stage.set("Run cancelled")
        self.write_log("Run cancelled; no result was saved.")

    def _run_finished(self, path: Path, result: dict[str, np.ndarray]) -> None:
        self.progress.stop()
        self._running = False
        self.run_stage.set(f"Completed in {int(time.monotonic() - self._progress_started)} s")
        self.stop_button.configure(state="disabled")
        self._update_run_button()
        self._set_result(path, result)
        self.write_log(f"Completed successfully.\nSaved: {path.resolve()}")
        messagebox.showinfo("Calculation completed", f"Result saved to:\n{path.resolve()}", parent=self.root)

    def _run_failed(self, details: str) -> None:
        self.progress.stop()
        self._running = False
        self.stop_button.configure(state="disabled")
        self._clear_solver_run()
        self.run_stage.set("Calculation failed — see Messages")
        self._update_run_button()
        self.write_log("Calculation failed:\n" + details)
        messagebox.showerror("Calculation failed", "See the message panel for details.", parent=self.root)

    def open_result(self) -> None:
        path = filedialog.askopenfilename(title="Open solver result", filetypes=(("NumPy result", "*.npz"), ("All files", "*.*")))
        if not path:
            return
        try:
            result_path = Path(path).resolve()
            result = load_result_archive(result_path)
            self._restore_input_for_result(result_path)
            self._set_result(result_path, result)
            self.write_log(f"Opened result:\n{result_path}", clear=True)
        except Exception as error:
            messagebox.showerror("Cannot open result", str(error), parent=self.root)

    def _restore_input_for_result(self, result_path: Path) -> None:
        """Restore the matching DAT and input controls when opening an archive.

        Result archives intentionally contain numerical arrays only.  For the
        curated examples their stem identifies the source DAT (``FOO_result``
        corresponds to ``FOO.DAT``), so restore that case before plotting.
        If no matching deck can be found, retain the currently selected input.
        """
        stem = result_path.stem
        for suffix in ("_result", "-result"):
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        candidates = [case for case in self.examples if case.path.stem.lower() == stem.lower()]
        if not candidates:
            # Validation archives may add a descriptive suffix (for example
            # ``INFPL1_hoop_validated.npz``); use the longest matching DAT
            # stem so that the corresponding input controls are still restored.
            matching = [case for case in self.examples
                        if stem.lower().startswith(case.path.stem.lower() + "_")]
            if matching:
                candidates = [max(matching, key=lambda item: len(item.path.stem))]
        if not candidates:
            search_roots = [result_path.parent, example_directory()]
            for root in search_roots:
                for dat_path in root.glob(f"{stem}.[dD][aA][tT]"):
                    try:
                        deck = read_input_deck(dat_path)
                    except (OSError, ValueError):
                        continue
                    candidates = [ExampleCase(dat_path, deck, "", dat_path.name, "", "")]
                    break
                if candidates:
                    break
        if not candidates:
            return
        case = candidates[0]
        self.input_path.set(str(case.path))
        self.output_path.set(str(result_path))
        self.active_case = case
        self.active_deck = case.deck
        self.active_thesis_section = case.thesis_section
        self._load_inputs(case.deck, case.path, case)
        self.case_filter.set(case.family)
        self.case_subcategory.set(case.subcategory)
        self.case_choice.set(case.label)
        self._refresh_thesis_graph_choices()
        self._refresh_result_choices()

    def _set_result(self, path: Path, result: dict[str, np.ndarray]) -> None:
        self.result_path = path
        self.result = result
        self.node_box.configure(to=len(result["coordinates"]))
        if self.active_case is None and self.active_deck is not None:
            selected = tuple(int(n) for n in self.active_deck.node_plots
                             if 1 <= int(n) <= len(result["coordinates"]))
            selected_nodes = selected or tuple(range(1, len(result["coordinates"]) + 1))
            self.nodes_text.set(", ".join(map(str, selected_nodes)))
            self.node.set("")
        self.step_box.configure(to=max(1, len(result["time"])))
        self.node.set("1")
        self.report_button.configure(state="normal")
        self.csv_button.configure(state="normal")
        self._refresh_thesis_graph_choices()
        self.plot_result()

    def _history_nodes(self, node_count: int) -> list[int]:
        if self.thesis_history_nodes is not None:
            selected = [node for node in self.thesis_history_nodes if 1 <= node <= node_count]
            if selected:
                return selected
        # Manual controls take priority over the DAT or thesis preset.  This
        # lets All nodes, DAT nodes, and a single Node be switched at runtime.
        if self.nodes_text.get().strip():
            requested = tuple(int(token) for token in re.findall(r"\d+", self.nodes_text.get()))
            selected = [node for node in requested if 1 <= node <= node_count]
            if selected:
                return selected
        if self.node.get().strip().isdigit():
            node = int(self.node.get())
            if 1 <= node <= node_count:
                return [node]
        # Archived thesis cases may define a fixed NODE.PLO group.  A custom
        # DAT loaded by the user should always honor the manual Node control.
        requested = (self.active_deck.node_plots
                     if self.active_deck is not None else ())
        nodes = [int(node) for node in requested if 1 <= int(node) <= node_count]
        return nodes or [int(self.node.get())]

    def _caption(self, view: str, quantity: str, component: int) -> str:
        selected_graph = getattr(self, "thesis_graph", None)
        if selected_graph is not None and selected_graph.get():
            return selected_graph.get()
        if quantity == "traction" and view != "Deformed geometry":
            return "Applied boundary traction (MPa) — no separate PhD figure"
        figures = chapter3_figures(self.active_deck) if self.active_thesis_section else None
        if figures is None:
            return f"Figure unavailable in PhD - {quantity.title()} {view.lower()}"
        if view == "Time history" and quantity == "displacement":
            item = figures.x_displacement if component == 0 else figures.y_displacement
            if item is not None:
                return f"Figure {item[0]}  {item[1]}"
        if view == "Deformed geometry":
            item = figures.deformed_shape
            return f"Figure {item[0]}  {item[1]}"
        if quantity == "hoop stress":
            item = figures.boundary_profile if view == "Boundary profile" else figures.stress_history
            return f"Figure {item[0]}  {item[1]}"
        return f"{quantity.title()} {view.lower()} (no separate PhD figure)"

    def _node_inset_bounds(self) -> tuple[float, float, float, float]:
        """Use the outside margin when available, otherwise find clear space inside."""
        position = self.axes.get_position()
        if position.x1 < 0.76:
            return (0.40, 0.54, 0.22, 0.22)
        self.axes.get_xlim()
        self.axes.get_ylim()
        to_axes = self.axes.transAxes.inverted()
        samples = []
        for line in self.axes.lines:
            path = line.get_path().interpolated(32)
            points = to_axes.transform(line.get_transform().transform(path.vertices))
            samples.extend(points[np.isfinite(points).all(axis=1)])
        # Reserve a square inset; this keeps circular boundaries circular when
        # the complete graph is embedded in a PDF page.
        # Keep the geometry indicator large enough to read while preserving a
        # square box so circles are shown with the same scale in both axes.
        width = height = 0.38
        candidates = [(float(x), float(y), width, height)
                      for y in np.linspace(0.60, 0.03, 7)
                      for x in np.linspace(0.69, 0.03, 7)]
        if not samples:
            return candidates[0]
        points = np.asarray(samples)

        def clearance(bounds: tuple[float, float, float, float]) -> tuple[float, float]:
            x, y, w, h = bounds
            distances = np.maximum(np.maximum((x, y) - points, points - (x + w, y + h)), 0)
            squared = np.sum(distances ** 2, axis=1)
            return -float(np.count_nonzero(squared == 0)), float(np.min(squared))

        return max(candidates, key=clearance)

    def _add_node_inset(self, coordinates: np.ndarray) -> None:
        """Draw a compact geometry key with the boundary shape and node numbers."""
        if len(coordinates) < 3:
            return
        deck = self.active_deck
        # Keep the transient internal-pressure case stable in TkAgg; its
        # radial-pressure marker is shown in the report geometry page instead.
        if getattr(deck, "inside_pressure", 0.0) and self.view.get() == "Time history":
            return
        inset = self.axes.inset_axes(self._node_inset_bounds())
        inset.set_in_layout(False)
        self.node_inset = inset
        if self.nodes_text.get().strip():
            nodes = tuple(int(token) for token in re.findall(r"\d+", self.nodes_text.get()))
        elif self.node.get().strip().isdigit():
            nodes = (int(self.node.get()),)
        elif self.thesis_history_nodes:
            nodes = self.thesis_history_nodes
        elif len(coordinates) > 24:
            nodes = getattr(deck, "node_plots", ()) or tuple(np.linspace(1, len(coordinates), 8, dtype=int))
        else:
            nodes = range(1, len(coordinates) + 1)
        load = getattr(deck, "load", None)
        inside_pressure = getattr(deck, "inside_pressure", 0.0)
        draw_geometry_key(inset, coordinates, self.result["connectivity"], nodes,
                          load, inside_pressure=inside_pressure, fontsize=5)
        internal = self.result.get("internal_coordinates")
        if internal is None or not len(internal):
            internal = getattr(deck, "internal_points", ())
        if len(internal):
            pts = np.asarray(internal, dtype=float).reshape(-1, 2)
            labels = ("D", "F", "G", "K")
            inset.scatter(pts[:, 0], pts[:, 1], marker="x", s=18, color="#b3261e", zorder=6)
            for i, (x, y) in enumerate(pts):
                inset.annotate(labels[i] if i < len(labels) else f"I{i+1}", (x, y),
                               xytext=(3, 3), textcoords="offset points", fontsize=6,
                               color="#b3261e", fontweight="bold")
        if hasattr(load, "peak"):
            ramp = load
            xmin, xmax = float(np.min(coordinates[:, 0])), float(np.max(coordinates[:, 0]))
            ymin, ymax = float(np.min(coordinates[:, 1])), float(np.max(coordinates[:, 1]))
            height = max(0.08 * (ymax - ymin), 0.2)
            base = ymax
            half = ramp.plateau_length / 2.0
            xr = [ramp.center_x-half-ramp.ramp_length, ramp.center_x-half,
                  ramp.center_x+half, ramp.center_x+half+ramp.ramp_length]
            inset.plot(xr, [base, base+height, base+height, base], color="#b3261e", linewidth=1.2, zorder=5)
            inset.fill_between(xr, [base]*4, [base, base+height, base+height, base], color="#b3261e", alpha=.18, zorder=4)
            # Downward traction arrows, matching the reference schematic.
            for x_arrow in np.linspace(xr[0], xr[-1], 11):
                if x_arrow <= xr[1]:
                    h = height * (x_arrow-xr[0]) / ramp.ramp_length
                elif x_arrow <= xr[2]:
                    h = height
                else:
                    h = height * (xr[-1]-x_arrow) / ramp.ramp_length
                inset.annotate("", xy=(x_arrow, base), xytext=(x_arrow, base+h*.9),
                               arrowprops=dict(arrowstyle="-|>", color="#b3261e", lw=.7), zorder=6)
            inset.text((xr[0]+xr[1])/2, base+height*1.05, f"{ramp.ramp_length:g} m", fontsize=5, ha="center", color="#b3261e")
            inset.text((xr[1]+xr[2])/2, base+height*1.05, f"{ramp.plateau_length:g} m", fontsize=5, ha="center", color="#b3261e")
            inset.text((xr[2]+xr[3])/2, base+height*1.05, f"{ramp.ramp_length:g} m", fontsize=5, ha="center", color="#b3261e")
            inset.text(0.02, 0.98, f"ramp/plateau: {ramp.peak:g} MPa",
                       transform=inset.transAxes, va="top", fontsize=6, color="#b3261e")
        inset.set_aspect("equal", adjustable="box")

    def plot_result(self) -> None:
        self._refresh_result_choices()
        if self.result is None:
            return
        try:
            node_count = len(self.result["coordinates"])
            # Multi-node mode intentionally leaves Node empty. Use the first
            # listed node for the scalar controls and the complete list for
            # time-history drawing.
            if self.node.get().strip():
                node = int(self.node.get())
            else:
                listed = [int(token) for token in re.findall(r"\d+", self.nodes_text.get())]
                node = listed[0] if listed else 1
            if not 1 <= node <= node_count:
                raise ValueError(f"Node must be between 1 and {node_count}.")
            component = 0 if self.component.get() == "X" else 1
            column = 2 * (node - 1) + component
            quantity = self.quantity.get()
            step = min(max(int(self.step.get()), 1), len(self.result["time"])) - 1
            if getattr(self, "node_inset", None) is not None:
                self.node_inset.remove()
                self.node_inset = None
            self.axes.clear()
            # ``set_aspect('equal')`` is required for the deformed-geometry
            # view, but it must not leak into time-history/profile plots.
            # Otherwise a previous geometry plot can force the report graph
            # into a tall, very narrow frame and make the curves appear
            # vertically stretched.
            self.axes.set_aspect("auto")
            if quantity in STRESS_COMPONENTS and self.view.get() != "Deformed geometry":
                if "boundary_stress" not in self.result:
                    self.axes.text(0.5, 0.5, "This saved result does not contain boundary stresses.\nRun the case again to calculate them.",
                                   transform=self.axes.transAxes, ha="center", va="center")
                    self.figure.texts.clear()
                    self.canvas.draw()
                    return
                stress_values = np.asarray(self.result["boundary_stress"])
                if stress_values.ndim == 2:
                    expected_columns = 3 * node_count
                    if stress_values.shape != (len(self.result["time"]), expected_columns):
                        raise ValueError(f"boundary_stress has incompatible shape {stress_values.shape}")
                    stress_values = stress_values.reshape(len(self.result["time"]), node_count, 3)
                quantity_values = stress_values[:, :, STRESS_COMPONENTS[quantity]]
            else:
                quantity_values = self.result["displacement" if self.view.get() == "Deformed geometry" else quantity][:, component::2]
            if self.view.get() == "Time history":
                history_nodes = self._history_nodes(node_count)
                line_styles = ("-", (0, (2, 2)), (0, (7, 3, 1, 3)), (0, (7, 3)), (0, (1, 3)), (0, (10, 3)))
                line_colors = ("#1565c0", "#d32f2f", "#2e7d32", "#7b1fa2", "#ef6c00", "#00838f")
                step_numbers = np.arange(1, len(self.result["time"]) + 1)
                multiplier = 1000.0 if quantity == "displacement" else 1.0
                # Keep curve identifiers visible for the seven-node 5.25 case.
                show_node_labels = len(history_nodes) <= 12
                numerical_zero = None
                if quantity in {"radial stress", "shear stress"}:
                    selected_stress = quantity_values[:, np.asarray(history_nodes) - 1]
                    peak_stress = float(np.nanmax(np.abs(selected_stress)))
                    if peak_stress < 1e-10:
                        numerical_zero = peak_stress
                        show_node_labels = False
                for series_index, history_node in enumerate(history_nodes):
                    values = multiplier * quantity_values[:, history_node - 1]
                    if numerical_zero is not None:
                        values = np.zeros_like(values)
                    color = line_colors[series_index % len(line_colors)]
                    self.axes.plot(step_numbers, values, color=color, linestyle=line_styles[series_index % len(line_styles)], linewidth=1.7,
                                   label=f"Node {history_node}")
                    # Spread curve labels along the later part of the history
                    # so neighbouring series do not form an unreadable stack.
                    label_fraction = min(.88, .56 + .065 * series_index)
                    label_index = min(len(values) - 1, max(1, int(label_fraction * len(values))))
                    label_offset = (4, 5 + 3 * (series_index % 3))
                    if show_node_labels:
                        self.axes.annotate(str(history_node), (step_numbers[label_index], values[label_index]),
                            xytext=label_offset, textcoords="offset points", fontsize=11, color=color,
                            fontweight="bold", bbox=dict(facecolor="white", edgecolor="none",
                            alpha=0.75, pad=0.8), zorder=8)
                if history_nodes and numerical_zero is None:
                    self.axes.legend(loc="lower right", fontsize=7, ncol=min(4, len(history_nodes)),
                                     framealpha=0.85)
                # Also plot internal observation-point histories when the
                # DAT file contains INTERNAL.POINTS. These are not boundary
                # nodes, so they must be read from the dedicated result array.
                internal = internal_displacement_histories(self.result, quantity)
                if internal is not None:
                    labels = ("D", "F", "G", "K")
                    for point_index in range(internal.shape[1]):
                        values = 1000.0 * internal[:, point_index, component]
                        color = ("#8e24aa", "#ef6c00", "#00897b", "#6d4c41")[point_index % 4]
                        if point_index < len(labels) and "internal_coordinates" in self.result:
                            px, py = np.asarray(self.result["internal_coordinates"])[point_index]
                            label = f"{labels[point_index]} ({px:g},{py:g})"
                        else:
                            label = labels[point_index] if point_index < len(labels) else f"I{point_index + 1}"
                        self.axes.plot(step_numbers, values, linestyle=(0, (5, 2)),
                                       linewidth=1.5, color=color, label=label)
                        index = min(len(values) - 1, max(1, int(.82 * len(values))))
                        self.axes.annotate(label, (step_numbers[index], values[index]),
                                           xytext=(4, 4), textcoords="offset points",
                                           fontsize=10, color=color, fontweight="bold")
                self.axes.axhline(0, color="black", linewidth=.8)
                if numerical_zero is not None:
                    self.axes.set_ylim(-1, 1)
                    self.axes.text(0.5, 0.92,
                                   f"Numerical zero at selected nodes (max |stress| = {numerical_zero:.2e} MPa)",
                                   transform=self.axes.transAxes, ha="center", va="top", fontsize=8)
                self.axes.spines["top"].set_visible(False)
                self.axes.spines["right"].set_visible(False)
                self.axes.set_xlim(0, len(step_numbers))
                # Keep the two-line step/time labels readable at every zoom level.
                tick_stride = max(1, int(np.ceil(len(step_numbers) / 6)))
                tick_steps = np.arange(0, len(step_numbers) + 1, tick_stride)
                if tick_steps[-1] != len(step_numbers):
                    tick_steps = np.r_[tick_steps, len(step_numbers)]
                total_time = self.active_deck.total_time if self.active_deck is not None else float(self.result["time"][-1])
                tick_seconds = tick_steps * total_time / len(step_numbers)
                self.axes.set_xticks(tick_steps)
                self.axes.set_xticklabels([f"{step:g}\n{seconds:.5g}" for step, seconds in zip(tick_steps, tick_seconds, strict=True)], fontsize=7)
                self.axes.set_xlabel("time step\nelapsed time (sec) (X)", fontweight="bold")
                component_name = (r"$\mathbf{u}_{\mathrm{x}}$" if component == 0 else r"$\mathbf{u}_{\mathrm{y}}$") if quantity == "displacement" else (r"$\mathbf{t}_{\mathrm{x}}$" if component == 0 else r"$\mathbf{t}_{\mathrm{y}}$")
                unit = " (mm)" if quantity == "displacement" else " (MPa)"
                self.axes.set_ylabel(f"{component_name} {quantity}{unit} (Y)", fontweight="bold")
                if quantity in STRESS_COMPONENTS:
                    symbol = {"hoop stress": r"$\mathbf{\sigma}_{22}$", "radial stress": r"$\mathbf{\sigma}_{11}$", "shear stress": r"$\mathbf{\sigma}_{12}$"}[quantity]
                    self.axes.set_ylabel(f"{symbol} (MPa) (Y)", fontweight="bold")
                self._add_node_inset(self.result["coordinates"])
            elif self.view.get() == "Boundary profile":
                if quantity == "hoop stress":
                    draw_stress_profile(self.axes, self.result["coordinates"],
                                        self.result["connectivity"], quantity_values[step])
                else:
                    nodes = np.arange(1, node_count + 1)
                    self.axes.plot(nodes, quantity_values[step], color="#9b3c2d", marker="o", markersize=3)
                    # Boundary profiles are defined over every boundary node,
                    # but the selected node set should remain visible and
                    # consistent with the time-history controls.  Highlight
                    # each selected node rather than letting the scalar Node
                    # field make the plot look like a single-node result.
                    selected_profile_nodes = self._history_nodes(node_count)
                    profile_colors = ("#1565c0", "#d32f2f", "#2e7d32", "#7b1fa2", "#ef6c00", "#00838f")
                    for series_index, selected_node in enumerate(selected_profile_nodes):
                        value = quantity_values[step, selected_node - 1]
                        color = profile_colors[series_index % len(profile_colors)]
                        self.axes.plot(selected_node, value, marker="o", markersize=5,
                                       color=color, markeredgecolor="white", markeredgewidth=.7,
                                       linestyle="None", zorder=4)
                        self.axes.annotate(str(selected_node), (selected_node, value),
                                           xytext=(3, 4), textcoords="offset points",
                                           fontsize=8, color=color, fontweight="bold",
                                           zorder=5)
                    max_node = int(np.nanargmax(np.abs(quantity_values[step])) + 1)
                    max_value = float(quantity_values[step, max_node - 1])
                    self.axes.plot(max_node, max_value, marker="*", markersize=11,
                                   color="#111827", markeredgecolor="white", zorder=6)
                    self.axes.annotate(f"max: node {max_node}", (max_node, max_value),
                                       xytext=(6, -12), textcoords="offset points",
                                       fontsize=8, fontweight="bold", color="#111827", zorder=7)
                    self.axes.axhline(0, color="black", linewidth=.7)
                    self.axes.set_xlabel("Boundary node (X)")
                    self.axes.set_ylabel(f"{quantity.title()} ({'m' if quantity == 'displacement' else 'MPa'}) (Y)")
                self.axes.set_title(f"{quantity.title()} around boundary - t={self.result['time'][step]:g} s", fontsize=9)
                self._add_node_inset(self.result["coordinates"])
            else:
                coordinates = self.result["coordinates"]
                displacement = self.result["displacement"][step].reshape(-1, 2)
                span = max(float(np.ptp(coordinates[:, 0])), float(np.ptp(coordinates[:, 1])), 1.0)
                peak = max(float(np.max(np.linalg.norm(displacement, axis=1))), 1e-15)
                scale = 0.08 * span / peak
                deformed = coordinates + scale * displacement
                draw_boundary(self.axes, coordinates, self.result["connectivity"], "--", color="#777777", label="Original")
                draw_boundary(self.axes, deformed, self.result["connectivity"], color="#1f5a91", marker="o", markersize=3, label=f"Deformed (enlargement: x{scale:.0f})")
                if self.active_deck is not None:
                    draw_explosion(self.axes, self.active_deck.load)
                self.axes.set_aspect("equal", adjustable="box")
                self.axes.set_title(f"Deformed boundary - t={self.result['time'][step]:g} s")
                self.axes.set_xlabel("Horizontal coordinate (X)")
                self.axes.set_ylabel("Vertical coordinate (Y)")
                self.axes.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
                                 borderaxespad=0.0, frameon=True)
                self._add_node_inset(coordinates)
            if self.view.get() == "Time history":
                self.axes.grid(False)
            else:
                self.axes.grid(True, alpha=0.25)
            self.axes.xaxis.label.set_fontsize(8)
            self.axes.yaxis.label.set_fontsize(8)
            self.axes.tick_params(axis="both", labelsize=7)
            self.axes.xaxis.label.set_alpha(0.72)
            self.axes.yaxis.label.set_alpha(0.72)
            for tick in (*self.axes.get_xticklabels(), *self.axes.get_yticklabels()):
                tick.set_alpha(0.68)
            self.figure.texts.clear()
            if getattr(self, "_report_rendering", False):
                current_view = self.view.get()
                time_value = float(self.result["time"][step])
                if current_view == "Time history":
                    if quantity == "displacement":
                        component_name = "X" if component == 0 else "Y"
                        caption = f"{component_name} displacements versus time steps"
                    elif quantity == "traction":
                        component_name = "X" if component == 0 else "Y"
                        caption = f"{component_name} tractions versus time steps"
                    else:
                        caption = f"{quantity.title()} versus time steps"
                elif current_view == "Boundary profile":
                    profile_values = quantity_values[step]
                    max_node = int(np.nanargmax(np.abs(profile_values)) + 1)
                    if quantity == "hoop stress":
                        label = "Tangential (hoop) stresses"
                    else:
                        label = quantity.title()
                    max_value = float(profile_values[max_node - 1])
                    unit = "MPa" if quantity in STRESS_COMPONENTS or quantity == "traction" else "m"
                    caption = f"{label} around boundary\nMax value: {max_value:.4g} {unit} at step {int(self.step.get())}, node {max_node}"
                else:
                    disp_step_values = np.linalg.norm(self.result["displacement"][step].reshape(-1, 2), axis=1)
                    max_node = int(np.argmax(disp_step_values) + 1)
                    caption = f"Displacement\nMax value: {float(disp_step_values[max_node - 1]):.4g} m at step {int(self.step.get())}, node {max_node}"
            elif self.view.get() == "Time history" and quantity == "traction":
                component_name = "X" if component == 0 else "Y"
                caption = f"{component_name} tractions versus time steps"
            elif self.view.get() == "Boundary profile" and quantity == "traction":
                caption = f"{quantity.title()} around boundary - t={float(self.result['time'][step]):g} s"
            else:
                caption = self._caption(self.view.get(), quantity, component)
            internal_coordinates = self.result.get("internal_coordinates")
            if (self.view.get() == "Time history" and quantity != "displacement"
                    and internal_coordinates is not None and len(internal_coordinates)):
                labels = ("D", "F", "G", "K")
                point_labels = [
                    f"{labels[index] if index < len(labels) else f'I{index + 1}'} ({x:g},{y:g})"
                    for index, (x, y) in enumerate(np.asarray(internal_coordinates).reshape(-1, 2))
                ]
                caption += (
                    "\nInternal points " + ", ".join(point_labels)
                    + ": displacement histories are shown on displacement plots; "
                    + "internal-point traction/stress histories are not calculated."
                )
            if self._modified_inputs:
                caption += " (modified inputs)"
            self.figure.text(0.5, 0.012, caption, ha="center", fontsize=10)
            self.figure.tight_layout(rect=(0, 0.06, 1, 1))
            self.canvas.draw()
        except Exception as error:
            messagebox.showerror("Plot error", str(error), parent=self.root)

    def export_csv(self) -> None:
        if self.result is None:
            return
        default = (self.result_path.stem if self.result_path else "result") + ".csv"
        path = filedialog.asksaveasfilename(title="Export readable results", initialfile=default, defaultextension=".csv", filetypes=(("CSV table", "*.csv"),))
        if not path:
            return
        try:
            export_result_csv(self.result, path)
            messagebox.showinfo("CSV exported", f"Saved:\n{Path(path).resolve()}", parent=self.root)
        except Exception as error:
            messagebox.showerror("CSV export failed", str(error), parent=self.root)

    def export_report(self) -> None:
        """Write a professional, browser-printable HTML report."""
        if self.result is None:
            return
        default = (self.result_path.stem if self.result_path else "result") + "_report.pdf"
        path = filedialog.asksaveasfilename(title="Export professional PDF report",
                                            initialdir=str(Path(__file__).resolve().parents[2] / "examples/results/pdf"),
                                            initialfile=default,
                                            defaultextension=".pdf", filetypes=(("PDF report", "*.pdf"),))
        if not path:
            return
        self.progress.configure(mode="indeterminate")
        self.progress.start(10)
        self.write_log("Generating PDF report…", clear=False)
        self.root.update()
        try:
            result = self.result
            title = f"BEM - {self.active_deck.title}" if self.active_deck else "BEM - Boundary Element Underground result"
            project_root = Path(__file__).resolve().parents[2]
            def display_path(value: str | Path) -> str:
                if not value:
                    return "Not recorded"
                try:
                    return ".\\" + str(Path(value).resolve().relative_to(project_root))
                except (ValueError, OSError):
                    return str(value)
            input_file = display_path(self.input_path.get())
            result_file = display_path(self.result_path) if self.result_path else "Not recorded"
            coords = result["coordinates"]; disp = result["displacement"]; traction = result["traction"]
            peak_u = float(np.max(np.linalg.norm(disp.reshape(-1, 2), axis=1)))
            peak_t = float(np.max(np.abs(traction)))
            input_details = _format_input_details(asdict(self.active_deck)) if self.active_deck else "No parsed input details available."
            result_ranges = [
                f"Displacement range: {float(np.min(disp)):.6g} to {float(np.max(disp)):.6g} m",
                f"Traction range: {float(np.min(traction)):.6g} to {float(np.max(traction)):.6g} MPa",
                f"Peak displacement: {peak_u:.6g} m",
                f"Peak traction: {peak_t:.6g} MPa",
            ]
            if self.active_deck and isinstance(self.active_deck.load, ExplosionPoint):
                # Keep the computed boundary traction separate from the
                # analytical incident stress specified by the explosion law.
                traction_vectors = np.asarray(traction).reshape(len(result["time"]), len(coords), 2)
                node1_traction = traction_vectors[:, 0, :]
                direct_compressive = float(np.min(node1_traction))
                result_ranges.append(
                    f"Computed boundary traction at node 1: {direct_compressive:.6g} MPa"
                )
                try:
                    if not isinstance(self.active_deck.load, ExplosionPoint):
                        raise ValueError("analytical explosion stress applies only to point-explosion loading")
                    from .loads import _explosion_pressure
                    load = self.active_deck.load
                    material = self.active_deck.material
                    shear = material.young / (2.0 * (1.0 + material.poisson))
                    lame = material.young * material.poisson / ((1.0 + material.poisson) * (1.0 - 2.0 * material.poisson))
                    p_speed = np.sqrt((lame + 2.0 * shear) / material.density)
                    distance = float(np.linalg.norm(np.asarray(coords[0]) - np.array([load.x, load.y])))
                    incident = float(np.real(_explosion_pressure(
                        np.array([distance]), p_speed * load.p_amplitude, load, material.density
                    )[0]))
                    result_ranges.append(
                        f"Analytical incident compressive stress at node 1: {-incident:.6g} MPa"
                    )
                except Exception:
                    result_ranges.append("Analytical incident compressive stress at node 1: unavailable")
            if "boundary_stress" in result:
                stress = result["boundary_stress"]
                result_ranges.append(f"Boundary stress range: {float(np.min(stress)):.6g} to {float(np.max(stress)):.6g} MPa")
            rows = "".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k,v in [
                ("Input file", input_file), ("Result file", result_file),
                ("Domain steps", len(result["time"])), ("Boundary nodes", len(coords)),
                ("Boundary elements", len(result["connectivity"])),])
            from matplotlib.backends.backend_pdf import PdfPages
            from matplotlib.figure import Figure
            from matplotlib.lines import Line2D
            generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
            page_number = 0
            summary = [f"{k}: {v}" for k,v in [("Input file", input_file), ("Result file", result_file),
                ("Domain steps", len(result["time"])), ("Boundary nodes", len(coords)),
                ("Boundary elements", len(result["connectivity"]))]]
            if self.active_deck:
                description_lines = _deck_summary(self.active_deck, Path(self.input_path.get())).splitlines()
                description_lines = [line for line in description_lines if not line.startswith("Input:")]
                if getattr(self.active_deck, "inside_pressure", 0.0):
                    description_lines += [
                        "Internal explosion interpretation:",
                        ("The explosion is located at the centre of the cylindrical cavity "
                         "(0, 0) m. An internal radial pressure of "
                         f"{self.active_deck.inside_pressure:g} MPa is applied suddenly and "
                         "maintained constant."),
                    ]
                if isinstance(self.active_deck.load, ExplosionPoint):
                    load = self.active_deck.load
                    description_lines += [
                        "Explosion-wave interpretation:",
                        (f"Wave speed={load.wave_velocity:g} m/s; "
                         f"source=({load.x:g}, {load.y:g}) m; formula F={load.formula}; "
                         f"B={load.coefficient:g}; N={load.exponent:g}; decay={int(load.decay)}; "
                         f"P multiplier={load.p_amplitude:g}; SV multiplier={load.sv_amplitude:g}."),
                        ("Analytical incident stress is calculated before solving the boundary system. "
                         "The boundary traction is reported separately."),
                        ("For a cubic-root explosion, the calculation is: "
                         "P0 = B * F * rho * cP * (R / W^(1/3))^(-N), and sigma_incident = -P0."),
                        ("PhD GP-bomb reference at node 1: R=5.5 m, W=188.7 kg explosive, "
                         "B=160, F=1, N=2.75, with P=1 and SV=0. Substitution in the published "
                         "calculation gives sigma_incident approximately -33.7 MPa."),
                        ("Numerical substitution: W^(1/3)=188.7^(1/3)=5.737; "
                         "R/W^(1/3)=5.5/5.737=0.9585; "
                         "(R/W^(1/3))^(-2.75)=0.9585^(-2.75)=1.122. "
                         "Using the published wave-factor normalization, "
                         "P0=30.02*1.122=33.7 MPa, so sigma_incident=-P0=-33.7 MPa."),
                    ]
                if getattr(self.active_deck, "internal_points", ()):
                    description_lines += ["Internal displacement-history points:"]
                    for index, (x_point, y_point) in enumerate(self.active_deck.internal_points, 1):
                        description_lines.append(f"Point {index}: ({x_point:g}, {y_point:g}) m")
            else:
                description_lines = ["Input deck details are not available."]
            with PdfPages(path) as pdf:
                def save_page(figure):
                    nonlocal page_number
                    # Keep Tk responsive while a multi-page PDF is assembled.
                    self.root.update()
                    page_number += 1
                    figure.text(.08, .025, f"Generated {generated} | Boundary Element Underground", fontsize=8, color="#52606d")
                    figure.text(.92, .025, f"Page {page_number}", fontsize=8, color="#52606d", ha="right")
                    pdf.savefig(figure)
                    self.root.update()
                # Dedicated report cover: the title is always at the top of
                # page 1, independent of the length of the input summary.
                cover = Figure(figsize=(8.27, 11.69)); cover_ax = cover.add_subplot(111); cover_ax.axis("off")
                # Keep long DAT titles inside the A4 text margins.  Character
                # wrapping is paired with a size cap so long explosion titles
                # cannot run past the right edge of the cover.
                title_lines = textwrap.wrap(title, width=42, break_long_words=False, break_on_hyphens=False) or [""]
                title_size = 16 if len(title_lines) <= 2 else 13
                cover_ax.text(.08, .96, "\n".join(title_lines),
                              fontsize=title_size, weight="bold", color="#1f4e79",
                              va="top", ha="left", wrap=True)
                cover_ax.text(.08, .86, "Boundary Element Underground numerical report",
                              fontsize=11, color="#52606d", va="top")
                save_page(cover)
                fig = Figure(figsize=(8.27, 11.69)); ax = fig.add_subplot(111); ax.axis("off")
                y = .985
                y = .95
                ax.text(.08, y, "Data Input", fontsize=13, weight="bold", color="#1f4e79", va="top"); y -= .045
                for line in summary: ax.text(.1, y, line, fontsize=10, va="top"); y -= .028
                # Start the input-description section on its own page so the
                # section heading stays directly above its geometry figure.
                geometry_page = Figure(figsize=(8.27, 11.69)); geometry_ax = geometry_page.add_subplot(111)
                geometry_page.text(.08, .95, "Data Input", fontsize=13, weight="bold", color="#1f4e79", va="top")
                geometry_ax.set_title("Input geometry and loading", fontsize=15, weight="bold", color="#1f4e79", pad=18)
                geometry_nodes = tuple(range(1, len(coords) + 1))
                draw_geometry_key(geometry_ax, coords, result["connectivity"], geometry_nodes,
                                  self.active_deck.load if self.active_deck else None,
                                  inside_pressure=self.active_deck.inside_pressure if self.active_deck else 0.0, fontsize=7)
                if self.active_deck is not None:
                    draw_liner_materials(geometry_ax, self.active_deck, fontsize=8)
                low = np.min(coords, axis=0); high = np.max(coords, axis=0)
                span = max(float(np.max(high - low)), 1e-9); y_dim = low[1] - 0.16 * span; x_dim = low[0] - 0.16 * span
                geometry_ax.annotate("", xy=(high[0], y_dim), xytext=(low[0], y_dim), arrowprops=dict(arrowstyle="<->", color="#1f4e79", lw=0.9))
                geometry_ax.text(float((low[0] + high[0]) / 2), y_dim, f"X extent = {high[0] - low[0]:g} m", ha="center", va="bottom", fontsize=8, color="#1f4e79")
                geometry_ax.annotate("", xy=(x_dim, high[1]), xytext=(x_dim, low[1]), arrowprops=dict(arrowstyle="<->", color="#1f4e79", lw=0.9))
                geometry_ax.text(x_dim, float((low[1] + high[1]) / 2), f"Y extent = {high[1] - low[1]:g} m", ha="right", va="center", rotation=90, fontsize=8, color="#1f4e79")
                geometry_ax.set_xlabel("X (m)"); geometry_ax.set_ylabel("Y (m)"); geometry_ax.grid(True, alpha=0.2); geometry_ax.set_aspect("equal", adjustable="box")
                geometry_page.text(.08, .04, "Boundary dimensions, all node labels, and loading are taken directly from the DAT input.", fontsize=8, color="#52606d")
                save_page(geometry_page)
                for line in description_lines:
                    for wrapped in textwrap.wrap(line, width=92) or [""]:
                        ax.text(.1, y, wrapped, fontsize=9, va="top"); y -= .025
                y -= .03; ax.text(.08, y, "Parsed input data", fontsize=13, weight="bold", color="#1f4e79", va="top"); y -= .045
                for line in input_details.splitlines():
                    for wrapped in textwrap.wrap(line, width=92) or [""]:
                        ax.text(.1, y, wrapped, fontsize=8, va="top"); y -= .021
                        if y < .12:
                            save_page(fig); fig.clear()
                            fig = Figure(figsize=(8.27, 11.69)); ax = fig.add_subplot(111); ax.axis("off"); y = .95
                            ax.text(.08, y, "Parsed input data continued", fontsize=13, weight="bold", color="#1f4e79", va="top"); y -= .05
                y -= .03; ax.text(.08, y, "Numerical result ranges", fontsize=13, weight="bold", color="#1f4e79", va="top"); y -= .045
                for line in result_ranges: ax.text(.1, y, line, fontsize=9, va="top"); y -= .028
                y -= .03; ax.text(.08, y, "Interpretation", fontsize=13, weight="bold", color="#1f4e79", va="top"); y -= .045
                ax.text(.1, y, "The reported values summarize the saved numerical result. Use the GUI plots and exported CSV for node-by-node and time-step inspection.", fontsize=9, va="top", wrap=True)
                save_page(fig); fig.clear()
                saved_view, saved_quantity, saved_component = self.view.get(), self.quantity.get(), self.component.get()
                self._report_rendering = True
                boundary_stress = result.get("boundary_stress")
                stress_steps = {}
                if boundary_stress is not None:
                    boundary_stress = np.asarray(boundary_stress)
                    if boundary_stress.ndim == 2:
                        boundary_stress = boundary_stress.reshape(len(result["time"]), len(coords), 3)
                    for stress_name, stress_index in STRESS_COMPONENTS.items():
                        stress_steps[stress_name] = int(np.unravel_index(np.nanargmax(np.abs(boundary_stress[:, :, stress_index])), boundary_stress[:, :, stress_index].shape)[0] + 1)
                displacement_norm = np.linalg.norm(result["displacement"].reshape(len(result["time"]), len(result["coordinates"]), 2), axis=2)
                displacement_step = int(np.unravel_index(np.nanargmax(displacement_norm), displacement_norm.shape)[0] + 1)
                quantities = ("displacement", "traction")
                if "boundary_stress" in result:
                    quantities += ("hoop stress", "radial stress", "shear stress")
                for view_name in RESULT_VIEWS:
                    for quantity_name in quantities:
                        # The thesis report uses the published 14th-step
                        # boundary figures; omit generic t=0 displacement and
                        # traction snapshots from the report.
                        if view_name == "Boundary profile" and quantity_name != "hoop stress":
                            continue
                        if view_name == "Deformed geometry" and quantity_name != "displacement":
                            continue
                        components = ("X", "Y") if quantity_name in {"displacement", "traction"} and view_name == "Time history" else (saved_component,)
                        for component_name in components:
                            if view_name in {"Boundary profile", "Deformed geometry"}:
                                profile_step = stress_steps.get(quantity_name, displacement_step)
                                # Show deformation at the hoop-stress peak so
                                # the two geometry pages describe one instant.
                                geometry_step = stress_steps.get("hoop stress", displacement_step)
                                self.step.set(str(profile_step if view_name == "Boundary profile" else geometry_step))
                            self.view.set(view_name); self.quantity.set(quantity_name); self.component.set(component_name)
                            self._refresh_result_choices()
                            self.plot_result()
                            self.root.update()
                            # Report pages should remain crisp even though the
                            # interactive canvas uses softened axis annotations.
                            self.axes.xaxis.label.set_alpha(1.0)
                            self.axes.yaxis.label.set_alpha(1.0)
                            for tick in (*self.axes.get_xticklabels(), *self.axes.get_yticklabels()):
                                tick.set_alpha(1.0)
                            page_number += 1
                            graph_buffer = io.BytesIO()
                            # Use a high-resolution intermediate so axis labels,
                            # node numbers, and the inset legend stay readable in
                            # the PDF after zooming or printing.
                            self.figure.savefig(graph_buffer, format="png", dpi=720, bbox_inches="tight", pad_inches=0.02)
                            graph_buffer.seek(0)
                            graph_page = Figure(figsize=(8.27, 11.69)); graph_ax = graph_page.add_subplot(111); graph_ax.axis("off")
                            from matplotlib.image import imread
                            graph_image = imread(graph_buffer)
                            # Use one fixed report frame for every graph. The
                            # interactive plots have different natural aspect
                            # ratios (especially equal-aspect hoop profiles),
                            # which otherwise makes their apparent zoom and
                            # typography vary from page to page.
                            frame_left, frame_right = 0.08, 0.92
                            frame_bottom, frame_top = 0.10, 0.88
                            frame_width, frame_height = frame_right - frame_left, frame_top - frame_bottom
                            image_ratio = graph_image.shape[1] / graph_image.shape[0]
                            frame_ratio = frame_width / frame_height
                            if image_ratio >= frame_ratio:
                                image_width = frame_width
                                image_height = image_width / image_ratio
                            else:
                                image_height = frame_height
                                image_width = image_height * image_ratio
                            image_left = (frame_left + frame_right - image_width) / 2
                            image_bottom = (frame_bottom + frame_top - image_height) / 2
                            graph_ax.imshow(graph_image, interpolation="none", aspect="equal",
                                            extent=(image_left, image_left + image_width,
                                                    image_bottom, image_bottom + image_height))
                            graph_page.text(.08, .935, f"Displayed result graph — {view_name}, {quantity_name}", fontsize=12, weight="bold", color="#1f4e79", va="top")
                            graph_page.text(.08, .985, title, fontsize=7.5, color="#697784", va="top")
                            graph_page.lines.append(Line2D([.08, .92], [.965, .965], transform=graph_page.transFigure,
                                                               color="#c8d0d7", linewidth=0.7))
                            graph_page.text(.08, .025, f"Generated {generated} | Boundary Element Underground", fontsize=8, color="#52606d")
                            graph_page.text(.92, .025, f"Page {page_number}", fontsize=8, color="#52606d", ha="right")
                            pdf.savefig(graph_page)
                self._report_rendering = False
                self.view.set(saved_view); self.quantity.set(saved_quantity); self.component.set(saved_component); self._refresh_result_choices(); self.plot_result()
                # Numerical appendix: compact tables matching the CSV data,
                # split by node groups so every value remains printable.
                table_specs = [("Displacement", disp, "m", 2), ("Traction", traction, "MPa", 2)]
                if "boundary_stress" in result:
                    table_specs.append(("Boundary stress", result["boundary_stress"].reshape(len(result["time"]), -1), "MPa", 3))
                for table_title, values, unit, components_per_node in table_specs:
                    node_count = len(coords)
                    # Two nodes per page keeps scientific notation and headers
                    # readable on A4; coordinates are listed beneath the title.
                    for first_node in range(0, node_count, 2):
                        last_node = min(node_count, first_node + 2)
                        cols = ["Time (s)"]
                        for node_index in range(first_node, last_node):
                            x, y_coord = coords[node_index]
                            components = ("X", "Y") if components_per_node == 2 else ("Radial", "Hoop", "Shear")
                            cols.extend(f"Node {node_index + 1}\n{component}\n({x:g}, {y_coord:g})" for component in components)
                        body = [[f"{float(result['time'][row]):.6g}"] +
                                [f"{float(value):.6g}" for value in values[row, first_node * components_per_node:last_node * components_per_node]]
                               for row in range(len(result["time"]))]
                        table_page = Figure(figsize=(8.27, 11.69)); table_ax = table_page.add_subplot(111); table_ax.axis("off")
                        table_ax.text(.08, .935, f"Numerical data — {table_title} ({unit})", fontsize=13, weight="bold", color="#1f4e79", va="top")
                        table_ax.text(.08, .985, title, fontsize=7.5, color="#697784", va="top")
                        table_page.lines.append(Line2D([.08, .92], [.965, .965], transform=table_page.transFigure,
                                                           color="#c8d0d7", linewidth=0.7))
                        table = table_ax.table(cellText=body, colLabels=cols, loc="upper center", cellLoc="right", colLoc="center", bbox=(.04, .08, .92, .80))
                        table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1, 1.3)
                        for cell in table.get_celld().values():
                            cell.set_edgecolor("#b8c4ce")
                        for column in range(len(cols)):
                            table[0, column].get_text().set_wrap(True)
                            table[0, column].set_height(0.075)
                        save_page(table_page)
                        self.root.update()

            messagebox.showinfo("Report exported", f"Saved:\n{Path(path).resolve()}", parent=self.root)
            webbrowser.open(Path(path).resolve().as_uri())
        except Exception as error:
            messagebox.showerror("Report export failed", str(error), parent=self.root)
        finally:
            self.progress.stop()
            self.root.update_idletasks()


def main() -> None:
    root = Tk()
    BemApplication(root)
    root.mainloop()


if __name__ == "__main__":
    main()









