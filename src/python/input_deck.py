"""Reader for the keyword-based PREBEM input deck.

This first reader covers the modules exercised by the supplied INFPL1 and
TINFPL1 reference cases. Unknown modules fail explicitly instead of being
silently misinterpreted.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CompoundSegment:
    """One straight or circular/elliptical ``GEOCOMB`` segment."""

    shape: int
    elements: int
    geometry: tuple[float, ...]


@dataclass(frozen=True)
class Section:
    shape: int
    region: int
    x_elements: int
    geometry: tuple[float, ...]
    y_elements: int | None = None
    segments: tuple[CompoundSegment, ...] = ()


@dataclass(frozen=True)
class MaterialInput:
    young: float
    poisson: float
    density: float
    viscosity: float = 0.0
    viscosity_model: int = 1


@dataclass(frozen=True)
class ZoneInput:
    """Material and boundary assignment for one MULTI-ZONE region."""

    number: int
    outer_boundary: int
    inner_boundary: int
    material: MaterialInput
    boundaries: tuple[int, ...] = ()


@dataclass(frozen=True)
class PlaneWave:
    longitudinal_stress: float
    transverse_shear: float
    angle_degrees: float
    vertical_stress: float = 0.0
    frequency: float = 0.0
    xdistance: float | None = None
    ydistance: float | None = None

    @property
    def wavefront_point(self) -> tuple[float, float] | None:
        if self.xdistance is None or self.ydistance is None:
            return None
        return (self.xdistance, self.ydistance)

@dataclass(frozen=True)
class SurfaceRamp:
    peak: float
    ramp_length: float
    plateau_length: float
    center_x: float = 0.0

@dataclass(frozen=True)
class TractionNodes:
    """Direct nodal tractions, matching PREBEM TRACTION.NODES."""
    values: tuple[tuple[int, float, float], ...]


@dataclass(frozen=True)
class InternalPressure:
    pressure: float
    boundary: int = 1


@dataclass(frozen=True)
class ExplosionPoint:
    x: float
    y: float
    wave_velocity: float
    exponent: float
    frequency: float
    decay: float
    sv_amplitude: float = 0.0
    p_amplitude: float = 1.0
    formula: int = 1
    coefficient: float = 0.0
    secondary_exponent: float = 0.0


@dataclass(frozen=True)
class InputDeck:
    title: str
    domain: int
    analysis: int
    units: int | None
    sections: tuple[Section, ...]
    material: MaterialInput
    steps: int
    total_time: float
    load: PlaneWave | ExplosionPoint | SurfaceRamp | TractionNodes
    plots: tuple[int, ...]
    node_plots: tuple[int, ...]
    inside_pressure: float = 0.0
    free_field_stresses: int | None = None
    simulation: int | None = None
    element_nodes: tuple[int, ...] = ()
    internal_points: tuple[tuple[float, float], ...] = ()
    prescribed_conditions: tuple[tuple[int, int, int], ...] = ()
    gaussian_order: int = 6
    draw_zones: tuple[int, ...] = ()
    zones: tuple[ZoneInput, ...] = ()


def _numbers(line: str) -> list[float]:
    try:
        return [float(value.replace("D", "E").replace("d", "e")) for value in line.split()]
    except ValueError as error:
        raise ValueError(f"expected numerical data, got {line!r}") from error


def read_input_deck(path: str | Path) -> InputDeck:
    """Parse a PREBEM ``.DAT`` file used by the included reference cases."""

    raw_lines = Path(path).read_text(encoding="ascii").splitlines()
    lines = [line.strip() for line in raw_lines if line.strip() and not line.lstrip().startswith("*")]
    position = 0

    def take(expected: str | None = None) -> str:
        nonlocal position
        if position >= len(lines):
            raise ValueError("unexpected end of input deck")
        line = lines[position]
        position += 1
        if expected is not None and line.upper() != expected:
            raise ValueError(f"expected {expected!r}, got {line!r}")
        return line

    take("TITLE")
    title = take()
    take("TYPE.PROBLEM")
    type_labels = take().upper().split()
    type_values = _numbers(take())
    if len(type_labels) != len(type_values):
        raise ValueError("TYPE.PROBLEM labels and values have different lengths")
    type_data = dict(zip(type_labels, map(int, type_values), strict=True))

    take("SECTION.TYPE")
    line = take()
    if line.upper() == "BOUNDARIES":
        boundary_count = int(_numbers(take())[0])
        labels = take().upper().split()
    else:
        boundary_count = 1
        labels = line.upper().split()

    sections: list[Section] = []
    for _ in range(boundary_count):
        values = _numbers(take())
        if len(values) < len(labels):
            raise ValueError("SECTION.TYPE has fewer values than labels")
        fields = dict(zip(labels, values[: len(labels)], strict=True))
        geometry_index = labels.index("GEOMETRY")
        geometry = tuple(values[geometry_index:])
        shape = int(fields["SHAPE"])
        region = int(fields.get("REGION", 1))
        segments: tuple[CompoundSegment, ...] = ()
        if shape < 0:
            count = -shape
            # Current archived decks use an abbreviated form for a single
            # straight segment: -1 REGION ELEMENTS X0 Y0 X1 Y1.  Later
            # PREBEM sources use: -N REGION SEGMENT_SHAPE ELEMENTS ... .
            remainder = values[2:] if "REGION" in labels else values[1:]
            parsed: list[CompoundSegment] = []
            if count == 1 and (not remainder or int(remainder[0]) not in (1, 4)):
                parsed.append(CompoundSegment(4, int(fields["ELEMENTS"]), geometry))
            else:
                position_in_segment = 0
                first_remainder = list(remainder)
                for segment_index in range(count):
                    tokens = first_remainder if segment_index == 0 else _numbers(take())
                    if len(tokens) < 2:
                        raise ValueError("compound segment requires shape and element count")
                    segment_shape, segment_elements = int(tokens[0]), int(tokens[1])
                    required = 4 if segment_shape == 4 else 7 if segment_shape == 1 else 0
                    if required == 0 or len(tokens) < 2 + required:
                        raise ValueError("compound segment shape must be 4 (line) or 1 (ellipse)")
                    parsed.append(CompoundSegment(segment_shape, segment_elements, tuple(tokens[2:2 + required])))
                    position_in_segment += 1
            if any(segment.elements <= 0 for segment in parsed):
                raise ValueError("compound segment element counts must be positive")
            segments = tuple(parsed)
            section = Section(shape, region, sum(segment.elements for segment in segments), (), segments=segments)
        else:
            section = Section(
                shape=shape,
                region=region,
                x_elements=int(fields["ELEMENTS"]),
                geometry=geometry,
                y_elements=int(fields["YELEMENTS"]) if "YELEMENTS" in fields else None,
            )
        if section.shape == 0 or section.x_elements <= 0 or (not geometry and not section.segments):
            raise ValueError("invalid section shape, element count, or geometry")
        sections.append(section)

    take("MATERIAL.PROPERTIES")
    material_labels = take().upper().split()
    material_values = _numbers(take())
    material_fields = dict(zip(material_labels, material_values, strict=True))
    material = MaterialInput(
        young=material_fields["E"],
        poisson=material_fields["NU"],
        density=material_fields["RO"],
        viscosity=material_fields.get("VISC", 0.0),
        viscosity_model=int(material_fields.get("MODE", 1)),
    )

    take("STEP.NUMBER")
    steps = int(_numbers(take())[0])
    take("TIME.TOTAL")
    total_time = _numbers(take())[0]
    take("LOAD.STRESSES")
    load_kind = take().upper()
    if load_kind == "TRACTION.NODES":
        ncases = int(_numbers(take())[0])
        labels = take().upper().split()
        if not {"XTRACTION", "YTRACTION", "NODE"}.issubset(labels):
            raise ValueError("TRACTION.NODES requires XTRACTION YTRACTION NODE")
        values: list[tuple[int, float, float]] = []
        for _ in range(ncases):
            nums = _numbers(take())
            if len(nums) < 3:
                raise ValueError("TRACTION.NODES case requires x traction, y traction, and node list")
            xtraction, ytraction = nums[:2]
            for node in nums[2:]:
                inode = int(node)
                if inode == 0:
                    break
                values.append((inode, xtraction, ytraction))
        load = TractionNodes(tuple(values))
        # No additional load-value line is consumed for this load kind.
    elif load_kind == "SURFACE.RAMP":
        load_labels = take().upper().split()
        load_values = _numbers(take())
        if len(load_labels) != len(load_values):
            raise ValueError(f"{load_kind} labels and values have different lengths")
        fields = dict(zip(load_labels, load_values, strict=True))
        required = {"PEAK", "RAMP.LENGTH", "PLATEAU.LENGTH"}
        if not required.issubset(fields):
            raise ValueError(f"SURFACE.RAMP is missing {sorted(required - fields.keys())}")
        load = SurfaceRamp(fields["PEAK"], fields["RAMP.LENGTH"], fields["PLATEAU.LENGTH"], fields.get("CENTER.X", 0.0))
    elif load_kind in ("PLANE.WAVE", "HARMONIC.WAVE"):
        load_labels = take().upper().split()
        load_values = _numbers(take())
        if len(load_labels) != len(load_values):
            raise ValueError(f"{load_kind} labels and values have different lengths")
        fields = dict(zip(load_labels, load_values, strict=True))
        allowed = {"LONGITUDINAL.STRESS", "TRANSVERSE.SHEAR", "TRANSVERSE.STRESS",
                   "VERTICAL.STRESS", "ANGLE", "FREQUENCY", "XDISTANCE", "YDISTANCE"}
        if fields.keys() - allowed:
            raise ValueError(f"unknown plane-wave fields: {sorted(fields.keys() - allowed)}")
        if ("XDISTANCE" in fields) != ("YDISTANCE" in fields):
            raise ValueError("XDISTANCE and YDISTANCE must be provided together")
        load = PlaneWave(fields.get("LONGITUDINAL.STRESS", 0.0),
                         fields.get("TRANSVERSE.SHEAR", fields.get("TRANSVERSE.STRESS", 0.0)),
                         fields.get("ANGLE", 0.0), fields.get("VERTICAL.STRESS", 0.0),
                         fields.get("FREQUENCY", 0.0),
                         fields.get("XDISTANCE"), fields.get("YDISTANCE"))
        if load.frequency < 0:
            raise ValueError("wave frequency must be nonnegative")
    elif load_kind == "EXPLOSION.POINT":
        load_labels = take().upper().split()
        load_values = _numbers(take())
        if len(load_labels) != len(load_values):
            raise ValueError(f"{load_kind} labels and values have different lengths")
        fields = dict(zip(load_labels, load_values, strict=True))
        # Appendix F makes all explosion parameters except XE/YE/W optional;
        # formula-specific defaults are applied below.
        required = {"XE", "YE", "W"}
        if not required.issubset(fields):
            raise ValueError(f"EXPLOSION.POINT is missing {sorted(required - fields.keys())}")
        formula = int(fields.get("FORMULA", 1))
        if formula not in (1, 2, 3):
            raise ValueError("EXPLOSION.POINT FORMULA must be 1, 2, or 3")
        defaults = {1: (160.0, 2.75, 0.0), 2: (2.0, 1.57, 0.0), 3: (1.4, 0.70, 1.50)}
        default_b, default_n, default_n1 = defaults[formula]
        coefficient = fields.get("B", default_b)
        exponent = fields.get("N", default_n)
        secondary_exponent = fields.get("N1", default_n1)
        # PREMOD2 converts the empirical coefficient and charge to the SI
        # quantities consumed by ZEXPLOSQ/TEXPLOSQ.  The supplied decks omit
        # UNITS, for which the original program uses SI input (kgf and metres).
        units = type_data.get("UNITS", 1)
        charge = fields["W"]
        x, y = fields["XE"], fields["YE"]
        if units == 1:
            charge *= 9.8066
        elif units == 2:
            charge *= 4.4482
            x *= 0.3048; y *= 0.3048
        if formula == 1:
            coefficient *= 0.3048 ** (exponent + 1.0) / 4.4482 ** (exponent / 3.0)
        elif formula == 2:
            coefficient *= 0.3048 ** (exponent + 1.0) / 4.4482 ** (exponent / 2.0)
        else:
            coefficient *= 9.8066 * exponent
        load = ExplosionPoint(
            x, y, charge, exponent, fields.get("F", 1.0), fields.get("DECAY", 1.0),
            fields.get("SV", 0.0), fields.get("P", 1.0), formula,
            coefficient, secondary_exponent,
        )
    else:
        raise ValueError(f"unsupported load module {load_kind!r}")

    plots: tuple[int, ...] = ()
    node_plots: tuple[int, ...] = ()
    free_field_stresses: int | None = None
    simulation: int | None = None
    zones: tuple[ZoneInput, ...] = ()
    gaussian_order = 6
    draw_zones: tuple[int, ...] = ()
    element_nodes: tuple[int, ...] = ()
    internal_points: tuple[tuple[float, float], ...] = ()
    prescribed_conditions: tuple[tuple[int, int, int], ...] = ()
    inside_pressure = 0.0
    while True:
        module = take().upper()
        if module == "END":
            break
        if module == "INSIDE.PRESSURE":
            count = int(_numbers(take())[0])
            if count != 1:
                raise ValueError("INSIDE.PRESSURE currently supports one pressure case")
            if take().upper() != "BOUNDARY":
                raise ValueError("INSIDE.PRESSURE requires BOUNDARY")
            boundary = int(_numbers(take())[0])
            if take().upper() != "INSIDE.PRESSURE":
                raise ValueError("INSIDE.PRESSURE requires a pressure value")
            inside_pressure = float(_numbers(take())[0])
            continue
        if module == "INTERNAL.POINTS":
            count = int(_numbers(take())[0])
            labels = take().upper().split()
            if count < 0 or count > 10 or set(labels) != {"XI", "YI"}:
                raise ValueError("INTERNAL.POINTS requires up to 10 XI/YI coordinate rows")
            points = []
            for _ in range(count):
                values = _numbers(take())
                if len(values) != 2:
                    raise ValueError("each INTERNAL.POINTS row requires two coordinates")
                fields = dict(zip(labels, values, strict=True))
                scale = 0.3048 if type_data.get("UNITS") == 2 else 1.0
                points.append((fields["XI"] * scale, fields["YI"] * scale))
            internal_points = tuple(points)
            continue
        if module == "PRESCRIBED.CONDITIONS":
            labels = take().upper().split()
            if not labels or labels[0] != "TOTAL.NODES":
                raise ValueError("PRESCRIBED.CONDITIONS must start with TOTAL.NODES")
            count = int(_numbers(take())[0])
            if count < 0:
                raise ValueError("TOTAL.NODES must be non-negative")
            entry_labels = labels[1:] if len(labels) > 1 else take().upper().split()
            required = {"POSITION.NODE", "AXIS.NUMBER", "CONDITION.CODE"}
            if set(entry_labels) != required:
                raise ValueError("PRESCRIBED.CONDITIONS requires POSITION.NODE AXIS.NUMBER CONDITION.CODE")
            conditions = []
            for _ in range(count):
                values = _numbers(take())
                if len(values) != 3:
                    raise ValueError("each prescribed-condition row requires three integer values")
                fields = dict(zip(entry_labels, map(int, values), strict=True))
                node, axis, condition = fields["POSITION.NODE"], fields["AXIS.NUMBER"], fields["CONDITION.CODE"]
                if node <= 0 or axis not in (1, 2) or condition not in (0, 1):
                    raise ValueError("invalid prescribed boundary condition")
                conditions.append((node, axis, condition))
            prescribed_conditions = tuple(conditions)
            continue
        if module == "ZONES":
            count = int(_numbers(take())[0])
            labels = take().upper().split()
            required = {"NUMBER", "E", "NU", "RO", "VISC", "MODE"}
            boundary_labels = [label for label in labels if label.startswith("BOU") or label == "BOUNDARY"]
            if not required.issubset(labels) or not boundary_labels:
                raise ValueError("ZONES requires zone number, boundary labels, E, NU, RO, VISC, and MODE")
            parsed_zones = []
            for _ in range(count):
                values = _numbers(take())
                if len(values) != len(labels):
                    raise ValueError("each ZONES row requires eight values")
                fields = dict(zip(labels, values, strict=True))
                raw_boundaries = tuple(int(fields[label]) for label in boundary_labels)
                boundaries = tuple(value for value in raw_boundaries if value != 0)
                if not boundaries:
                    raise ValueError("each ZONES row must identify at least one boundary")
                parsed_zones.append(ZoneInput(int(fields["NUMBER"]), raw_boundaries[0], raw_boundaries[1] if len(raw_boundaries) > 1 else 0,
                                              MaterialInput(fields["E"], fields["NU"], fields["RO"], fields["VISC"], int(fields["MODE"])), boundaries))
            zones = tuple(parsed_zones)
            continue
        if module == "GAUSSIAN.QUADRATURE":
            gaussian_order = int(_numbers(take())[0])
            if gaussian_order < 1:
                raise ValueError("GAUSSIAN.QUADRATURE order must be positive")
            continue
        if module == "DRAW.ZONES":
            draw_zones = tuple(map(int, _numbers(take())))
            continue
        values = _numbers(take())
        if module == "PLOTS":
            plots = tuple(map(int, values))
        elif module == "NODE.PLOTS":
            node_plots = tuple(map(int, values))
        elif module == "FREE.FIELD.STRESSES":
            free_field_stresses = int(values[0])
        elif module == "SIMULATION":
            simulation = int(values[0])
        elif module == "ELEMENT.NODES":
            element_nodes = tuple(map(int, values))
        else:
            raise ValueError(f"unsupported trailing module {module!r}")
    if position != len(lines):
        raise ValueError(f"unexpected data after END: {lines[position]!r}")

    return InputDeck(
        title=title,
        domain=type_data["DOMAIN"],
        analysis=type_data["ANALYSIS"],
        units=type_data.get("UNITS"),
        sections=tuple(sections),
        material=material,
        steps=steps,
        total_time=total_time,
        load=load,
        plots=plots,
        node_plots=node_plots,
        inside_pressure=inside_pressure,
        free_field_stresses=free_field_stresses,
        simulation=simulation,
        element_nodes=element_nodes,
        internal_points=internal_points,
        prescribed_conditions=prescribed_conditions,
        zones=zones,
        gaussian_order=gaussian_order,
        draw_zones=draw_zones,
    )
