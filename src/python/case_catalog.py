"""Discovery and presentation metadata for the historical example decks."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from .input_deck import ExplosionPoint, InputDeck, PlaneWave, SurfaceRamp, TractionNodes, read_input_deck

@dataclass(frozen=True)
class ExampleCase:
    path: Path
    deck: InputDeck
    family: str
    label: str
    thesis_section: str
    subcategory: str


_THESIS_EXAMPLES = {
    "THESIS_RECT_EXPLOSION": ("2.3.2", "2.3.2.2  Rectangular Cavity Loaded by an Explosion P Wave - Figures 2.40-2.43"),
    "INFPL1": ("2.3.1", "2.3.1.1  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave"),
    "INFPL391": ("2.3.1", "2.3.1.2a  Square Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave — Equal Elements"),
    "INFPL392": ("2.3.1", "2.3.1.2b  Square Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave — Unequal Elements"),
    "IFPLREE1": ("2.3.1", "2.3.1.3  Rectangular Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave"),
    "INFPLSV1": ("2.3.1", "2.3.1.4  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane SV Wave"),
    "INFPLPS1": ("2.3.1", "2.3.1.5  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane P+SV Wave"),
    "INFEXP1": ("2.3.2", "2.3.2.1a  Circular Cross Section Cavity in a Full Space Loaded by an Explosion Wave - Cubic Root Formula Without Decay in Time"),
    "TINFPL1": ("3.3.1", "3.3.1.1  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave"),
    "TIFPLSV1": ("3.3.1", "3.3.1.2  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane SV Wave"),
    "TIFPLPS1": ("3.3.1", "3.3.1.3  Circular Cross Section Cavity in an Infinite Medium Swept by a Plane P+SV Wave"),
}
# The published unequal-square graphs use duplicate corner nodes and refer to
# nodes through 84.  INFPL394 is the corresponding discontinuous-corner deck.
_THESIS_EXAMPLES.pop("INFPL392")
_THESIS_EXAMPLES["INFPL394"] = (
    "2.3.1", "2.3.1.2b  Square Cross Section Cavity in an Infinite Medium Swept by a Plane P Wave - Unequal Elements",
)

_THESIS_CATEGORIES = {
    "2.3.1": "2. LAPLACE DOMAIN",
    "2.3.2": "2. LAPLACE DOMAIN",
    "3.3.1": "3. TIME DOMAIN",
}

def project_root() -> Path:
    return Path(__file__).resolve().parents[2]

def example_directory(root: Path | None = None) -> Path:
    base = root or project_root()
    return base / "examples" / "in"

def thesis_figure_directory(root: Path | None = None) -> Path:
    """Return the directory containing the complete thesis-figure DAT aliases."""
    base = root or project_root()
    return example_directory(base)

def thesis_figure_dat_files(root: Path | None = None) -> tuple[Path, ...]:
    """List all thesis-figure DAT aliases, including corrected ``b`` files."""
    return tuple(sorted(thesis_figure_directory(root).glob("*.DAT")))

def _thesis_section(deck: InputDeck, stem: str) -> str:
    """Return the thesis application number for recognized examples."""
    entry = _THESIS_EXAMPLES.get(stem.upper())
    return entry[1].split()[0].rstrip("ab") if entry else ""

def _family(deck: InputDeck, section: str = "") -> str:
    entry = _THESIS_EXAMPLES.get(section.upper().replace(".", ""))
    medium = "Semi-infinite" if any(s.region != 1 for s in deck.sections) or "SEMI" in deck.title.upper() else "Infinite"
    load = "Explosion" if isinstance(deck.load, ExplosionPoint) else "Plane wave"
    domain = "Laplace" if deck.domain == 1 else "CQ time" if deck.domain == 5 else "Direct time"
    if section:
        paragraph = ".".join(section.split(".")[:3])
        return _THESIS_CATEGORIES.get(paragraph, f"{paragraph} {domain} / {load} / {medium}")
    return f"Reference case  {domain} / {load} / {medium}"

def discover_examples(root: Path | None = None) -> tuple[ExampleCase, ...]:
    cases = []
    base = root or project_root()
    paths = {path.stem.upper(): path for path in sorted(example_directory(base).glob("*.DAT"))}
    rectangle = example_directory(base) / "THESIS_RECT_EXPLOSION.DAT"
    if rectangle.is_file():
        paths[rectangle.stem.upper()] = rectangle
    # TINFPL1 is the exact Chapter 3.3.1.1 thesis deck but is stored with the
    # retained executable examples, not in the adopted 1992 corpus.
    canonical = example_directory(base) / "TINFPL1.DAT"
    if canonical.is_file():
        paths[canonical.stem.upper()] = canonical
    # Present examples in their thesis order, rather than file-name order.
    for path in sorted(paths.values(), key=lambda item: _THESIS_EXAMPLES.get(item.stem.upper(), ("9", item.stem.upper()))[1]):
        try:
            deck = read_input_deck(path)
        except (OSError, ValueError):
            continue
        mapping = _THESIS_EXAMPLES.get(path.stem.upper())
        if not mapping:
            continue
        paragraph, label = mapping
        section = label.split()[0].rstrip("ab")
        subcategory = f"{paragraph} {'Plane Wave' if isinstance(deck.load, PlaneWave) else 'Explosion at a Point'}"
        cases.append(ExampleCase(path, deck, _THESIS_CATEGORIES[paragraph], label, section, subcategory))
    return tuple(cases)

def identify_input_case(path: Path, deck: InputDeck, examples: tuple[ExampleCase, ...]) -> ExampleCase:
    """Identify browsed decks, including renamed copies of catalog examples."""
    known = next((case for case in examples if case.path.resolve() == path.resolve()), None)
    if known is None or known.deck != deck:
        known = next((case for case in examples if case.deck == deck), None)
    if known is not None:
        label = known.label if known.path.resolve() == path.resolve() else f"{path.name} — {known.label}"
        return ExampleCase(path, deck, known.family, label, known.thesis_section, known.subcategory)

    load = "Plane Wave" if isinstance(deck.load, PlaneWave) else "Explosion at a Point"
    family = {0: "1. ELASTOSTATIC DOMAIN", 1: "2. LAPLACE DOMAIN", 2: "3. TIME DOMAIN",
              5: "5. MULTI-ZONAL BEM"}.get(
        deck.domain, "CUSTOM INPUT FILE")
    paragraph = {1: "2.3", 2: "3.3"}.get(deck.domain)
    subcategory = f"{paragraph}.{1 if isinstance(deck.load, PlaneWave) else 2} {load}" if paragraph else load
    if deck.zones:
        family, subcategory = "5. MULTI-ZONAL BEM", f"Multi-zone / {load}"
    title = " ".join(deck.title.split())
    label = f"{path.name} — {title}" if title else path.name
    return ExampleCase(path, deck, family, label, "", subcategory)


def case_details(case: ExampleCase) -> str:
    deck, load = case.deck, case.deck.load
    description = ""
    if case.thesis_section and isinstance(load, PlaneWave):
        if len(deck.sections) == 1 and deck.sections[0].shape == 1:
            description += f"The radius of the cavity is {deck.sections[0].geometry[0]:.1f} m.\n\n"
        if load.longitudinal_stress and load.transverse_shear:
            description += (
                "The incident waves are shock compressional and shear (SV) plane waves "
                f"of {load.longitudinal_stress:.1f} MPa and {load.transverse_shear:.1f} MPa "
                "magnitude, respectively, that are suddenly applied and maintained constant in time. "
            )
        elif load.longitudinal_stress or load.transverse_shear:
            wave = "compressional" if load.longitudinal_stress else "shear (SV)"
            magnitude = load.longitudinal_stress or load.transverse_shear
            description += (
                f"The incident wave is a shock {wave} plane wave of {magnitude:.1f} MPa "
                "magnitude that is suddenly applied and maintained constant in time. "
            )
        description += (
            f"The angle of incidence of the incoming wave is {load.angle_degrees:g}°, "
            "measured from the global horizontal x axis."
        )
    if isinstance(load, PlaneWave):
        loading = f"Plane wave: P={load.longitudinal_stress:g}, SV={load.transverse_shear:g}, angle={load.angle_degrees:g} deg"
    elif isinstance(load, SurfaceRamp):
        loading = (f"Surface ramp: peak={load.peak:g} MPa, ramp length={load.ramp_length:g} m, "
                   f"plateau length={load.plateau_length:g} m, centre x={load.center_x:g} m")
    elif isinstance(load, TractionNodes):
        loading = f"Direct nodal tractions: {len(load.values)} node assignments"
    else:
        loading = f"Explosion at ({load.x:g}, {load.y:g}); P={load.p_amplitude:g}, SV={load.sv_amplitude:g}, decay={load.decay:g}"
    section_line = f"Thesis test: {case.thesis_section}\n" if case.thesis_section else ""
    return (f"{case.path.name}\n{deck.title.strip()}\n\n{section_line}Category: {case.family}\n"
            f"Boundaries: {len(deck.sections)}   Steps: {deck.steps}   Total time: {deck.total_time:g}\n"
            f"Material: E={deck.material.young:g}, nu={deck.material.poisson:g}, density={deck.material.density:g}\n{loading}"
            + (f"\n\n{description}" if description else ""))
