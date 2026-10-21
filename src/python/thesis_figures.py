"""Figure numbers and captions for circular plane-wave thesis cases."""
from __future__ import annotations
from dataclasses import dataclass
from .input_deck import InputDeck, PlaneWave

@dataclass(frozen=True)
class ThesisFigures:
    stress_history: tuple[str, str]
    boundary_profile: tuple[str, str]
    deformed_shape: tuple[str, str]
    x_displacement: tuple[str, str]
    y_displacement: tuple[str, str] | None
    featured_step: int


@dataclass(frozen=True)
class ThesisGraph:
    """One graph prescribed by a thesis application example."""
    figure: str
    caption: str
    view: str
    quantity: str
    component: str = "X"
    step: int | None = None
    nodes: tuple[int, ...] | None = None


def thesis_graphs_for_case(stem: str, deck: InputDeck | None) -> tuple[ThesisGraph, ...]:
    """Exact graph selections documented for the translated plane-wave examples."""
    stem = stem.upper()
    graphs: dict[str, tuple[ThesisGraph, ...]] = {
        "IN-5.10A-LD": (
            ThesisGraph("5.10", "Tangential stresses, stiff medium-stiffer liner (LD)",
                        "Time history", "hoop stress", nodes=(1, 9, 33, 41)),
        ),
        "IN-5.10B-TD": (
            ThesisGraph("5.10", "Tangential stresses, stiff medium-stiffer liner (TD comparison)",
                        "Time history", "hoop stress", nodes=(1, 9, 33, 41)),
        ),
        "IN-5.10B-CQ": (
            ThesisGraph("5.10", "Tangential stresses, stiff medium-stiffer liner (CQ comparison)",
                        "Time history", "hoop stress", nodes=(1, 9, 33, 41)),
        ),
        "IN-5.11-LD": (
            ThesisGraph("5.11", "Tangential stresses, soft medium-stiff liner",
                        "Time history", "hoop stress", nodes=(1, 9, 33, 41)),
        ),
        "IN-5.13-TD": (
            ThesisGraph("5.13", "Tangential stresses, stiff liner under point explosion",
                        "Time history", "hoop stress", nodes=(5, 9, 37, 41)),
        ),
        "IN-5.13-LD": (
            ThesisGraph("5.13", "Tangential stresses, stiff liner under point explosion (LD)",
                        "Time history", "hoop stress", nodes=(5, 9, 37, 41)),
        ),
        "IN-5.13-CQ": (
            ThesisGraph("5.13", "Tangential stresses, stiff liner under point explosion (BDF2 CQ)",
                        "Time history", "hoop stress", nodes=(5, 9, 37, 41)),
        ),
        "THESIS_RECT_EXPLOSION": (
            ThesisGraph("2.40", "Tangential stresses at the 6th time step", "Boundary profile", "hoop stress", step=6),
            ThesisGraph("2.41", "Deformed shape at the 6th time step", "Deformed geometry", "displacement", step=6),
            ThesisGraph("2.42", "Tangential stresses versus time steps", "Time history", "hoop stress", nodes=(68, 70, 71, 72, 84, 73, 74)),
            ThesisGraph("2.43", "X displacements versus time steps", "Time history", "displacement", nodes=(68, 70, 71, 72, 84, 73, 74)),
        ),
        "INFPL1": (
            ThesisGraph("2.1", "Variation of tangential stresses in time", "Time history", "hoop stress"),
            ThesisGraph("2.2", "Tangential stresses at time step 14, for all nodes", "Boundary profile", "hoop stress", step=14),
            ThesisGraph("2.3", "Deformed shape of the mesh at time step 14", "Deformed geometry", "displacement", step=14),
            ThesisGraph("2.4", "X displacements versus time steps", "Time history", "displacement"),
        ),
        "INFPL391": (
            ThesisGraph("2.5", "Variation of tangential stresses in time", "Time history", "hoop stress"),
            ThesisGraph("2.6", "Hoop stresses at time step 7, for all nodes", "Boundary profile", "hoop stress", step=7),
            ThesisGraph("2.7", "Deformed shape of the mesh at time step 7", "Deformed geometry", "displacement", step=7),
            ThesisGraph("2.8", "X displacements in time", "Time history", "displacement"),
        ),
        "INFPL394": (
            ThesisGraph("2.10", "Tangential stresses at time step 7, for all nodes", "Boundary profile", "hoop stress", step=7),
            ThesisGraph("2.11", "Deformed shape of the mesh at time step 7", "Deformed geometry", "displacement", step=7),
            ThesisGraph("2.12", "Variation of hoop stress in time", "Time history", "hoop stress", nodes=(1, 41, 51, 61, 71, 81, 84)),
            ThesisGraph("2.13", "Variation of hoop stress in time", "Time history", "hoop stress", nodes=(49, 50, 51, 52, 53, 81)),
            ThesisGraph("2.14", "Variation of hoop stress in time", "Time history", "hoop stress", nodes=(69, 70, 71, 72, 73, 84)),
            ThesisGraph("2.15", "Variation of X displacement in time", "Time history", "displacement", nodes=(51, 81, 84, 71, 1, 41, 61)),
        ),
        "IFPLREE1": (
            ThesisGraph("2.17", "Distribution of tangential stresses at time step 7", "Boundary profile", "hoop stress", step=7),
            ThesisGraph("2.18", "Variation of hoop stress in time", "Time history", "hoop stress", nodes=(1, 41, 50, 72, 81, 84)),
            ThesisGraph("2.19", "Variation of X displacements in time", "Time history", "displacement", nodes=(1, 41, 50, 72, 81, 84)),
        ),
        "INFPLSV1": (
            ThesisGraph("2.20", "Variation of tangential stresses in time", "Time history", "hoop stress"),
            ThesisGraph("2.21", "Tangential stresses on the cavity boundary at time step 14", "Boundary profile", "hoop stress", step=14),
            ThesisGraph("2.22", "X displacements in time", "Time history", "displacement", "X"),
            ThesisGraph("2.23", "Y displacements in time", "Time history", "displacement", "Y"),
        ),
        "INFPLPS1": (
            ThesisGraph("2.24", "Variation of hoop stresses in time", "Time history", "hoop stress"),
            ThesisGraph("2.25", "Tangential stresses at time step 14", "Boundary profile", "hoop stress", step=14),
            ThesisGraph("2.26", "X displacements in time", "Time history", "displacement", "X"),
            ThesisGraph("2.27", "Y displacements in time", "Time history", "displacement", "Y"),
        ),
        "INFEXP1": (
            ThesisGraph("2.28", "Tangential stresses at the 6th time step", "Boundary profile", "hoop stress", step=6),
            ThesisGraph("2.29", "Deformed shape of the mesh at the 6th time step", "Deformed geometry", "displacement", step=6),
            ThesisGraph("2.30", "Variation of hoop stresses and displacements in time at nodes 1, 9, 11, 13, 15", "Time history", "hoop stress", nodes=(1, 9, 11, 13, 15)),
            ThesisGraph("2.31", "Peak values of displacements at the 9th time step at node 1", "Time history", "displacement", nodes=(1,)),
        ),
    }
    if stem in graphs:
        return graphs[stem]
    figures = chapter3_figures(deck)
    if figures is None:
        return ()
    result = [
        ThesisGraph(figures.stress_history[0], figures.stress_history[1], "Time history", "hoop stress"),
        ThesisGraph(figures.boundary_profile[0], figures.boundary_profile[1], "Boundary profile", "hoop stress", step=figures.featured_step),
        ThesisGraph(figures.deformed_shape[0], figures.deformed_shape[1], "Deformed geometry", "displacement", step=figures.featured_step),
        ThesisGraph(figures.x_displacement[0], figures.x_displacement[1], "Time history", "displacement", "X"),
    ]
    if figures.y_displacement:
        result.append(ThesisGraph(figures.y_displacement[0], figures.y_displacement[1], "Time history", "displacement", "Y"))
    return tuple(result)

def chapter3_figures(deck: InputDeck | None) -> ThesisFigures | None:
    """Return the matching Chapter 2 or Chapter 3 figure set."""
    # The published figure sets below are only for the circular-cavity cases.
    if deck is None or not isinstance(deck.load, PlaneWave) or len(deck.sections) != 1 or deck.sections[0].shape != 1:
        return None
    p = deck.load.longitudinal_stress != 0.0
    sv = deck.load.transverse_shear != 0.0
    if deck.domain == 1:
        if p and sv:
            return ThesisFigures(("2.24", "Tangential stresses versus time steps"),
                ("2.25", "Tangential stresses at 14th time step"),
                ("2.25", "Deformed shape at 14th time step"),
                ("2.26", "x displacements versus time steps"),
                ("2.27", "y displacements versus time steps"), 14)
        if sv:
            return ThesisFigures(("2.20", "Tangential stresses versus time steps"),
                ("2.21", "Tangential stresses at 14th time step"),
                ("2.21", "Deformed shape at 14th time step"),
                ("2.22", "x displacements versus time steps"),
                ("2.23", "y displacements versus time steps"), 14)
        return ThesisFigures(("2.1", "Tangential stresses versus time steps"),
            ("2.2", "Tangential stresses at 14th time step"),
            ("2.3", "Deformed shape at 14th time step"),
            ("2.4", "x displacements versus time steps"), None, 14)
    if p and sv:
        return ThesisFigures(("3.10", "Tangential stresses versus time steps"),
            ("3.11", "Tangential stresses at 11th time step"),
            ("3.11", "Deformed shape at 11th time step"),
            ("3.12", "x displacements versus time steps"),
            ("3.13", "y displacements versus time steps"), 11)
    if sv:
        return ThesisFigures(("3.6", "Tangential stresses versus time steps"),
            ("3.7", "Tangential stresses at 12th time step"),
            ("3.7", "Deformed shape at 12th time step"),
            ("3.8", "x displacements versus time steps"),
            ("3.9", "y displacements versus time steps"), 12)
    return ThesisFigures(("3.2", "Tangential stresses versus time steps"),
        ("3.3", "Tangential stresses at 12th time step"),
        ("3.4", "Deformed shape at 12th time step"),
        ("3.5", "x displacements versus time steps"), None, 12)
