"""Transactional editor for the simple boundary shapes supported by PREBEM."""

from dataclasses import replace
from tkinter import Canvas, StringVar, Toplevel, messagebox, ttk

import numpy as np

from .input_deck import InputDeck, Section
from .mesh import mesh_section


SHAPES = {1: "Ellipse / circle", 2: "Square", 3: "Rectangle", 4: "Straight line", -1: "Combined shape (-n segments)"}
# Labels and defaults follow the geometry tuple consumed by mesh_section.
GEOMETRY_FIELDS = {
    1: (("Radius / semiaxis", 5), ("Axis ratio b/a", 1), ("Centre X", 0),
        ("Centre Y", 0), ("Rotation (degrees)", 0)),
    2: (("Side length", 5), ("Centre X", 0), ("Centre Y", 0),
        ("Rotation (degrees)", 0), ("X grading", 1), ("Y grading", 1),
        ("Corner nodes (DIS: 1 shared, 2 separate)", 2)),
    3: (("X-side length", 5), ("Y-side length", 4), ("Centre X", 0),
        ("Centre Y", 0), ("Rotation (degrees)", 0), ("X grading", 1),
        ("Y grading", 1), ("Corner nodes (DIS: 1 shared, 2 separate)", 2)),
    4: (("Start X", 0), ("Start Y", 0), ("End X", 10), ("End Y", 0),
        ("Grading", 1)),
}
REGIONS = ("1: Cavity", "2: Inclusion", "3: Open boundary", "4: Free surface",
           "5: Enclosing open boundary", "6: Enclosing closed boundary")


def geometry_deck(deck: InputDeck, sections: tuple[Section, ...]) -> InputDeck:
    """Validate edits without silently invalidating indexed DAT assignments."""
    if not sections:
        raise ValueError("Keep at least one boundary.")
    if sections == deck.sections:
        return deck
    meshes = []
    for number, section in enumerate(sections, 1):
        if section.x_elements < 1 or (section.y_elements is not None and section.y_elements < 1):
            raise ValueError(f"Boundary {number}: element counts must be positive integers.")
        if section.shape == 1 and section.x_elements < 3:
            raise ValueError(f"Boundary {number}: use at least 3 elements around a circle or ellipse.")
        if section.region not in range(1, 7) or not np.all(np.isfinite(section.geometry)):
            raise ValueError(f"Boundary {number}: invalid region or non-finite geometry.")
        if section.shape == 4 and section.region not in (3, 4):
            raise ValueError(f"Boundary {number}: straight lines require REGION=3 or REGION=4.")
        if section.shape in (1, 2, 3) and section.geometry[0] <= 0:
            raise ValueError(f"Boundary {number}: dimensions must be positive.")
        if section.shape in (1, 3) and len(section.geometry) > 1 and section.geometry[1] <= 0:
            raise ValueError(f"Boundary {number}: axis ratio or height must be positive.")
        corner_index = {2: 6, 3: 7}.get(section.shape)
        if corner_index is not None and len(section.geometry) > corner_index and section.geometry[corner_index] not in (1, 2):
            raise ValueError(f"Boundary {number}: corner nodes must be 1 or 2.")
        mesh = mesh_section(section)
        if not np.all(np.isfinite(mesh.coordinates)):
            raise ValueError(f"Boundary {number}: the generated mesh is not finite.")
        # Separate boundaries must not share nodes (separate corner nodes
        # within an individual boundary are intentional in some DAT files).
        nodes = {tuple(point) for point in np.round(mesh.coordinates, 10)}
        if any(nodes.intersection(previous) for previous in meshes):
            raise ValueError(f"Boundary {number} touches another boundary. Adjust its position or dimensions.")
        meshes.append(nodes)
    if deck.zones or deck.prescribed_conditions or deck.simulation is not None or deck.element_nodes:
        raise ValueError("This DAT has zone, simulation, or node assignments. Edit its geometry and assignments together in the DAT file.")
    # Node numbers and thesis presets refer to the old mesh. Use generic
    # output selection after geometry changes, without altering the DAT file.
    return replace(deck, sections=sections, node_plots=())


class GeometryEditor:
    def __init__(self, parent, deck: InputDeck, on_apply):
        self.deck = deck
        self.sections = list(deck.sections)
        self.on_apply = on_apply
        self.index = None
        self.window = Toplevel(parent)
        self.window.title("Edit boundaries")
        self.window.transient(parent)
        self.window.geometry("760x700")
        self.window.minsize(650, 600)
        footer = ttk.Frame(self.window, padding=8)
        footer.pack(side="bottom", fill="x")
        ttk.Button(footer, text="Apply geometry", command=self.apply).pack(side="right")
        ttk.Button(footer, text="Cancel", command=self.window.destroy).pack(side="right", padx=6)
        ttk.Label(self.window, text="Add boundaries, then edit each shape and position. Lengths use the DAT length unit.",
                  wraplength=700, padding=8).pack(fill="x")
        listing = ttk.Frame(self.window, padding=(8, 0))
        listing.pack(fill="x")
        self.tree = ttk.Treeview(listing, columns=("shape", "region", "elements"), show="tree headings", height=4)
        self.tree.heading("#0", text="Boundary")
        self.tree.column("#0", width=80, stretch=False)
        for key, title in (("shape", "Shape"), ("region", "Region"), ("elements", "X / Y elements")):
            self.tree.heading(key, text=title)
            self.tree.column(key, width=170)
        self.tree.pack(side="left", fill="x", expand=True)
        scroll = ttk.Scrollbar(listing, command=self.tree.yview)
        scroll.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.bind("<<TreeviewSelect>>", self.select)
        actions = ttk.Frame(self.window, padding=8)
        actions.pack(fill="x")
        ttk.Button(actions, text="Add boundary", command=self.add).pack(side="left")
        ttk.Button(actions, text="Remove selected", command=self.remove).pack(side="left", padx=6)
        form_area = ttk.Frame(self.window)
        form_area.pack(fill="both", expand=True)
        canvas = Canvas(form_area, highlightthickness=0)
        form_scroll = ttk.Scrollbar(form_area, command=canvas.yview)
        form_scroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        canvas.configure(yscrollcommand=form_scroll.set)
        self.form = ttk.Frame(canvas, padding=8)
        form_window = canvas.create_window(0, 0, window=self.form, anchor="nw")
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(form_window, width=event.width))
        self.form.bind("<Configure>", lambda event: canvas.configure(scrollregion=canvas.bbox("all")))
        self.refresh(0)
        self.window.grab_set()

    def refresh(self, index):
        self.tree.delete(*self.tree.get_children())
        for i, section in enumerate(self.sections):
            self.tree.insert("", "end", iid=str(i), text=str(i + 1),
                             values=(SHAPES.get(section.shape, "Compound (DAT)"), section.region,
                                     f"{section.x_elements} / {section.y_elements or section.x_elements}"))
        self.index = index
        self.tree.selection_set(str(index))
        self.tree.see(str(index))
        self.show(self.sections[index])

    def show(self, section):
        for widget in self.form.winfo_children():
            widget.destroy()
        self.fields = {}
        self.shape = StringVar(self.window, SHAPES.get(section.shape, "Compound (DAT)"))
        if section.shape < 0:
            ttk.Label(self.form, text=f"Combined shape: {-section.shape} segments",
                      font=("TkDefaultFont", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
            ttk.Label(self.form, text="Combined shapes consist of elliptical/circular and straight-line segments.\n"
                                      "Their segment definitions are preserved when geometry is applied; edit the DAT file\n"
                                      "to change segment coordinates, angles, radii, or segment element counts.",
                      justify="left").grid(row=1, column=0, columnspan=2, sticky="w", pady=(8, 8))
            for number, segment in enumerate(section.segments, 1):
                kind = "Straight line" if segment.shape == 4 else "Ellipse / circle"
                ttk.Label(self.form, text=f"Segment {number}: {kind}, {segment.elements} elements\n"
                                         f"Geometry: {', '.join(f'{value:g}' for value in segment.geometry)}",
                          justify="left").grid(row=number + 1, column=0, columnspan=2, sticky="w", pady=2)
            return
        ttk.Label(self.form, text="Shape").grid(row=0, column=0, sticky="w")
        shape = ttk.Combobox(self.form, textvariable=self.shape,
                             values=tuple(name for code, name in SHAPES.items() if code > 0), state="readonly")
        shape.grid(row=0, column=1, sticky="ew")
        shape.bind("<<ComboboxSelected>>", self.change_shape)
        self.region = StringVar(self.window, REGIONS[section.region - 1])
        ttk.Label(self.form, text="REGION").grid(row=1, column=0, sticky="w")
        ttk.Combobox(self.form, textvariable=self.region, values=REGIONS, state="readonly").grid(row=1, column=1, sticky="ew")
        values = [("ELEMENTS (X-side for rectangles)", section.x_elements)]
        if section.shape == 3:
            values.append(("YELEMENTS", section.y_elements or section.x_elements))
        geometry_values = list(section.geometry)
        values.extend((label, geometry_values[i] if i < len(geometry_values) else default)
                      for i, (label, default) in enumerate(GEOMETRY_FIELDS[section.shape]))
        for row, (label, value) in enumerate(values, 2):
            ttk.Label(self.form, text=label).grid(row=row, column=0, sticky="w", pady=2)
            variable = StringVar(self.window, f"{value:g}")
            ttk.Entry(self.form, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=2)
            self.fields[label] = variable
        self.form.columnconfigure(1, weight=1)

    def change_shape(self, _event=None):
        code = next(code for code, name in SHAPES.items() if name == self.shape.get())
        self.show(Section(code, 3 if code == 4 else 1, 8,
                          tuple(default for _, default in GEOMETRY_FIELDS[code])))

    def save_current(self):
        if self.index is None or self.shape.get() not in SHAPES.values():
            return True
        try:
            code = next(code for code, name in SHAPES.items() if name == self.shape.get())
            values = {key: float(value.get()) for key, value in self.fields.items()}
            if not all(np.isfinite(value) for value in values.values()):
                raise ValueError("Enter finite numbers in every field.")
            counts = [values["ELEMENTS (X-side for rectangles)"]]
            if code == 3:
                counts.append(values["YELEMENTS"])
            if any(count < 1 or not count.is_integer() for count in counts):
                raise ValueError("Element counts must be positive whole numbers.")
            geometry = tuple(values[label] for label, _ in GEOMETRY_FIELDS[code])
            section = Section(code, int(self.region.get().split(":")[0]), int(counts[0]), geometry,
                              int(counts[1]) if code == 3 else None)
            original = self.sections[self.index]
            if original.shape == code:
                expanded = tuple(original.geometry[i] if i < len(original.geometry) else default
                                 for i, (_, default) in enumerate(GEOMETRY_FIELDS[code]))
                if (expanded == geometry and original.region == section.region
                        and original.x_elements == section.x_elements
                        and (code != 3 or (original.y_elements or original.x_elements) == section.y_elements)):
                    return True
            # Validate the edited boundary before navigating to another row.
            geometry_deck(replace(self.deck, sections=(), zones=(), prescribed_conditions=(),
                                  simulation=None, element_nodes=()), (section,))
            self.sections[self.index] = section
            return True
        except (ValueError, OverflowError) as error:
            messagebox.showerror("Geometry error", str(error), parent=self.window)
            return False

    def select(self, _event=None):
        selected = self.tree.selection()
        if not selected or int(selected[0]) == self.index:
            return
        index = int(selected[0])
        if not self.save_current():
            self.tree.selection_set(str(self.index))
            return
        self.refresh(index)

    def add(self):
        if not self.save_current():
            return
        # Place the new default circle beyond the existing mesh extent.
        right = max(float(mesh_section(section).coordinates[:, 0].max()) for section in self.sections)
        self.sections.append(Section(1, 1, 8, (5., 1., right + 15., 0., 0.)))
        self.refresh(len(self.sections) - 1)

    def remove(self):
        if len(self.sections) == 1:
            messagebox.showerror("Geometry error", "Keep at least one boundary.", parent=self.window)
            return
        self.sections.pop(self.index)
        self.refresh(min(self.index, len(self.sections) - 1))

    def apply(self):
        if not self.save_current():
            return
        try:
            deck = geometry_deck(self.deck, tuple(self.sections))
            self.on_apply(deck)
        except ValueError as error:
            messagebox.showerror("Geometry error", str(error), parent=self.window)
            return
        self.window.destroy()
