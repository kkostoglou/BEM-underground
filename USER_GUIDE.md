# Boundary Element Underground — User Guide

Boundary Element Underground is a Python translation of the historical 2-D boundary-element program. It provides a graphical input form, command-line DAT execution, and figure-numbered DAT examples. The former automated test scripts are absent from this checkout.

## Requirements

- Windows 10 or later
- Python 3.13 or newer; `environment.yml` specifies Python 3.13
- A working desktop session for the Tk/Matplotlib GUI

For a fresh installation, follow [INSTALLATION.md](INSTALLATION.md). The launcher uses `.conda-env` when present, or `.venv` otherwise. To create the recommended Conda environment, open PowerShell in the project directory and run:

```powershell
conda env create --prefix .\.conda-env --file environment.yml
& '.\.conda-env\python.exe' --version
& '.\.conda-env\python.exe' -c "import numpy, scipy, matplotlib, tkinter; print('Python environment OK')"
```

If Conda is unavailable, install Python 3.13 and use a virtual environment. The same CMD launcher will detect it:

```powershell
py -3.13 -m venv .venv
& '.\.venv\Scripts\python.exe' -m pip install -e .
.\RUN_BEM_GUI.cmd
```

## Project layout

```text
src/python/                  Python package and GUI
src/fortran-original/        Preserved original Fortran disks
examples/in              Runnable DAT examples; 1992 decks are in Prev1992
examples/reference/          Thesis-figure index and Appendix F guidance
examples/results/npz         Saved Python solver results
examples/results/pdf         Saved example PDF reports
Documentation/               Validation notes and translation scope
```

The former `examples/in/py` test and helper-script directory is absent from this checkout. The runtime application does not depend on it, but those checks must be restored from a known backup before they can be rerun.

## Start the GUI

From PowerShell in the repository root:

```powershell
.\RUN_BEM_GUI.cmd
```

In the form:

1. The circular-cavity plane-wave case opens with editable defaults and RUN enabled. To change cases, select **Category**, **Sub-Category**, and **Example**, or click **Browse DAT...**. There is no separate “Curated example” control.
2. Choose a file from `examples\in` for maintained examples.
3. The same folder contains the figure-numbered PhD cases.
4. Review or edit the input values shown in the input cells.
5. Click **RUN**.
   The status line beneath the buttons shows elapsed time and the current
   calculation stage. Direct-time (DOMAIN 2) and convolution-quadrature
   (DOMAIN 5) runs use different numerical methods. The moving
   bar indicates that the calculation is still active.
   **STOP** terminates the solver process promptly. A cancelled run discards
   its temporary output and leaves any previously saved result unchanged.
6. Use **Open existing result...** to load an NPZ result from `examples\results\npz`.
7. Use **Export CSV...** to export displayed result data.

The Browse dialog remembers the last input folder used by the application.

Use **Edit boundaries…** beside the boundary count to add or remove boundaries,
or change a selected boundary's shape, dimensions, position, rotation, role, and
element counts. Circles/ellipses, squares, rectangles, and straight lines are
supported. Select **Apply geometry** to use the changes, or **Cancel** to discard
them. The count is calculated from the boundary list. With multiple boundaries,
edit dimensions in this dialog; the single-boundary size fields are inactive.

Geometry edits preserve the current material and loading values. **RUN** uses
the edited geometry, and **Reset inputs** restores the original case. Changes
remain in the current session and do not overwrite the source DAT. Changing
geometry clears node selections and thesis graph presets associated with the
original mesh. Compound segments and cases with zone, simulation, or explicit
node assignments must be edited together with their assignments in the DAT file.

Square and rectangular meshes use shared corner nodes when the DAT omits the
corner setting, matching the PREBEM default (`DIS=1`). Set `DIS=2` to insert
four duplicate corner nodes so each adjoining side retains its own stress
value. The geometry editor also defaults to shared corners. Rerun a case after
changing this setting; previously saved NPZ files retain their old mesh.

## Run a DAT file from PowerShell

Use the project Python environment from the `BEM-underground` directory:

```powershell
$env:PYTHONPATH = 'src'
& '.\.conda-env\python.exe' -m python `
  examples\in\In-5.27-LD.DAT `
  --output examples\results\npz\In-5.27-LD_check.npz
```

A result NPZ contains time values, mesh coordinates and connectivity, displacements, tractions, and available stress results.

## DAT examples

Additional notes about selected DAT cases are in:

`examples\reference\test_data_notes\curated\README.md`

The machine-readable list of selected cases is:

`examples\in\manifest.toml`

The GUI reads DAT files from `examples\in`; it is not limited to `manifest.toml`. The manifest's former automated validator was in `examples\in\py` and is not present in this checkout.

Important examples include:

- `In-5.27-LD.DAT` and `In-5.27-TD.DAT` — Chapter 5 lined half-space cavity.
- `In-5.31-LD.DAT` and `In-5.31-TD.DAT` — Chapter 5 two steel-lined cavities.
- `In-5.13-CQ.DAT` — separate convolution-quadrature comparison.
- `ch6_8_1_1_halfspace_plane_wave.dat` — Chapter 6 P-wave case.
- `THESIS_RECT_EXPLOSION.DAT` in `Prev1992` — PhD Figures 2.40–2.43.

Use the historical directory only when reproducing original 1992 input decks.

`DOMAIN=2` is the original direct time-domain formulation for both single-zone
and `ZONES` decks. It assembles the translated BEMTIM time-lag kernels and is
the default TD case.

`DOMAIN=5` selects BDF2 convolution quadrature and requires `ZONES`. It
samples the coupled multizone Laplace boundary matrices at the BDF2 complex
frequencies, inverse-transforms those matrices into time-domain influence
weights, and then solves one time step at a time. With displacement and
interface traction in the unknown vector, each step solves
`A[0] x[n] = sum(B[j] p[n-j]) - sum(A[j] x[n-j], j>=1)`.
The zone equations enforce displacement continuity and opposite interface
tractions during every step. The CQ method supports the current elastic
multi-zone plane-wave and explosion loading.
It does not implement the quasi-static Perzyna yield law or initial-stress
domain cells from the separate viscoplastic paper.

## Tests and validation scripts

The former `examples\in\py` directory is absent, so the `pyproject.toml` pytest path cannot currently be used. Restore those scripts from a known project backup before running the old automated and thesis-comparison checks. The interpreter and package-import checks in **Requirements** still verify a fresh environment, but they do not validate solver accuracy.

## PhD validation

Use these documents for application-specific comparisons:

- `Documentation\CHAPTER5_MULTIZONE_VALIDATION.md`
- `Documentation\CHAPTER6_WAVE_VALIDATION.md`
- `Documentation\TRANSLATION_GAPS.md`
- `validation\thesis_figure_runs\VALIDATION.md`

The figure-numbered DAT files map directly to the PhD figure numbers, for example `FIG_2_40.DAT` through `FIG_2_43.DAT`.

## Troubleshooting

**The GUI does not start.** Run `RUN_BEM_GUI.cmd` from the repository root and confirm that the Python environment has Tk/Matplotlib installed.

**The DAT file is not listed.** Use **Browse DAT...** and select it from `examples\in`.

**The terminal says NumPy or Matplotlib cannot import.** Use `.conda-env\python.exe` and recreate the environment from `environment.yml` if needed.

**A test path is missing.** The old `examples\in\py` directory has been removed from this checkout; restore it before running pytest.

**A historical file is missing.** Check `examples\in\Prev1992` and the archived disks.

## Current scope

The translated and validated scope includes single-zone Laplace analysis, plane P/SV waves, supported point-explosion cases, multi-zone workflows, time-domain cases, PhD figure mappings, result readers, and GUI input/result handling.

Direct-time enclosing-boundary equations, reflected explosion loading, some specialized legacy records, and several unimplemented GUI output quantities remain documented in `Documentation\TRANSLATION_GAPS.md`.
