# Install and run BEM Underground

This guide starts from a fresh download. The **Conda route** is recommended on Windows because it installs Python 3.13 and the numerical libraries together. A standard Python virtual environment is also supported. Use one route, not both. The application needs a desktop session for its graphical interface; command-line DAT runs do not need a display.

## 1. Get the project

Clone the GitHub repository or download its source ZIP and extract it. Open PowerShell **in the directory containing** `RUN_BEM_GUI.cmd`, `environment.yml`, `examples`, and `src`. For example:

```powershell
cd 'C:\path\to\BEM-underground'
Test-Path .\RUN_BEM_GUI.cmd
Test-Path .\examples\in\In-5.27-LD.DAT
```

Both commands should print `True`. Keep the project folder in place after creating an environment: Python environments may contain absolute paths and should be recreated if the folder is moved to another computer or location.

## 2A. Windows: Conda installation (recommended)

1. Install [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) or an existing Anaconda installation. Open an **Anaconda Prompt**, or a PowerShell session where `conda --version` works.
2. In the project directory, create the environment:

   ```powershell
   conda env create --prefix .\.conda-env --file environment.yml
   ```

   This reads `environment.yml` and installs Python 3.13, NumPy, SciPy, Matplotlib, and pytest into the local `.conda-env` folder. It needs internet access the first time. The folder is excluded from GitHub because each visitor creates it for their own system.
3. Verify the interpreter and GUI libraries:

   ```powershell
   & '.\.conda-env\python.exe' --version
   & '.\.conda-env\python.exe' -c "import numpy, scipy, matplotlib, tkinter; print('Environment OK')"
   ```

4. Start the GUI:

   ```powershell
   .\RUN_BEM_GUI.cmd
   ```

The launcher finds `.conda-env\python.exe`, sets `PYTHONPATH` to `src`, and opens the graphical application. You do not need to activate Conda before using this launcher.

## 2B. Windows: standard Python virtual environment

Install Python **3.13 or newer** with Tk support from [python.org](https://www.python.org/downloads/), then verify `py -3.13 --version` in PowerShell. In the project directory:

```powershell
py -3.13 -m venv .venv
& '.\.venv\Scripts\python.exe' -m pip install --upgrade pip
& '.\.venv\Scripts\python.exe' -m pip install -e .
& '.\.venv\Scripts\python.exe' -c "import numpy, scipy, matplotlib, tkinter; print('Environment OK')"
.\RUN_BEM_GUI.cmd
```

`pip install -e .` installs the package and its declared numerical dependencies. The launcher uses `.venv` when `.conda-env` is absent; if both exist, it chooses `.conda-env`. Directly invoking `.venv\Scripts\python.exe` avoids PowerShell activation-policy changes.

## 2C. Linux or macOS

The Windows CMD launcher and preserved Fortran executables are Windows-specific. The Python command-line solver can be started with a Python 3.13+ installation; the GUI additionally requires a desktop session and a working Tk installation. This route has not been validated as extensively as Windows.

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e .
PYTHONPATH=src .venv/bin/python -m python.gui
```

If `import tkinter` fails, install the Tk package supplied by your operating system or Python distribution, then repeat the verification command. On a headless server, use the command-line solver instead of the GUI.

## 3. Run an example

In the GUI, choose **Category → Sub-Category → Example**, or use **Browse DAT...** to select a file from `examples/in`. Review the problem's material values, geometry, load, DOMAIN, time steps, and plot nodes, then press **RUN**. The GUI saves an NPZ archive under `examples/results/npz` by default. **Open existing result...** loads a saved archive without recalculating it. PDF reports can be exported from the results view.

To run a DAT file from PowerShell without opening the GUI, use the interpreter you installed:

```powershell
$env:PYTHONPATH = (Resolve-Path '.\src').Path
& '.\.conda-env\python.exe' -m python '.\examples\in\In-5.27-LD.DAT' --output '.\examples\results\npz\In-5.27-LD_result.npz'
```

For a `.venv` setup, substitute `.\.venv\Scripts\python.exe`. The command-line solver requires the output directory to exist; on a fresh download run `New-Item -ItemType Directory -Force .\examples\results\npz` first. The GUI creates its output directory automatically.

DAT names ending in `-LD`, `-TD`, and `-CQ` normally select Laplace-domain, direct time-domain, and convolution-quadrature variants. The `DOMAIN` record in the DAT file is authoritative: DOMAIN 1 is Laplace, DOMAIN 2 is direct TD, and DOMAIN 5 is CQ. Some historical or thesis-example decks use methods not fully translated to Python. A solver error for one such case is distinct from an installation failure; consult `Documentation/TRANSLATION_GAPS.md` and the status notes under `Documentation/`.

## 4. Check or repair an installation

Run these commands from the project directory. Replace the interpreter path with `.\.venv\Scripts\python.exe` if you chose the venv route.

```powershell
& '.\.conda-env\python.exe' --version
& '.\.conda-env\python.exe' -c "import numpy, scipy, matplotlib, tkinter; print('Imports OK')"
$env:PYTHONPATH = (Resolve-Path '.\src').Path
& '.\.conda-env\python.exe' -c "from python.case_catalog import project_root; print(project_root())"
```

The final command should print this project's directory. If the launcher says no environment was found, check that `.conda-env\python.exe` or `.venv\Scripts\python.exe` exists beside `RUN_BEM_GUI.cmd`. If Python reports `No module named python`, run from the project directory and set `PYTHONPATH` to `src` as shown above. If `tkinter` fails, use a Python distribution with Tk support. If NumPy, SciPy, or Matplotlib fails, repeat the installation command for the selected environment. If you moved the project after creating `.conda-env` or `.venv`, recreate that environment in the new location rather than copying its files.

The GitHub source bundle omits generated NPZ and PDF results. The historical validation scripts mentioned in the status notes are not included in this checkout. Start with `USER_GUIDE.md` for DAT editing and GUI controls, and `examples/README.md` for the input and result layout.

The environment commands follow the official [Conda environment-file reference](https://docs.conda.io/projects/conda/en/stable/commands/env/create.html), [Python `venv` guide](https://docs.python.org/3.13/library/venv.html), and [pip local-project installation guide](https://pip.pypa.io/en/stable/topics/local-project-installs/).
