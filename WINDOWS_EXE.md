# Standalone Windows application

Double-click `BEMUnderground.exe`. Python and the numerical/GUI libraries are
bundled; no Python or Conda installation is required. The executable targets
64-bit Windows and may take a few seconds to unpack at startup.

On first launch, example DAT files are copied to
`%LOCALAPPDATA%\BEMUnderground\examples\in`. Results are saved by default under
`%LOCALAPPDATA%\BEMUnderground\examples\results`. Existing inputs are not
overwritten on subsequent launches. Custom DAT files can be opened using Browse.
The thesis PDFs and original Fortran archive remain in the GitHub repository.

The executable is not digitally signed. Verify that your copy came from the
project owner before running it. Numerical validation limitations documented
in the source repository also apply to this build.

## Rebuild from source

On Windows, create `.venv` using INSTALLATION.md, then run `BUILD_EXE.cmd`.
The output is `dist/BEMUnderground.exe`; the build directories are ignored by Git.
The build uses PyInstaller 6.22.3. Python is required to build, but not to run
this executable. Run the build on Windows for a Windows executable.

A packaged smoke check can be run with `BEMUnderground.exe --self-test REPORT`,
where REPORT is an absolute JSON path in an existing writable directory. It
creates a GUI without displaying it, discovers examples, runs the In-2.1 solver
in a spawned process, and writes NPZ/PDF test outputs beside the report.

## Build verification

The Windows x64 build passed GUI initialization, discovery of 59 examples,
a 20-step In-2.1 solve in a spawned worker, finite displacement checks, and
PDF export. The executable was launched outside the project directory with
Python removed from PATH. A separate check verified that existing user inputs
are preserved. Full interactive and numerical regression testing was not run.
