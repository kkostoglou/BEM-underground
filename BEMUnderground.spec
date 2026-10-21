# Build with: python -m PyInstaller --clean --noconfirm BEMUnderground.spec
from pathlib import Path
root = Path(SPECPATH)
a = Analysis(
    [str(root / "packaging/desktop_entry.py")],
    pathex=[str(root / "src")],
    binaries=[],
    datas=[(str(root / "examples/in"), "examples/in")],
    hiddenimports=["matplotlib.backends.backend_pdf"],
    hookspath=[], hooksconfig={"matplotlib": {"backends": ["TkAgg"]}},
    runtime_hooks=[], excludes=[], noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, a.binaries, a.datas, [],
          name="BEMUnderground", debug=False, bootloader_ignore_signals=False,
          strip=False, upx=False, console=False)
