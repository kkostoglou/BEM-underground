"""Persistent user data for frozen builds; source checkouts retain their layout."""
from pathlib import Path
import os
import shutil
import sys


def project_directory() -> Path:
    if getattr(sys, "frozen", False):
        return Path(os.environ.get("BEM_USER_HOME") or
                    Path(os.environ.get("LOCALAPPDATA", Path.home())) / "BEMUnderground")
    return Path(__file__).resolve().parents[2]


def prepare_desktop_data() -> None:
    if not getattr(sys, "frozen", False):
        return
    destination = project_directory()
    source = Path(sys._MEIPASS) / "examples" / "in"
    for original in source.rglob("*"):
        if original.is_file():
            target = destination / "examples" / "in" / original.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                shutil.copy2(original, target)
    for kind in ("npz", "pdf"):
        (destination / "examples" / "results" / kind).mkdir(parents=True, exist_ok=True)
    os.chdir(destination)
