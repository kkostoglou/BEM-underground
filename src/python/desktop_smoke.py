"""Explicit packaged-app smoke check, invoked only with --self-test REPORT."""
from pathlib import Path
import json
import multiprocessing
import traceback


def run_smoke(report_path):
    report = Path(report_path).resolve()
    status = {"passed": False}
    try:
        import numpy as np
        from tkinter import Tk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_pdf import PdfPages
        from .case_catalog import project_root
        from .gui import BemApplication
        from .gui_run_worker import run_gui_case
        from .input_deck import read_input_deck
        window = Tk()
        window.withdraw()
        try:
            app = BemApplication(window)
            cases = app.examples
            assert cases, "No bundled examples discovered"
            window.update()
        finally:
            window.destroy()
        output = report.with_suffix(".npz")
        error = report.with_suffix(".error.txt")
        deck = read_input_deck(project_root() / "examples/in/In-2.1.DAT")
        context = multiprocessing.get_context("spawn")
        progress = context.Queue()
        worker = context.Process(target=run_gui_case, args=(deck, output, error, progress))
        worker.start()
        worker.join(120)
        if worker.is_alive():
            worker.terminate()
            worker.join()
            raise RuntimeError("Packaged solver exceeded 120 seconds")
        assert worker.exitcode == 0, error.read_text() if error.exists() else worker.exitcode
        progress.close()
        with np.load(output) as result:
            assert result["time"].size > 0
            assert np.isfinite(result["displacement"]).all()
            status["time_steps"] = int(result["time"].size)
        fig = Figure()
        fig.add_subplot().plot([0, 1], [0, 1])
        with PdfPages(report.with_suffix(".pdf")) as pdf:
            pdf.savefig(fig)
        status.update(passed=True, examples=len(cases))
    except Exception:
        status["error"] = traceback.format_exc()
    report.write_text(json.dumps(status, indent=2), encoding="utf-8")
    if not status["passed"]:
        raise SystemExit(1)
