"""Isolated GUI solver process so STOP can terminate native numerical work."""

from __future__ import annotations

import time
import traceback
from pathlib import Path

from .workflow import run_laplace_case, run_time_case, save_result


def run_gui_case(deck, temporary_path: Path, error_path: Path, progress_queue) -> None:
    last_update = 0.0
    last_stage = ""

    def progress(stage: str) -> None:
        nonlocal last_update, last_stage
        now = time.monotonic()
        # Kernel callbacks may occur thousands of times. Keep the newest
        # visible position without turning progress reporting into I/O work.
        category = stage.split(":", 1)[0]
        if category != last_stage or now - last_update >= 0.15:
            progress_queue.put(stage)
            last_update = now
            last_stage = category

    try:
        if deck.domain in (1, 4):
            result = run_laplace_case(deck, progress=progress)
        elif deck.domain in (2, 5):
            result = run_time_case(deck, progress=progress)
        else:
            raise ValueError(f"DOMAIN {deck.domain} is not supported end-to-end")
        progress_queue.put("Saving result")
        temporary_path.parent.mkdir(parents=True, exist_ok=True)
        save_result(result, temporary_path)
    except BaseException:
        error_path.write_text(traceback.format_exc(), encoding="utf-8")
        raise
