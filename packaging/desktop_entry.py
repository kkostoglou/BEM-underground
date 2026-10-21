"""Frozen desktop entry point with multiprocessing bootstrap support."""
import multiprocessing


def main():
    # Must run before importing GUI/numerical modules in spawned workers.
    multiprocessing.freeze_support()
    import sys
    from python.runtime_paths import prepare_desktop_data
    prepare_desktop_data()
    if len(sys.argv) == 3 and sys.argv[1] == "--self-test":
        from python.desktop_smoke import run_smoke
        run_smoke(sys.argv[2])
    else:
        from python.gui import main as gui_main
        gui_main()


if __name__ == "__main__":
    main()
