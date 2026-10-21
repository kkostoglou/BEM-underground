from argparse import ArgumentParser
from pathlib import Path

from .input_deck import read_input_deck
from .workflow import run_laplace_case, run_time_case, save_result


def main() -> None:
    parser = ArgumentParser(description="Run the translated 2-D Boundary Element solver")
    parser.add_argument("input", type=Path, help="PREBEM keyword .DAT file")
    parser.add_argument("--output", type=Path, default=Path("examples/results/npz/python_bem_result.npz"))
    arguments = parser.parse_args()
    deck = read_input_deck(arguments.input)
    if deck.domain in (1, 4):
        result = run_laplace_case(deck)
    elif deck.domain in (2, 5):
        result = run_time_case(deck)
    else:
        raise SystemExit(f"DOMAIN {deck.domain} is not yet supported end-to-end")
    save_result(result, arguments.output)
    print(f"Saved {len(result.time)} time steps to {arguments.output}")


if __name__ == "__main__":
    main()
