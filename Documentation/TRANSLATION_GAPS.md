# Translation scope and remaining gaps

The Python translation is validated for the curated single-zone, multi-zone, plane-wave, point-explosion, and PhD figure workflows. The following original-program capabilities remain outside the complete end-to-end translation.

## Priority 1: end-to-end workflow gaps

- Direct-time simulation with enclosing-boundary equations in `workflow.py`.
- Reflected explosion loading in Laplace and time workflows.
- Multi-zone load types other than the currently supported plane-wave cases.
- Simulation and enclosing-boundary matrix equations in `assembly.py`.

These require new numerical workflows and DAT acceptance cases.

## Priority 2: numerical and input controls

- BEMZ custom quadrature controls (`NGAUS`/`NSNG`) in the multi-zone kernels.
- Remaining legacy load and trailing-module records rejected by `input_deck.py`.
- Exact equivalence for every historical deck, including output formatting and all optional records.

## Priority 3: GUI output quantities

The GUI labels the following quantities as unavailable:

- Cartesian stress X, Y, and XY;
- Traction X and Y.

The underlying displacement, hoop-stress, boundary-stress, and supported traction quantities remain available.

## Validation policy

A feature is considered translated when it has:

1. an input parser path;
2. a solver/workflow implementation;
3. a curated DAT example;
4. a deterministic expected-result check; and
5. a corresponding automated test.

The following command and result record an earlier checkout. The `tests-final`
suite is not included in this source bundle, so this historical result is not a
passing test claim for the bundle:

```powershell
python -m pytest tests-final -p no:cacheprovider
```

Current result: 217 passed and one opt-in slow validation skipped. The skipped test is enabled with `BEU_RUN_SLOW_VALIDATION=1`.
