# Chapter 5 multi-zone validation

This note records the validation cases transcribed from `Documentation/parts/Chap5_MultiZone-Doc_KostasPhD-3.pdf`.

## Cases executed

### §5.3.1.1 Equal-material lined circular cavity

- Inner radius: 2.50 m; outer radius: 2.70 m
- Two boundaries with 16 quadratic elements each (64 boundary nodes)
- E = 20,000 MPa, nu = 0.25, density = 2.55 Mg/m3
- Plane compressional wave: -10 MPa, incidence 0 degrees
- 20 steps over 0.015 s

Input: [ch5_3_1_1_equal_liner.dat](../examples/in/Prev1992/ch5_3_1_1_equal_liner.dat)  
Result: [ch5_3_1_1_equal_liner.npz](../results/ch5_3_1_1_equal_liner.npz)

Observed: 64 coordinates, 20 time values, and finite displacement and traction histories. The maximum absolute displacement is 0.0025804 and the maximum absolute stored traction is 19.5120 (solver units).

### §5.3.1.2 Stiff liner in stiff medium

- Inner radius: 5.18 m; outer radius: 5.38 m
- Two boundaries with 16 quadratic elements each (64 boundary nodes)
- Surrounding medium: E = 62,055 MPa, nu = 0.25, density = 2.67 Mg/m3
- Liner: E = 206,750 MPa, nu = 0.25, density = 7.50 Mg/m3
- Plane compressional wave: -10 MPa, incidence 0 degrees
- 20 steps over 0.015 s

Input: [ch5_3_1_2_stiff_liner.dat](../examples/in/Prev1992/ch5_3_1_2_stiff_liner.dat)  
Result: [ch5_3_1_2_stiff_liner.npz](../results/ch5_3_1_2_stiff_liner.npz)

Observed: 64 coordinates, 20 time values, and finite displacement and traction histories. The maximum absolute displacement is 0.0016062 and the maximum absolute stored traction is 19.0186 (solver units).

## Scope and remaining comparison work

These executions validate input parsing, multi-zone geometry/material assignment, time stepping, and finite numerical histories. They do not yet constitute a point-for-point reproduction of the thesis stress plots. The current `MultiZoneResult` stores displacement and traction, but does not expose the reconstructed boundary stress/hoop-stress quantity used by the thesis figures. Therefore the reported thesis peaks (for example, the §5.3.1.4 explosion peaks of -108.3 MPa at node 37, step 7 and -76.9 MPa at node 41, step 10) cannot be compared automatically yet.

The next implementation required for graph-level validation is multi-zone boundary stress post-processing (including tangential/hoop stress and the thesis node/step conventions), followed by running the explosion and additional liner cases in §5.3.1.3–§5.3.1.6.

## Additional deck check

The generalized `examples/multizone_thesis_example.dat` deck was also executed. It currently stops during coupled-system assembly with `zone membership does not produce a square coupled system`. This identifies an input/deck-definition issue in that synthetic five-boundary example; it is not counted as a passing validation case until the boundary-to-zone mapping is corrected.
