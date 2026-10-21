# Multi-zone translation reference

The thesis Chapter 5 defines the multi-zone method for both Laplace and time
domains. Each zone is assembled independently, then common interface
conditions are added:

* displacement compatibility: the displacement on both sides of an interface
  is equal;
* traction equilibrium: interface tractions are equal and opposite after the
  boundary orientation is accounted for;
* free-boundary conditions are applied to the external boundary, while the
  interface jump term is handled separately for each zone.

The resulting global system contains the zone boundary operators and the
non-zone boundary operators. In the original implementation these stages are
represented by `GENRF`, `GENRF5`, `GENRGL`, `GENRG3`, and `GENRGR1`, followed by
`BOUND` and `SOLVE`.

The Python implementation currently provides the parsed zone records,
zone-specific materials, deterministic block layout, boundary node/element
ranges, and interface grouping. The next numerical stage must use these ranges
to place each zone's `[F]` and `[G]` blocks and then add compatibility and
equilibrium rows. The same block construction is reusable for time-domain
convolution; only the per-step kernels and history right-hand side change.

No multi-zone system is enabled through `run_laplace_case` or `run_time_case`
until those interface rows are assembled, preventing a zone deck from being
silently solved as a single homogeneous material.
