# Multi-zone validation status

## In-5.3 lining orientation

The LD and TD In-5.3 decks each define two concentric boundaries at radii
2.70 m and 2.50 m, with an exterior zone (`BOU1=0, BOU2=1`) and a finite
liner zone (`BOU1=1, BOU2=-2`). The mesh starts each circle clockwise.
`prepare_zone_connectivity` retains that order for exterior interface c and
liner inner surface b, and reverses the outer liner interface a. At the
rightmost point, the resulting material normals are c = (-1, 0),
a = (+1, 0), and b = (-1, 0). This is checked for both decks by
`validation/check_lining_orientation.py`.

Consequently, the implementation follows b clockwise and a counter-clockwise,
but c clockwise, **not** the counter-clockwise direction stated in one
paragraph of the thesis. The c arrow in the supplied Figure 5.4 and the
separate section 5.3.1.1 description indicate clockwise traversal. The
original BEMZ `CSHAPE` normal is the right-hand normal `(dy, -dx)`, so c
clockwise points outward from the exterior material into the liner, while a
counter-clockwise points outward from the liner into the exterior material.
These opposite normals are required at their common interface. No orientation
change was made on the basis of the conflicting sentence alone.

The translated multi-zone implementation is checked against the retained
thesis formulation and the November 1993 source routines (`BEMZ1` through
`BEMZ4`, `GENRF`, `GENRGL`, and `GENRGR1`). The repository contains matrix,
connectivity, kernel, Laplace workflow, and direct-time workflow tests.

An independent numerical comparison is now available for the Laplace-domain
In-5.3 lined cavity. A local GNU Fortran build of the November 1993
`BEMZ1.FOR`--`BEMZ4.FOR` sources solves the full 64-node, 20-step case.
The original source is untouched; build copies only translate Microsoft
allocatable syntax, remove the obsolete `BLOCKSIZE` argument, and allocate
one-element dummy arrays when there are no internal points. The numeric
`BEMLTDAT`, `BEMLTTRA`, and `BEMLTSTR` inputs are produced by
`validation/export_bemz_reference.py`, because the available older PREBEM
rejects the `ZONES` card. The 1993 tree retains `PREWRITE.FOR` and
`PREWRIT1.FOR`, but not a complete matching PREBEM main program and parser.
Thus the comparison independently validates BEMZ's
LD coupled solver with identical supplied incident loads, not the 1993
preprocessor or incident-load generation.

`validation/compare_in53_bemz.py` compares all 20 complex displacement
spectra. At the thesis-stated 2.50/2.70 m radii, the maximum absolute
difference is 4.780e-12 m, or 1.218e-6 of the reference peak. At the separate
5.00/5.20 m figure-comparison radii, it is 4.857e-12 m, or 8.758e-7 of the
peak. These differences are within the six-significant-figure precision of
the Fortran `.PRO` text output. The comparison does not validate hoop-stress
recovery: the generated Fortran incident-stress records are zero placeholders
because only displacement was compared. No original multizone direct-time
solver has been found in this source set, so the Python TD lined-cavity result
remains unvalidated against multizone Fortran.

The checked build and reference outputs are under
`validation/reference_build/bemz` and `validation/reference_runs/in53ld_bemz_*`.
On the validation Windows installation, rebuild with
`powershell.exe -NoProfile -ExecutionPolicy Bypass -File validation/build_bemz_reference.ps1`,
then rerun `validation/compare_in53_bemz.py` with `PYTHONPATH=src`. The
`export_bemz_reference.py` manifest in each run directory records the deck
and mesh used for its Fortran input.

## Check of the Disk2 BEMZONL.F source

`Disk2.ExeAndSource/SOURCE/MULTIZONE/BEMZONL.F` is an earlier, distinct
Laplace-domain solver (SHA-256
`EC9F4EF6CB48186F1642D32F065EE833D2F5092598E45921C8AA4D1EFEB932DB`).
It is not the November 1993 `BEMZ1.FOR`--`BEMZ4.FOR` program. Its numeric
`BEMLTDAT` format has a single `COMMON/MATER` material block and an `nffs`
simulation-surface switch, rather than the later solver's per-zone material
and `KZON`/`NZTYPE` records. It cannot check the general multi-material
implementation or the Python direct-time implementation.

An isolated GNU Fortran build of the Disk2 source is at
`validation/reference_build/bemzonl_1991/BEMZONL.F`. The untouched source
does not compile with modern GNU Fortran because one `IMPLICIT` declaration
follows a `COMMON` declaration and it references the unavailable Unix
`ieee_handler` routine. The validation copy reorders those declarations and
omits the floating-point trap. With bounds checking, the first run then
found an uninitialized `i` used in `konn(jcl,i)` in the displacement-history
storage loop (original lines 1608--1610). The validation copy adds an
explicit `i=1..3` loop there; it does not alter the kernel or equation
assembly routines. The source file in `Disk2.ExeAndSource` remains unchanged.

`validation/export_bemzonl_1991.py` reconstructs the older numeric input
for In-5.3 with the outer circle as a type-3 simulation surface and the
inner circle as a type-2 free surface. The original matching 1991 PREBEM
output for this case was unavailable, so these node/element codes are
inferred from source comments and conditions. Incident tractions still
come from the Python load preprocessor. The patched program completes a
64-node, 20-step run under bounds checking at
`validation/reference_runs/in53ld_bemzonl_1991_full/IN53.PRO`, but the
result **does not match** the Python/1993 result: maximum difference
`1.052239175e-4 m` versus a Python peak of `3.925073045e-6 m`.
Its final second-stage record gives zero `u_x` at inner node 33, so this
input reconstruction cannot be accepted as a numerical reference for the
lined cavity. `validation/compare_in53_bemzonl_1991.py` prints the
diagnostics. The match to the later 1993 BEMZ solver remains the stronger
independent LD check; the TD lined-cavity discrepancy and comparison with
the thesis figures remain unresolved.

As a separate control, the Python single-zone TD `TINFPL1.DAT` result was
compared with the retained original BEMTIM Fortran `.PRO` output. All 20 steps
and 32 displacement components agree to a maximum 5.072e-9 m, or 1.003e-6
of the reference peak (again limited by printed precision). This validates
the single-zone TD path, but it does not establish correctness of the newer
Python coupled-zone TD extension used for In-5.3.

The following reference locations were searched for a retained historical
multi-zone result or executable:

- `G:\KostasMy\PhD\Dynamics.FinalDisks\FinalNov93`
- `G:\KostasMy\PhD\prev\Bak3.10.93`
- `G:\KostasMy\PhD\prev\bak93\Bak3.10.93`
- `G:\KostasMy\PhD\other\iFortran\i1\test1\test1`
- `D:\KostasMy\PhD\Bak3.10.93`
- `G:\KostasMy\PhD\other\Conferences\susi`
- `G:\KostasMy\PhD\other\Conferences\Conf9.1.92`
- the archived `Disk2.ExeAndSource` tree in this project

Those locations contain Fortran source, make files, executables, input decks,
and postprocessor programs, but no saved `BEMZ1`/`BEMZONL` displacement or
traction result. The original multi-zone executable is also not present in the
retained FinalNov93 executable directory; the current comparison uses a local
build from its preserved source.

The `susi` directory contains thesis and conference working documents
(`THESIS.WP`, `BEMLAPL.THE`, `REFLECT.THE`, and related files). They confirm the
Laplace-transform/BEM formulation, multi-boundary capability, and Durbin/FFT
inversion described by the implementation, but they contain no multi-zone
numeric output or `BEMZONL` run files.

`Conf9.1.92` contains an earlier standalone multi-zone source set:
`BEMZONL.F` (SHA-256
`CB9BDB51852528D7C2E52722BED18749D58F1B48E717CFA6740C7D23F5C73340`),
`PREBEM.F`, `ELASTOST.F`, and `POSTBEM.F`. It has no DAT input, scratch
files, displacement history, traction history, or executable result. Its
`BEMZONL.F` is a distinct 1991 source revision from the archived 1993 source
(`EC9F4EF6CB48186F1642D32F065EE833D2F5092598E45921C8AA4D1EFEB932DB`), so it
is useful for historical source comparison but cannot serve as a numerical
reference output.

The `G:\KostasMy\PhD\other\iFortran\i1\test1\test1` directory contains the
1993 multi-zone source tree and legacy make files (`BEMZ1.MAK`, `BEMZ.MAK`, and
`POSTBEMZ.MAK`). Its `BEMZ1.FOR` and `BEMZ3.FOR` hashes exactly match the
preserved originals in `src/fortran-original/Disk2.ExeAndSource/SOURCE/other/Bak3.10.93`, confirming the source
provenance. It contains no compiled `BEMZ1.EXE`, input deck, scratch file, or
numeric result output.

Available reproducible checks are:

- `examples/multizone_lined_cavity.dat` — 20-step Laplace run;
- `examples/multizone_lined_cavity_time.dat` — coarse six-step direct-time run;
- finite, nonzero displacement and traction histories from both runs;
- automated tests covering the translated operators and workflows.

The direct-time multi-zone example is a Python extension of the thesis
multi-zone formulation. It still requires comparison with an independent
multizone direct-time implementation or an original TD result; the rebuilt
`BEMZ1` program solves only in the Laplace domain.

The Unix examples directory
`G:\KostasMy\PhD\prev\Bak15.8.92.Unix.Examples` does contain historical
numerical artifacts, but they are from the ordinary BEMLAPL/BEMTIM workflows,
not the `BEMZONL` program. In particular:

- `TSEMEXP1.DAT` describes a 20-step direct-time semi-infinite, two-boundary
  explosion case and is paired with `TSEMEXP1.CEL` (229,080 bytes);
- `TIFPLSV1.DAT` is paired with `TIFPLSV1.CEL` (117,248 bytes);
- `TINFPLC2.DAT` and `TINFPLC3.DAT` are paired with `.CEL` files of 103,360
  and 106,424 bytes.

The `.CEL` files are legacy binary result/postprocessor files. They are useful
for future validation of ordinary multi-boundary simulation cases, but they do
not contain a `ZONES`/`BEMZONL` run and are not directly comparable with the
translated coupled-zone result arrays without first decoding the CEL format.
