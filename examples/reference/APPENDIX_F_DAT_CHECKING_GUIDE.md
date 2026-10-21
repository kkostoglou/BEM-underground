# Appendix F DAT checking guide

Use `APPENDIX_F_CHECK_TEMPLATE.DAT` as a numbered reference when checking input files.

1. Confirm module order: `TITLE`, `TYPE.PROBLEM`, `SECTION.TYPE`, material, steps, time, loading, plots, and `END`.
2. Confirm shape and region codes against the comments in the DAT template.
3. Confirm geometry field order. For squares and rectangles, the final fields are `TYPX`, `TYPY`, and `DIS`.
   For a rectangle, `XSIDE` and `YSIDE` are positive local side lengths; `XC,YC` are the centre; `ROT` is the counter-clockwise angle from global +X to local +X. Generate local corners `(+-XSIDE/2, +-YSIDE/2)`, rotate by `ROT`, then translate by `(XC,YC)`.
   Region 1 cavities use the cavity boundary orientation; region/open-surface codes can reverse traversal and normals.
4. Confirm `DIS=2` when two one-sided corner stresses are required. The four duplicate corner nodes increase the node count by four.
5. Confirm that rectangle elements are per side: `X elements / Y elements`.
6. Confirm loading type and units before comparing results.
   For `EXPLOSION.POINT`, check `XE, YE, W, N, F, DECAY`; `F=1` is the cubic-root formula and `DECAY=0` means no decay in time.
7. Confirm every `NODE.PLOTS` number belongs to the generated mesh; for duplicated corners, choose the side-specific duplicate.
8. Appendix F permits comments only when the first character is `*`; do not place comments inside a numeric block.
