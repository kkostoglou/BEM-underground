# Example inputs and results

`in/` contains the maintained DAT decks, including thesis figure inputs and
step-size variants. `in/Prev1992/` preserves historical decks.
`in/manifest.toml` describes a selected subset, not the complete inventory.

`reference/` contains the thesis-figure index PDF and Appendix F DAT checking
template and guide. The complete thesis is in `../Documentation/Doc_KostasPhD.pdf`.

Generated solver archives in `results/npz/` and reports in `results/pdf/` are
excluded from Git. Existing local results are retained on disk. After cloning,
run a DAT case to create an archive; the GUI can reopen it without solving again.
See `../INSTALLATION.md` for commands and `../USER_GUIDE.md` for GUI controls.

The historical automated comparison scripts and reference-data directory are
not included. The notes under `../Documentation/` record earlier validation
work and its limitations; references there to missing scripts describe that
historical working project.
