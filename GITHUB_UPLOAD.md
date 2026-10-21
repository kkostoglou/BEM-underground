# Upload to GitHub

This folder is the repository root. Generated results and local environments
are ignored; original source, DAT inputs, and documentation are included.
See BUNDLE_CONTENTS.md for the exact scope. Keep this checkout and use an
editable installation (`pip install -e .`), as the GUI discovers examples
relative to the source tree.

Create an empty repository in your GitHub account without an initial README,
license, or .gitignore. From PowerShell in this project folder, review and commit:

```powershell
git status --short
git add .
git diff --cached --stat
git commit -m "Initial import of BEM Underground"
```

Then replace YOUR-ACCOUNT and YOUR-REPOSITORY with your actual values:

```powershell
git remote add origin https://github.com/YOUR-ACCOUNT/YOUR-REPOSITORY.git
git push -u origin main
```

Git may ask you to configure your name and email before committing and to
sign in to GitHub before pushing. Use your own identity and account.

No reuse license has been selected. Choose distribution terms separately for
the code and preserved historical documents if you want to grant reuse rights.
The smoke-check workflow checks installation, imports, and CLI startup; it does
not establish numerical agreement with the thesis or validate the GUI visually.

## Git ownership warning

If Git reports "detected dubious ownership" for this folder, use a per-command
exception for this trusted project. This avoids changing global Git settings:

```powershell
git -c safe.directory="$((Get-Location).Path.Replace('\', '/'))" status --short
```

Apply the same `-c safe.directory=...` option immediately after `git` in each
commit, remote, and push command above if the warning persists.

## Preparation checks

Local checks passed on Windows with Python 3.14: editable installation, source
compilation, NumPy/SciPy/Matplotlib/Tk imports, source-root discovery, parsing
`In-5.27-LD.DAT`, CLI help, and `pip check`. Git exclusions and the main guides'
local links were checked. A scan for common credential patterns found no matches.
The GitHub workflow targets Python 3.13 and has not yet run on GitHub.
Full numerical regression and interactive GUI testing were not performed.
