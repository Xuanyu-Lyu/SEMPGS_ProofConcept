"""Check that the active environment can run this repository's Python and R code.

Run inside the env:  python check_env.py   (exit code 0 = everything needed is present)
"""
import importlib
import shutil
import subprocess
import sys

ok = True


def report(name, good, detail=""):
    global ok
    ok &= good
    print(f"  [{'ok' if good else 'FAIL'}] {name}{'  ' + detail if detail else ''}")


print(f"Python {sys.version.split()[0]}")
for mod in ["numpy", "pandas", "scipy", "joblib", "matplotlib", "IPython", "nbconvert", "pytest"]:
    try:
        m = importlib.import_module(mod)
        report(mod, True, getattr(m, "__version__", ""))
    except ImportError as e:
        report(mod, False, str(e))

pd_major = int(importlib.import_module("pandas").__version__.split(".")[0])
report("pandas < 3", pd_major < 3)

rscript = shutil.which("Rscript")
report("Rscript on PATH", rscript is not None, rscript or "")
if rscript:
    r_code = r"""
    cat("R", R.version$major, ".", R.version$minor, "\n", sep = "")
    for (p in c("OpenMx", "data.table", "stringr", "MASS", "polycor", "parallel")) {
      v <- tryCatch(format(packageVersion(p)), error = function(e) NA)
      cat("pkg", p, if (is.na(v)) "MISSING" else v, "\n")
    }
    if (requireNamespace("OpenMx", quietly = TRUE)) {
      suppressMessages(library(OpenMx))
      cat("flag NPSOL", imxHasNPSOL(), "\n")
      cat("flag OpenMP", imxHasOpenMP(), "\n")
      mxOption(NULL, "Default optimizer", "NPSOL")
      set.seed(1); d <- data.frame(x = rnorm(200))
      m <- mxModel("t", type = "RAM", manifestVars = "x", mxPath("x", arrows = 2, free = TRUE, values = .5),
                   mxPath("one", "x"), mxData(d, "raw"))
      f <- tryCatch(mxRun(m, silent = TRUE), error = function(e) NULL)
      cat("flag NPSOL_fit", !is.null(f) && f$output$status$code == 0, "\n")
    }
    """
    out = subprocess.run([rscript, "-e", r_code], capture_output=True, text=True).stdout
    for line in out.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0].startswith("R") and len(parts) == 1:
            print(parts[0].replace("R", "R ", 1))
        elif parts[0] == "pkg":
            report(parts[1], parts[2] != "MISSING", parts[2])
        elif parts[0] == "flag":
            report(f"OpenMx {parts[1]}", parts[2] == "TRUE")

print("\nEnvironment OK" if ok else "\nEnvironment INCOMPLETE (see FAIL lines; ./setup_env.sh rebuilds)")
sys.exit(0 if ok else 1)
