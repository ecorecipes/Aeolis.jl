# Contributing to Aeolis.jl

## Local setup

```sh
git clone https://github.com/ecorecipes/Aeolis.jl
cd Aeolis.jl
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

`Pkg.test()` is the only run that counts as "the suite passes".

## Style

- Julia ≥ 1.12. Every dependency, including stdlibs and test-only extras, has a `[compat]` entry.
- Format with JuliaFormatter using the repository `.JuliaFormatter.toml` (`style = "yas"`).
- `snake_case` for functions and variables, `CamelCase` for types, a leading underscore for internal helpers;
  functions that change their arguments end in `!`.
- Public functions have docstrings with a signature line and a one-sentence summary.
- Exceptions are typed, carry the offending names, live in `src/errors.jl` and subtype `AeolisError`; invalid
  arguments and keywords raise `ArgumentError`.

## Agreement with AeoLiS

Aeolis.jl reproduces aeolis-python. A change to a process must keep `test/test_reference.jl` passing. New
reference cases are recorded from aeolis-python: run it on the case, record the state after every step as
float64 (variables with a leading `step` dimension, the Python `(ny+1, nx+1, ...)` order, plus `t` and `dt`),
and add the case directory with its `reference.nc` under `test/data/`.

## Adding a vignette

1. Create `vignettes/NN_short_name/NN_short_name.qmd` with front matter including `engine: julia` and
   `pdf: default` under `format:` (copy an existing vignette's front matter).
2. Add any new dependencies to `vignettes/Project.toml`.
3. `cd vignettes && quarto render && quarto render --to gfm` renders HTML, GitHub-flavoured Markdown and PDF
   (the PDF needs `lualatex` and the JuliaMono font in `../../fonts/JuliaMono/`; see `vignettes/README.md`).
   Commit the `.qmd`, the rendered `.md`, `.html` and `.pdf`, and `NN_short_name_files/`.
4. `julia scripts/sync_vignettes.jl` and commit `docs/src/tutorials/`.

## Pull requests

Branch from `main`, keep PRs focused, and make sure CI is green (tests on Julia 1.12 and latest, formatter,
vignette sync check, docs build).
