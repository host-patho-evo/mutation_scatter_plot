# Contributing to mutation_scatter_plot

We welcome pull requests that fix bugs, improve documentation, or add new
features. This document covers the practical aspects of contributing.

## Licensing

This project is licensed under
**[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)**.
By submitting a pull request you agree that your contributions will be
licensed under the same terms. No dual-licensing or CLA is required —
just make sure your contribution is compatible with CC BY 4.0.

Please keep the copyright header present in each source file:

```python
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International.
```

## Getting Started

1. Fork the repository on GitHub.
2. Clone your fork and create a **feature branch** with a short, descriptive
   name (do not work directly on `main`):

   ```bash
   git clone https://github.com/<you>/mutation_scatter_plot.git
   cd mutation_scatter_plot
   git checkout -b fix-colorbar-label
   ```

3. Install in editable mode with development dependencies:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -e .[dev]
   ```

## Coding Conventions

- Follow [PEP 8](https://peps.python.org/pep-0008/) and
  [PEP 257](https://peps.python.org/pep-0257/) for code style and docstrings.
- Prefer **f-strings** for string formatting. Legacy `%` and `.format()`
  exist in some older scripts but new code should use f-strings.
- Use `decimal.Decimal` for frequency values and any arithmetic that must be
  bit-identical across platforms — never bare `float`.
- Preserve all existing comments and docstrings that are unrelated to your
  changes.
- Import order: stdlib → third-party → local (relative imports within the
  package).

### Static Analysis

We enforce clean static analysis on every push via GitHub Actions:

- **pyflakes** — must report zero errors. Unused imports, undefined names,
  and redefined-but-unused variables will fail CI.
- **flake8 + flake8-bugbear** — via pre-commit (see below).
- **pylint** — via pre-commit; target score is 10.0/10.

### Pre-commit Hooks

We use [pre-commit](https://pre-commit.com/) to run linters automatically
before each `git commit`. Install it once:

```bash
pip install pre-commit
pre-commit install
```

The hooks configured in `.pre-commit-config.yaml` include:
- Trailing whitespace and end-of-file fixes
- YAML validation
- flake8 with flake8-bugbear
- pylint

## Testing

**All pull requests must pass the full test suite.** New features and bug
fixes should include tests wherever practical.

### Running Tests

```bash
pip install -e .[dev]
pytest tests/
```

The test suite contains **8 400+** regression tests and takes approximately
5–6 minutes on a modern machine. Please run the full suite before
submitting a pull request.

### Golden File Regression Tests

The test suite uses a "golden file" approach: rendered outputs (PNG, PDF,
HTML, TSV) are compared against committed baselines. If you intentionally
change rendering behaviour, regenerate the baselines:

```bash
REGENERATE_GOLDENS=1 pytest tests/test_mutation_scatter_plot.py
```

Review the regenerated files carefully with `git diff` before committing
them. For details see
[docs/testing_framework.md](docs/testing_framework.md).

### Continuous Integration

Once you submit your pull request, the **Pyflakes** workflow runs
automatically on Python 3.11, 3.12, and 3.13. CI must pass before
your pull request will be merged.

## Pull Request Guidelines

- **One logical change per pull request.** If you are fixing a bug and also
  refactoring nearby code, please split them into separate PRs.
- **Write a clear commit message.** The first line should summarise the
  change in ≤72 characters. Use the body to explain *why* the change is
  needed, not just *what* it does.
- **Reference issues.** If your PR fixes an open issue, include
  `Fixes #123` in the commit message body.
- **Do not mix formatting-only changes with functional changes.**

### Commit Message Convention

```
<type>: <short summary>

<optional body explaining motivation and design decisions>
```

Common types: `fix`, `feat`, `refactor`, `docs`, `test`, `perf`, `ci`.

Example:

```
fix: treat cmap_vmin=None as unset in colormap initialisation

The hasattr check was passing when the attribute existed but was set
to None (e.g. from CLI defaults or test fixtures), leaving cmap_vmin
as None and causing 'int - NoneType' TypeError in the continuous
cmap scoring path.
```

## Reporting Bugs

Open an issue on GitHub with:

1. **Steps to reproduce** — ideally a minimal command line invocation.
2. **Expected vs. actual behaviour.**
3. **Python version** (`python3 --version`) and **OS**.
4. **Input data** — if possible, attach or link a minimal TSV / FASTA file
   that triggers the bug. For GISAID data, describe the subset without
   sharing protected sequences.

## Project Structure

```
src/mutation_scatter_plot/
├── mutation_scatter_plot/      # Main scatter plot & timeline modules
│   ├── __init__.py             # mutation_scatter_plot entry point
│   ├── cli.py                  # scatter plot CLI
│   ├── core.py                 # shared BLOSUM scoring, colormap logic
│   ├── colorbar_helpers.py     # colorbar rendering (matplotlib + Bokeh)
│   ├── timeline.py             # timeline scatter plot renderer
│   └── timeline_cli.py         # timeline CLI
├── calculate_codon_frequencies/  # codon frequency calculator
└── scripts/                    # standalone helper scripts
tests/
├── inputs/                     # test input FASTA + TSV files
├── outputs/                    # golden baseline outputs
├── test_mutation_scatter_plot.py  # main regression suite
├── test_dynamic_colormap.py    # colormap unit tests
└── test_scaling_logic.py       # scaling mode tests
scripts/                        # pipeline shell scripts
docs/                           # technical documentation
```

## Documentation

- **README.md** — user-facing documentation, installation, usage examples.
- **docs/** — technical design documents (scaling, profiling, testing).
- **Docstrings** — all public functions should have PEP 257 docstrings.

If your change modifies user-visible behaviour (new CLI flag, changed
output format, etc.), please update the README accordingly.

## Citation

If you use this software in your research, please cite:

> Shoshany A., Tian R., Padilla-Blanco M., Hruška A., Baxová K., Zoler E.,
> Mokrejš M., Schreiber G., Zahradník J. (submitted) In Vitro and Viral
> Evolution Convergence Reveal the Selective Pressures Driving Omicron
> Emergence.
> [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.23.650148v1)

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).
Please read it before participating.

## Questions?

If you are unsure whether a change would be welcome, open an issue first
to discuss the idea before investing time in a pull request.
