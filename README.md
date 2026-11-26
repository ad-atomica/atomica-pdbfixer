[![PDBFixer Continuous Integration Workflow](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml/badge.svg)](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml)

PDBFixer
========

PDBFixer is an easy to use application for fixing problems in Protein Data Bank files in preparation for simulating them.  It can automatically fix the following problems:

- Add missing heavy atoms.
- Add missing hydrogen atoms.
- Build missing loops.
- Convert non-standard residues to their standard equivalents.
- Select a single position for atoms with multiple alternate positions listed.
- Delete unwanted chains from the model.
- Delete unwanted heterogens.
- Build a water box for explicit solvent simulations.

See our [manual](https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html)

## Installation

PDBFixer can be installed with conda or mamba.

```
conda install -c conda-forge pdbfixer
```

Alternatively you can install from source, as described in the manual.

## Local development

With [mise](https://mise.jdx.dev/) installed, run these commands from the repository root.

```sh
mise trust
mise install
mise run setup
mise run test
```

Mise installs `uv`, Ruff, and ty. The setup task uses `uv` to download Python 3.12,
create `.venv`, and install this checkout in editable mode alongside `pytest`
and `setuptools`. To run Python or the CLI directly, activate the environment
first.

```sh
source .venv/bin/activate
python --version
pdbfixer --help
```

Run Ruff and ty through mise. The type checker uses the dependencies in `.venv`.

```sh
mise run format
mise run lint
mise run typecheck
```

Use `mise run format --check` to check formatting without editing files.
Run `mise run fix` to apply Ruff fixes and formatting, then `mise run check` to
verify formatting, lint, and types together.
