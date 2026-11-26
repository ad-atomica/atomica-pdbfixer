# Typed PDBFixer and upstream output comparison

On November 26, 2026, all nine controlled PDB and mmCIF outputs matched upstream byte for byte. With default CLI settings, hydrogen coordinates varied across repeated runs of each implementation.

## Revisions

| Source | Revision |
| --- | --- |
| Upstream `openmm/pdbfixer` | [`ef678254109cda7fb7d9d949a2d4ac2dce9ecb7a`](https://github.com/openmm/pdbfixer/commit/ef678254109cda7fb7d9d949a2d4ac2dce9ecb7a) |
| Tested `atomica-pdbfixer` | `709126992f61199b8975a987cd643499d37429d3` |

The tested revision and the commit that records its SHA share the same `pdbfixer/` Git tree, `1d2c8b3763c84e53a46c755b330b032c65d48b17`. To verify the package tree for this report, run the following from the repository root.

```sh
revision=$(git log -1 --diff-filter=A --format=%H -- docs/2026-11-26-typed-pdbfixer.md)
git rev-parse "$revision:pdbfixer"
```

## Changes in this commit

- Added mise tasks for environment setup, tests, Ruff formatting and linting, and ty checks. Applied Ruff across the Python sources and tests.
- Added type declarations for the OpenMM units, arithmetic operators, and compiled residue matcher used by PDBFixer. Updated the CCD reader's return type to include `None` when a definition has no atom data.
- Changed constructor validation to raise `ValueError` for invalid source arguments or an empty topology. The no-argument call `PDBFixer()` raised `Exception` upstream and `ValueError` in Atomica, with the same message.
- Closed local input files and CCD download streams on error. The web UI now stops after a failed download or export, and its HTTP responses use the encoded byte length.
- Updated the CLI tests to pass arguments directly to the executable. The simulation tests now close their temporary output files on failure and log the failing PDB identifier and traceback.

## Controlled outputs

Both revisions ran in separate Python processes on Linux x86-64 with Python 3.12.16 and OpenMM 8.6.1. Each process verified which checkout it imported. The comparison used the complete file bytes, including headers.

| Input and operation | Output | SHA-256 shared by both revisions |
| --- | --- | --- |
| `1BHL.pdb`, CLI with `--add-atoms=none` | PDB | `8106007340a9ee2395c4a15b2463bd8fbff5883f8354645897a2ec075aff35c9` |
| `4JSV.pdb`, same CLI options | PDB | `cd921f95716ab33b6321a20d8a22abc6c3361f37abc4741e0d44c93b848fa9ec` |
| `CAS.pdb`, same CLI options | PDB | `9d84dc052c97b3492ae0bca6df8137788c10ec8c08325df1b83f00e9f569d511` |
| `alanine-dipeptide.pdb`, same CLI options | PDB | `0657467bb83f2dc36fda3b8aada7a2bd980a3c8607ad5eebbad3da2fa2b99f5e` |
| `test.pdb`, same CLI options | PDB | `73859adad5b3ee8466f997b96c9739e2836b76b8b25046be99b3d13b5c991dfa` |
| Alanine dipeptide with its `CB` atom removed, seeded API repair | PDB | `310972809bc293d618059ee7593535c29f439e3bcf3c265bc887cd5276995eee` |
| Same repair | mmCIF | `803bc209b1dbea825162731929c4ad4959a54d6c58a25fadcaf80a8ef2d01004` |
| N-terminal gap fixture, seeded API repair | PDB | `05b6ea0e78feeede71edd252a11ed49214b01fb5525f70e7ef70d156094b98ad` |
| Same repair | mmCIF | `6e1906623a31a2430fabf5e8b1bcfafdc7d37662e1c914d17f7432c99b7e0a94` |

For the five CLI comparisons, both processes received the same absolute input path under `/home/vercel-sandbox/atomica-pdbfixer/pdbfixer/tests/data/`. The CLI writes that path into the PDB header. OpenMM also records its version and the date in its PDB and mmCIF output, so those values affect the hashes.

For the API comparisons, both processes used the OpenMM `Reference` platform, Python's `random.seed(48271)`, and `addMissingAtoms(seed=48271)`. They called `findMissingResidues()` and `findMissingAtoms()` before adding atoms, then called `addMissingHydrogens()` and wrote both formats with `keepIds=True`. The environment used `PYTHONHASHSEED=0`, `OPENMM_CPU_THREADS=1`, and `OPENBLAS_NUM_THREADS=1`. Both versions found one missing heavy atom in the alanine case and two missing residues plus nine missing heavy atoms in the terminal-gap fixture. [Terminal-gap fixture](../pdbfixer/tests/test_add_missing_residues.py)

## Default CLI output

Each revision ran twice on the original `alanine-dipeptide.pdb` with default repair options and no random seed. All four files contained the same 22 atom records in the same order.

| Pair compared | Records with coordinate differences |
| --- | --- |
| Upstream run 1 and Atomica run 1 | 12 hydrogens |
| Upstream runs 1 and 2 | 12 hydrogens |
| Atomica runs 1 and 2 | 12 hydrogens |

In all three comparisons, only the coordinate fields of those hydrogens changed. OpenMM uses random initial positions when adding hydrogens, so even repeated upstream runs can produce different file hashes for this input.
