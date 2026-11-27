# Copy minimized coordinates to the corresponding rebuilt atom

## Intent

`addMissingAtoms()` builds one topology for relaxation and another for the
returned structure. A chain with only untemplated observed residues can be
omitted from relaxation even when a missing residue in that chain has a
template. Pairing the added atoms by list position can then assign minimized
coordinates to another chain. Each rebuilt atom that participates in relaxation
should receive its own minimized position. Retained atoms keep their input
positions, and atoms added only to the output keep their initial positions.

## References

The fork baseline is `2b0fe80`, which includes the chirality restraint and
copies coordinates by list position. The issue also reproduces the mismatch
on PDBFixer 1.12. [Fork issue #7](https://github.com/ad-atomica/atomica-pdbfixer/issues/7) · [Upstream release](https://github.com/openmm/pdbfixer/releases/tag/v1.12)

## Commit rationale

| Commit (subject or SHA) | What changed | Why it was needed |
| --- | --- | --- |
| `fix: copy minimized coordinates to the corresponding rebuilt atom` | Include chains with templated insertions in relaxation. Record each added atom's site in both topology builds and copy coordinates by site. Add a coordinate-level regression. | Filtering can give the two builds different chains and atom orders, so a positional copy can silently corrupt the returned structure. |

## Reapplying after an upstream update

1. Run the regression below against the updated `addMissingAtoms()`. Both
   repaired chains should receive their own minimized coordinates and the
   output-only terminal oxygen should retain its initial placement.
2. When `_addAtomsToTopology()` filters chains for relaxation, consider both
   observed residues and templated residues scheduled for insertion. Leave
   wholly untemplated chains outside relaxation when they have no templated
   insertion.
3. In both builds, record each added atom by its input chain index, residue
   index within the rebuilt chain, atom name, and occurrence of that name
   among the additions to that residue. The builds share residue order within
   an included chain even when filtering changes topology-wide indices.
   Copy each relaxed atom's position to its corresponding output atom.

## Verification

The regression extends the issue fixture with an untemplated chain between
the two repaired chains and requests a terminal oxygen for that chain. It
compares every relaxed, newly added atom's final position with its returned
position and also checks all retained positions and the output-only oxygen.

```sh
mise run setup
OPENMM_CPU_THREADS=2 .venv/bin/python -m pytest -q \
    pdbfixer/tests/test_add_missing_residues.py -k corresponding
mise run check
OPENMM_CPU_THREADS=2 mise run test
```

On the baseline, the returned nitrogen of B:ALA2 differed from its relaxed
position. Including both repaired chains alone still failed at that nitrogen
with the output-only oxygen present. Copying by site alone left A:ALA2
outside relaxation. With both changes, the regression passes.

On the original issue fixture, all 18 rebuilt atoms received their exact
minimized positions. The same check passed with duplicate chain IDs.

All 50 tests passed, and `mise run check` passed formatting, lint, and type
checks.
