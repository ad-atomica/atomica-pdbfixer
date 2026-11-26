# Preserve template handedness during missing-atom placement

## Intent

When `addMissingAtoms()` fits a residue template to the observed atoms, the
closest fit can be a reflection, especially when those atoms lie nearly in one
plane. A reflected template can place a rebuilt alanine Cβ on the wrong side of
the backbone. The intended fit is the closest rotation with determinant +1.
It must preserve the template's handedness and recover the original rotation
and translation when the point sets uniquely determine that rigid motion.

The baseline multiplies the two orthogonal factors from the singular value
decomposition. Their product can have determinant −1, which reverses the signed
volume around the Cα atom. When that happens, the fix flips the direction with
the smallest singular value. This gives determinant +1 with the smallest
possible increase in squared fitting error.

## References

The fork baseline was `2c5b006`. The correction and regression tests come from
upstream commit `9b6c8f2`. [Fork issue #1](https://github.com/ad-atomica/atomica-pdbfixer/issues/1) · [Upstream PR #355](https://github.com/openmm/pdbfixer/pull/355) · [Upstream patch](https://github.com/openmm/pdbfixer/commit/9b6c8f2775dbc5f95d0e7adfb387c88fde60da24)

## Commit rationale

| Commit (subject or SHA) | What changed | Why it was needed |
| --- | --- | --- |
| `docs: define fork contribution rules and a report template` | Added `AGENTS.md` and a template for reports covering each commit. | Keeping each minimal fix in a single commit and its own PR lets maintainers review, reapply, or revert it independently. The report records the intent and evidence for future upstream updates. |
| `fix: preserve template handedness when placing missing atoms` | Ported the upstream Kabsch correction and both regression tests, and recorded their results in this report. | The unconstrained fit could mirror a template. The tests verify the rotation and show its effect on rebuilt atom positions. |

## Reapplying after an upstream update

1. Run the regression cases below on the updated upstream checkout, porting the tests if needed. If all three pass and the fit guarantees a proper rotation, use the upstream implementation.
2. If the fit permits a reflection, follow `addMissingAtoms()` to the operation that aligns centered template coordinates to observed coordinates. In this revision it is `_overlayPoints()` in `pdbfixer/pdbfixer.py`.
3. For an SVD-based fit, reverse the direction with the smallest singular value when the unconstrained fit has negative determinant. This minimizes the fitting error among proper rotations. Preserve the updated caller's coordinate order, matrix convention, and centering translations.

## Verification

The same two upstream regression tests ran on both revisions with Python
3.12.16, NumPy 2.5.5, OpenMM 8.6.1, and `OPENMM_CPU_THREADS=2`. The rotation
test has two cases. [Regression tests](../pdbfixer/tests/test_find_positions.py)

| Measurement | Baseline (`2c5b006`) | Fix (commit adding this report) |
| --- | --- | --- |
| Rotated point set | det(R) = +1, recovered known rotation and translations | det(R) = +1, recovered known rotation and translations |
| Mirrored point set | det(R) = −1, returned a reflection | det(R) = +1, recovered the closest proper rotation |
| Rebuilt alanine Cβ atoms with opposite handedness in 1BHL | 10 of 14 | 0 of 14 |

The rebuild test deletes all 14 alanine Cβ atoms and calls `addMissingAtoms()`.
For each alanine it compares the sign of the scalar triple product of the
Cα→N, Cα→C, and Cα→Cβ vectors before deletion and after rebuilding.
A sign change places the rebuilt Cβ on the opposite side of the N–Cα–C plane.
The determinant measures the alignment step. The alanine count measures the
final coordinates after minimization on this fixture.

To reproduce the table, run the following from the PR checkout after
`mise run setup`. It extracts each revision into a temporary directory and
uses the same tests and input file for both. The two print statements expose
the measurements just before the original assertions.

```sh
python_bin="$PWD/.venv/bin/python"
scratch=$(mktemp -d)
test_file=pdbfixer/tests/test_find_positions.py
fixed_revision=$(git log -1 --diff-filter=A --format=%H -- docs/2026-11-26-preserve-template-handedness.md)
for revision in 2c5b006 "$fixed_revision"; do
    run_dir="$scratch/$revision"
    mkdir "$run_dir"
    git archive "$revision" pdbfixer | tar -x -C "$run_dir"
    git show "$fixed_revision:$test_file" | awk '
        /assert np.allclose\(rotate, known\)/ {
            print "    print(\"det(R) =\", round(np.linalg.det(rotate), 6))"
        }
        /assert _alanineHandedness\(fixer\) == originalHandedness/ {
            print "    rebuilt = np.array(_alanineHandedness(fixer))"
            print "    print(\"opposite handedness:\", np.count_nonzero(rebuilt != originalHandedness), \"of\", len(rebuilt))"
        }
        { print }
    ' > "$run_dir/$test_file"
    (cd "$run_dir" && OPENMM_CPU_THREADS=2 "$python_bin" -m pytest -vs --tb=no \
        "$test_file" -k 'overlayPointsRecoversRotation or addedCBChirality')
done
```

Expected results are one pass and two failures on the baseline, then three
passes with the fix. The full fork suite passed all 30 tests with
`OPENMM_CPU_THREADS=2 mise run test`. Formatting, lint, and type checks passed
with `mise run check`.
