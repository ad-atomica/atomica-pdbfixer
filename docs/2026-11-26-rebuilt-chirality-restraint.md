# Preserve Cα handedness during missing-atom relaxation

## Intent

`addMissingAtoms()` can place a residue with the correct handedness and still
invert its Cα during energy minimization or the dynamics used to resolve
clashes. The soft force field allows both orientations. A rebuilt center should
finish with the handedness of its own residue template, including a mirrored
template supplied through `registerTemplate()`. A center made entirely of
input atoms should keep the coordinates supplied by the caller. Warn and
identify any restrained center that finishes relaxation with the opposite
sign or an undefined sign (zero or `NaN` signed volume).

## References

The fork baseline is `bdf98bc`, which includes the correction to the initial
template fit. The fix ports the restraint from Sahil Athalye's second
patch, `02-chirality-restrained-minimization.patch`, attached to upstream
issue #352. [Fork issue #2](https://github.com/ad-atomica/atomica-pdbfixer/issues/2) · [Upstream issue and patch](https://github.com/openmm/pdbfixer/issues/352)

## Commit rationale

| Commit (subject or SHA) | What changed | Why it was needed |
| --- | --- | --- |
| `fix: preserve rebuilt Cα chirality during minimization` | Added `pdbfixer/_chirality.py`, which selects rebuilt centers, restrains them to their templates, and warns about any that finish without their template's sign. `addMissingAtoms()` calls it twice. Added a regression for both template signs and unit tests for the selection rules, restraint, and warning. | Bond lengths and angles can survive a reflection. The template supplies the intended sign. |

## Reapplying after an upstream update

1. Run the regression below on the updated `addMissingAtoms()`. Check both template signs and the fully retained alanine.
2. Keep `pdbfixer/_chirality.py`, which at runtime depends only on NumPy and OpenMM. `pdbfixer.py` imports it with `from . import _chirality`, so run PDBFixer through the `pdbfixer` command, `python -m pdbfixer.pdbfixer`, or the package, not as the file `pdbfixer/pdbfixer.py`.
3. In `addMissingAtoms()`, find the OpenMM system shared by the initial minimization, clash-recovery dynamics, and final minimization. Call `_chirality.addChiralityRestraints(system, newAtoms, self._getTemplate)` before creating its context. Passing the placement lookup makes relaxation honor the same template overrides. The module's docstrings and comments define which centers are restrained, the restraint energy, and the reasons for `restraintStrength` and `minTemplateVolume`.
4. Call `_chirality.warnAboutChirality(chiralCenters, state.getPositions(), stacklevel=2)` on the selected final state, so the warning names each residue and points at the caller of `addMissingAtoms()`.

## Verification

### Controlled peptide regression

The regression inserts one alanine between the acetyl cap and the existing
alanine in the bundled capped peptide. Both revisions start with the correct
template placement. The table shows normalized signed
volumes on Linux with OpenMM 8.6.1, the Reference platform, and seed 1.
[Regression and fixture setup](../pdbfixer/tests/test_find_positions.py)

| Template | At placement | Final, `parent` | Final, `fix` |
| --- | ---: | ---: | ---: |
| Standard | +0.714 | −0.929 | +0.719 |
| Mirrored | −0.714 | +0.605 | −0.714 |

Both centers inverted on `parent` and kept their template sign on `fix`. The
fully retained alanine stayed at +0.754 in all four runs, including the two
with a negative template target. The sign checks also passed on the Reference
and CPU platforms with seeds 1, 17, and 12345.

Unit tests check that marking any single quartet atom as new selects the
center, that several new atoms in one residue select it once, and that none is
selected when only O is new or there is no template.
They check that an incomplete or planar template, or one with a zero-length
bond, gives no target, and that the restraint force holds one bond per center.
They also check that the OpenMM energy matches the NumPy volume and that a
mirrored final structure produces a warning naming the residue.
[Restraint unit tests](../pdbfixer/tests/test_chirality.py)

A separate regression rebuilds Cβ on alanine with collinear, fixed N, Cα,
and C. Every Cβ position gives exactly zero signed volume. Previously the
repair returned zero silently. The test checks the returned coordinates, the
warning naming `B:ALA17A`, and that the warning points at the caller of
`addMissingAtoms()`.

```sh
.venv/bin/python -m pytest -q pdbfixer/tests/test_chirality.py -k zeroRebuilt
```

To reproduce the sign checks, run the following from this checkout after
`mise run setup`. It checks the standard template, a registered mirror, and a
mirror supplied only by a subclass's template lookup. Each case must fail on
`parent` at the final-handedness assertion and pass on this checkout. The
check runs in its own shell, so an unexpected pass or pytest invocation error
stops it without closing yours.

```sh
bash -eu <<'SH'
python_bin="$PWD/.venv/bin/python"
test_file=pdbfixer/tests/test_find_positions.py
parent=$(mktemp -d)
trap 'rm -rf "$parent"' EXIT
git archive bdf98bc pdbfixer | tar -x -C "$parent"
cp "$test_file" "$parent/$test_file"
# -P and PYTHONPATH make Python import the parent package, not the editable install.
for case in standard mirrored lookup; do
    status=0
    PYTHONPATH="$parent" "$python_bin" -P -m pytest -q \
        "$parent/$test_file::test_chiralityDuringMinimization[$case]" || status=$?
    if [ "$status" -ne 1 ]; then
        echo "Expected pytest exit 1 for parent $case, got $status" >&2
        exit 1
    fi
done
"$python_bin" -m pytest -q "$test_file" -k chiralityDuringMinimization
SH
```

`OPENMM_CPU_THREADS=2 mise run test` passed all 49 tests. `mise run check`
passed the formatting, lint, and type checks.
