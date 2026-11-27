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
| `docs: record the deposited-structure chirality comparison` | Recorded the four-structure comparison and added the scripts and per-residue data that reproduce it. | The comparison tests whether relaxation loses the template sign in real repairs, including residues with known experimental coordinates. |

## Reapplying after an upstream update

1. Run the regression below on the updated `addMissingAtoms()`. Check both template signs and the fully retained alanine.
2. Keep `pdbfixer/_chirality.py`, which at runtime depends only on NumPy and OpenMM. `pdbfixer.py` imports it with `from . import _chirality`, so run PDBFixer through the `pdbfixer` command, `python -m pdbfixer.pdbfixer`, or the package, not as the file `pdbfixer/pdbfixer.py`.
3. In `addMissingAtoms()`, find the OpenMM system shared by the initial minimization, clash-recovery dynamics, and final minimization. Call `_chirality.addChiralityRestraints(system, newAtoms, self._getTemplate)` before creating its context. Passing the placement lookup makes relaxation honor the same template overrides. The module's docstrings and comments define which centers are restrained, the restraint energy, and the reasons for `restraintStrength` and `minTemplateVolume`.
4. Call `_chirality.warnAboutChirality(chiralCenters, state.getPositions(), stacklevel=2)` on the selected final state, so the warning names each residue and points at the caller of `addMissingAtoms()`.

## Verification

### Deposited structures

We ran the same repairs on three revisions, labeled as in the CSV. `release`
is stock PDBFixer 1.12 (`94cfa4c`). `parent` is the fork after the placement
fix (`bdf98bc`, PR #5). `fix` adds the restraint from PR #6.
Comparing `parent` and `fix` isolates the restraint. The 72 repairs in the
table used Linux x86-64, Python 3.12.16, NumPy 2.5.5, and OpenMM 8.6.1 on the
Reference platform.
[Upstream release](https://github.com/openmm/pdbfixer/releases/tag/v1.12) · [Placement fix](https://github.com/ad-atomica/atomica-pdbfixer/pull/5)

**All 98 distinct sites started with the expected handedness on `parent` and
`fix`.** Their signed volumes at placement matched exactly. Across three
seeds, `release` returned 161 inverted centers in 294 measurements, `parent`
returned 135, and `fix` returned zero. On `parent`, 81 of the 98 sites
inverted at least once.

The final columns give inversion counts for seeds **1, 17, 12345**, in that
order. Each seed repeats the same sites. Placement is independent of seed.

| Input | Sites per seed | Inverted at placement, `release` | Final, `release` | Final, `parent` | Final, `fix` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6FNL, deposited gaps | 31 | 9 | 15, 14, 19 | 18, 17, 11 | 0, 0, 0 |
| 3DRD, deposited gaps | 50 | 19 | 26, 25, 36 | 24, 21, 19 | 0, 0, 0 |
| 1CRN, held-out coordinates | 8 | 0 | 5, 7, 1 | 5, 7, 1 | 0, 0, 0 |
| 1UBQ, held-out coordinates | 9 | 4 | 4, 4, 5 | 3, 5, 4 | 0, 0, 0 |

- **Gaps in unmodified deposits.** Human EphB4 kinase (6FNL) has gaps of
  11, 20, and 4 residues. The *B. subtilis* enzyme (3DRD) has gaps of 29, 1,
  29, and 1. Both are examples from upstream issue #352. The 31 and 50 rebuilt
  centers match the non-glycine residues the depositors listed as missing in
  `REMARK 465`. Their coordinates are unobserved, so we compared each rebuilt
  center to its residue template. All 245 and 776 measurable centers in the
  respective deposits had the same positive sign.
  [6FNL](https://www.rcsb.org/structure/6FNL) · [3DRD](https://www.rcsb.org/structure/3DRD) · [Upstream report](https://github.com/openmm/pdbfixer/issues/352)
- **Held-out experimental coordinates.** We removed three consecutive,
  fully observed residues per repair from crambin (1CRN) and ubiquitin (1UBQ).
  We chose the windows before measuring outcomes, starting at the quarter,
  half, and three-quarter indices of the first protein chain. The chain A
  starts were 12, 24, and 35 in crambin and 20, 39, and 58 in ubiquitin.
  Each repair began from a fresh copy of the deposit. We compared the output
  against the atoms we removed.
  [1CRN](https://www.rcsb.org/structure/1CRN) · [1UBQ](https://www.rcsb.org/structure/1UBQ)

The runner calls `findMissingResidues()`, `findMissingAtoms()`, and
`addMissingAtoms(seed=...)`. It also observes initial placement with
`_addAtomsToTopology()`. At both stages it computes the determinant
of the unit vectors from Cα toward N, C, and Cβ with NumPy. A sign opposite
to the reference counts as an inversion. Glycine has no Cβ and is excluded.
The runner checks that every fully retained quartet has identical coordinates
before and after repair. [Measurement code](../devtools/verify_rebuilt_chirality.py)

On `parent`, 110/243 measurements in the deposited gaps disagreed with the
template, and 25/51 held-out measurements disagreed with the experimental
sign. Of the 135 inversions, 70 had an absolute volume above 0.5. On `fix`,
the smallest absolute final volume across all 294 was 0.637. The CSV is the
combined runner output, one row per revision, platform, seed, and rebuilt
center. Rows labeled `CPU` record the additional platform check below.
[Per-residue measurements](data/rebuilt-chirality.csv)

Ubiquitin Tyr59 illustrates the sign change against an experimental
reference. Its deposited volume was +0.743. With seed 1, `parent` and `fix`
both placed it at +0.705. `parent` finished at −0.808 and `fix` finished at
+0.705. The rebuilt Cβ was 4.44 Å and 3.89 Å from the deposited position,
respectively. Repeating the three ubiquitin windows on OpenMM's CPU platform
with seed 1 also gave 4/9 inverted centers on `release`, 3/9 on `parent`, and
0/9 on `fix`.

Across the 51 held-out measurements on Reference, the median Cβ distance
from its deposited position was 3.89 Å on `parent` and 4.06 Å on `fix`. The
restraint recovered local handedness in these structures. Recovering the full
loop conformation is a separate accuracy measure.

### Reproduce the deposited-structure comparison

Run the script from this checkout after `mise run setup`. It needs network
access. In the directory it is given, it creates a virtual environment with the
versions listed above. The runner's `fetch` command then downloads the four
deposits and checks their SHA-256 hashes, and the script extracts `release` and
`parent` with `git archive`. It measures `fix` from this checkout.
[Comparison script](../devtools/verify_rebuilt_chirality.sh)

For a quick check of Tyr59, run it in `quick` mode. It prints the nine
ubiquitin sites for `parent`, then for `fix`. Compare the `final` column for
chain A, residue 59.

```sh
bench=$(mktemp -d)
devtools/verify_rebuilt_chirality.sh "$bench" quick
rm -rf "$bench"
```

Without `quick`, the script measures each revision with seeds 1, 17, and 12345
on Reference and with seed 1 on CPU for ubiquitin. It runs these 12
single-threaded measurements in parallel, which took 15 minutes (81 CPU
minutes) on an 8-core Linux machine. It writes the combined CSV to
`$bench/rebuilt-chirality.csv` in the committed row order. With the versions
listed above, that file matched the committed CSV byte for byte. It then runs
the runner's `summarize` command, which prints the rebuilt centers, inversions
at placement, and inversions in the final output per revision, platform, seed,
and PDB ID. Last, it runs `compare`. Different numerical libraries can change
the volumes and the outcomes of unrestrained minimization, so `compare` checks
only that every recorded site is present exactly once, that placement on
`parent` and `fix` matches the reference handedness, and that `fix` retains it
after relaxation. Rerunning with the same directory reuses the environment,
inputs, and extracted revisions.

```sh
bench=$(mktemp -d)
devtools/verify_rebuilt_chirality.sh "$bench"
rm -rf "$bench"
```

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
