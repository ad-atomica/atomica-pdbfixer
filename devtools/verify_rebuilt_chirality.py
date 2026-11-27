"""Fetch, measure, summarize, and compare rebuilt Cα handedness for the rebuilt chirality report.

devtools/verify_rebuilt_chirality.sh runs the whole comparison.  The measure command repairs the
structures with whichever pdbfixer is on the import path.
"""

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path
from typing import NamedTuple
from urllib.request import urlopen

import numpy as np
from openmm import Platform, app, unit

if not __debug__:
    sys.exit("The checks in this script are assertions, so run it without -O.")

# This runner also measures v1.12 and bdf98bc, which lack pdbfixer._chirality, so it defines
# its own quartet, template threshold, and volume instead of importing them.

QUARTET = ("N", "CA", "C", "CB")
CB = QUARTET.index("CB")
MIN_TEMPLATE_VOLUME = 0.1


class Structure(NamedTuple):
    """A deposited structure the report measured.

    sha256 identifies the measured file.  maskStarts holds the first residue number of each
    three-residue window to delete and rebuild.  None repairs the deposited gaps without
    deleting anything.
    """

    sha256: str
    maskStarts: tuple[int | None, ...]

    def requireMeasured(self, data, source):
        """Exit unless data, read from source, is the file the report measured."""
        if hashlib.sha256(data).hexdigest() != self.sha256:
            sys.exit(f"{source} differs from the file the report measured.")


STRUCTURES = {
    "6FNL": Structure(
        "5531cbbca15b288be2405c7674f1f09c2eaa8969dfcc3706b63c4ecca42123d5", (None,)
    ),
    "3DRD": Structure(
        "1ce2a95945536a2580d31b92a3d90e3ba2682e1e91f8be43f5a2e8d11738bf43", (None,)
    ),
    "1CRN": Structure(
        "42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811", (12, 24, 35)
    ),
    "1UBQ": Structure(
        "d4a6812d8951cf6594e6a0763f089e35f5a80b62acb3c117b2c5565228a7b161", (20, 39, 58)
    ),
}

# The labels for the measured revisions.  compare also checks placement on the revisions with
# the placement fix, and final handedness on those with the restraint.

REVISIONS = ("release", "parent", "fix")
WITH_PLACEMENT_FIX = ("parent", "fix")
WITH_RESTRAINT = ("fix",)


class Residue(NamedTuple):
    chain: str
    resid: str
    icode: str
    name: str


# The columns that identify one measured site, then the measurements.  summarize counts sites
# per GROUP.

SITE = ("revision", "platform", "pdb", "mask_start", "seed", *Residue._fields)
MEASUREMENTS = ("template", "deposited", "placed", "final", "cb_error_angstrom")
COLUMNS = SITE + MEASUREMENTS
GROUP = ("revision", "platform", "seed", "pdb")


def centers(topology, positions):
    """Map each residue with a complete quartet to the quartet's coordinates in nm."""
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    result = {}
    for residue in topology.residues():
        atoms = {atom.name: atom.index for atom in residue.atoms()}
        if all(name in atoms for name in QUARTET):
            key = Residue(
                residue.chain.id,
                residue.id,
                residue.insertionCode.strip(),
                residue.name,
            )
            assert key not in result, key
            result[key] = xyz[[atoms[name] for name in QUARTET]]
    return result


def volume(points):
    """Signed volume from unit CA-to-N, CA-to-C, and CA-to-CB vectors."""
    vectors = points[[0, 2, 3]] - points[1]
    value = np.linalg.det(vectors / np.linalg.norm(vectors, axis=1)[:, None])
    assert np.isfinite(value)
    return value


def measure(path, start, seed, platform):
    """Repair one structure and yield the measured columns for each rebuilt center.

    If start is not None, the three residues starting there are deleted first, so the rebuilt
    centers can be compared with their deposited coordinates.
    """
    # Imported here so the other commands run without a pdbfixer revision on the import path.

    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(path), platform=platform)
    deposited = centers(fixer.topology, fixer.positions)
    if start is not None:
        chain = next(fixer.topology.chains())
        residues = [r for r in chain.residues() if start <= int(r.id) < start + 3]
        assert len(residues) == 3
        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.delete(residues)
        fixer.topology, fixer.positions = modeller.topology, modeller.positions

    retained = centers(fixer.topology, fixer.positions)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    topology, positions, _, _ = fixer._addAtomsToTopology(
        heavyAtomsOnly=True, omitUnknownMolecules=True
    )
    placed = centers(topology, positions)
    fixer.addMissingAtoms(seed=seed)
    final = centers(fixer.topology, fixer.positions)

    assert final.keys() == placed.keys()
    for key, points in retained.items():
        np.testing.assert_array_equal(final[key], points)
    rebuilt = final.keys() - retained.keys()
    assert rebuilt
    if start is not None:
        assert rebuilt == deposited.keys() - retained.keys()
    for key in sorted(rebuilt):
        template = fixer.templates[key.name]
        (ideal_points,) = centers(template.topology, template.positions).values()
        ideal = volume(ideal_points)
        assert abs(ideal) >= MIN_TEMPLATE_VOLUME
        truth = deposited.get(key)
        yield {
            "pdb": path.stem,
            "mask_start": start,
            "seed": seed,
            **key._asdict(),
            "template": ideal,
            "deposited": "" if truth is None else volume(truth),
            "placed": volume(placed[key]),
            "final": volume(final[key]),
            "cb_error_angstrom": ""
            if truth is None
            else 10 * np.linalg.norm(final[key][CB] - truth[CB]),
        }


def reference(row):
    """Return the volume whose sign the site should have.

    That is the deposited volume if the site was observed, else the template's.
    """
    return float(row["deposited"] or row["template"])


def readSites(path):
    """Map each site in a CSV written by the measure command to its row."""
    rows = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = tuple(row[column] for column in SITE)
            assert key not in rows, ("duplicate site", key)
            rows[key] = row
    return rows


def runFetch(args):
    """Download any missing deposited structures and check every file's SHA-256."""
    args.pdb_dir.mkdir(parents=True, exist_ok=True)
    for pdb, structure in STRUCTURES.items():
        path = args.pdb_dir / f"{pdb}.pdb"
        if path.exists():
            structure.requireMeasured(path.read_bytes(), path)
            continue

        # Write only a verified download, and rename it into place so a partial write is not
        # mistaken for a complete file.

        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        with urlopen(url) as response:
            data = response.read()
        structure.requireMeasured(data, url)
        partial = path.with_name(f"{path.name}.partial")
        partial.write_bytes(data)
        partial.replace(path)


def runMeasure(args):
    """Write one revision's CSV to stdout."""
    if args.seed == 0:
        sys.exit("seed must be nonzero because OpenMM chooses a fresh seed for 0")
    platform = Platform.getPlatformByName(args.platform)
    writer = csv.DictWriter(sys.stdout, COLUMNS, lineterminator="\n")
    writer.writeheader()
    for pdb in args.pdb:
        for start in STRUCTURES[pdb].maskStarts:
            for row in measure(args.pdb_dir / f"{pdb}.pdb", start, args.seed, platform):
                writer.writerow(
                    {"revision": args.revision, "platform": platform.getName(), **row}
                )


def runSummarize(args):
    """Print rebuilt centers and inversions at placement and in the final output."""
    sites = Counter()
    placed = Counter()
    final = Counter()
    for row in readSites(args.csv).values():
        key = tuple(row[column] for column in GROUP)
        expected = reference(row)
        sites[key] += 1
        placed[key] += expected * float(row["placed"]) < 0
        final[key] += expected * float(row["final"]) < 0
    print(*GROUP, "centers", "inverted_at_placement", "inverted_final")
    for key in sorted(sites):
        print(*key, sites[key], placed[key], final[key])


def runCompare(args):
    """Check that a new measurement covers the recorded sites with the recorded handedness."""
    recorded = readSites(args.recorded)
    measured = readSites(args.measured)
    assert recorded.keys() == measured.keys(), recorded.keys() ^ measured.keys()
    for key, row in measured.items():
        expected = reference(recorded[key])
        assert expected * reference(row) > 0, ("reference", key)
        if row["revision"] in WITH_PLACEMENT_FIX:
            assert expected * float(row["placed"]) > 0, ("placement", key)
        if row["revision"] in WITH_RESTRAINT:
            assert expected * float(row["final"]) > 0, ("relaxation", key)
    print(f"{len(measured)} sites match the recorded handedness.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True)

    command = commands.add_parser("fetch", help=runFetch.__doc__)
    command.add_argument("pdb_dir", type=Path)
    command.set_defaults(run=runFetch)

    command = commands.add_parser("measure", help=runMeasure.__doc__)
    command.add_argument(
        "pdb_dir", type=Path, help="directory of the four deposited PDB files"
    )
    command.add_argument("seed", type=int, help="nonzero OpenMM random seed")
    command.add_argument(
        "--revision",
        required=True,
        choices=REVISIONS,
        help="label for the pdbfixer revision being measured",
    )
    command.add_argument("--platform", default="Reference", help="OpenMM platform name")
    command.add_argument(
        "--pdb",
        nargs="+",
        choices=STRUCTURES,
        default=list(STRUCTURES),
        help="PDB IDs to run (default: all four)",
    )
    command.set_defaults(run=runMeasure)

    command = commands.add_parser("summarize", help=runSummarize.__doc__)
    command.add_argument("csv", type=Path)
    command.set_defaults(run=runSummarize)

    command = commands.add_parser("compare", help=runCompare.__doc__)
    command.add_argument("recorded", type=Path)
    command.add_argument("measured", type=Path)
    command.set_defaults(run=runCompare)

    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
