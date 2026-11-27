import io
from pathlib import Path

import numpy as np
from openmm import Platform, unit

import pdbfixer
from pdbfixer import _chirality

# SEQRES starts two residues before the first resolved residue (ALA 1), so the
# residues built at the N terminus should be numbered -1 and 0.
PDB = """\
SEQRES   1 A    5  MET GLY ALA VAL GLU
ATOM      1  N   ALA A   1      -1.200   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       1.200   0.400   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.400   1.600   0.000  1.00 20.00           O
ATOM      5  N   VAL A   2       2.600   0.000   0.000  1.00 20.00           N
ATOM      6  CA  VAL A   2       3.800   0.000   0.000  1.00 20.00           C
ATOM      7  C   VAL A   2       5.000   0.400   0.000  1.00 20.00           C
ATOM      8  O   VAL A   2       5.200   1.600   0.000  1.00 20.00           O
ATOM      9  N   GLU A   3       6.400   0.000   0.000  1.00 20.00           N
ATOM     10  CA  GLU A   3       7.600   0.000   0.000  1.00 20.00           C
ATOM     11  C   GLU A   3       8.800   0.400   0.000  1.00 20.00           C
ATOM     12  O   GLU A   3       9.000   1.600   0.000  1.00 20.00           O
TER      13      GLU A   3
END
"""


def test_negative_residue_ids():
    fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(PDB))
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 0): ["MET", "GLY"]}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    ids = [residue.id for residue in fixer.topology.residues()]
    assert ids == ["-1", "0", "1", "2", "3"]


def test_copies_minimized_positions_to_the_corresponding_rebuilt_atoms(monkeypatch):
    fixture = Path(__file__).parent / "data" / "mixed-chain-gaps.pdb"
    fixer = pdbfixer.PDBFixer(
        filename=str(fixture), platform=Platform.getPlatformByName("Reference")
    )
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    spectator = next(res for res in fixer.topology.residues() if res.chain.id == "X")
    # Request a terminal oxygen on the untemplated chain that relaxation omits.
    fixer.missingTerminals[spectator] = ["OXT"]

    def atomKey(atom):
        return (atom.residue.chain.id, atom.residue.id, atom.name)

    original = {
        atomKey(atom): fixer.positions[atom.index].value_in_unit(unit.nanometer)
        for atom in fixer.topology.atoms()
    }
    minimized = {}
    warnAboutChirality = _chirality.warnAboutChirality

    def captureMinimizedPositions(centers, positions, *, stacklevel):
        # Observe the final relaxation state at the existing warning call.
        for center in centers:
            for atom in center.residue.atoms():
                minimized[atomKey(atom)] = positions[atom.index].value_in_unit(
                    unit.nanometer
                )
        warnAboutChirality(centers, positions, stacklevel=stacklevel + 1)

    monkeypatch.setattr(_chirality, "warnAboutChirality", captureMinimizedPositions)
    fixer.addMissingAtoms(seed=1)

    returned = {
        atomKey(atom): fixer.positions[atom.index].value_in_unit(unit.nanometer)
        for atom in fixer.topology.atoms()
    }
    for key, position in minimized.items():
        np.testing.assert_array_equal(returned[key], position, err_msg=str(key))
    outputOnly = ("X", "90", "OXT")
    # The retained CA–C bond is along x, so reflecting O changes its y offset.
    np.testing.assert_allclose(returned[outputOnly], (3.15, 2.88, 3.0))
    rebuilt = returned.keys() - original.keys() - {outputOnly}
    assert {key[0] for key in rebuilt} == {"A", "B"}
    assert rebuilt == minimized.keys()
    for key, position in original.items():
        np.testing.assert_array_equal(returned[key], position, err_msg=str(key))
