from io import StringIO
from pathlib import Path

import numpy as np
import openmm as mm
import pytest
from openmm import Platform, Vec3, app, unit

import pdbfixer
from pdbfixer._chirality import (
    addChiralityRestraints,
    findChiralCenters,
    normalizedChiralVolume,
    restraintStrength,
    templateVolume,
    warnAboutChirality,
)
from pdbfixer.pdbfixer import Template

_quartet = ("N", "CA", "C", "CB")


def _alanine():
    """Return a capped alanine fixer and its alanine's atoms by name."""
    fixer = pdbfixer.PDBFixer(
        filename=(Path(__file__).parent / "data" / "alanine-dipeptide.pdb").as_posix()
    )
    (residue,) = [r for r in fixer.topology.residues() if r.name == "ALA"]
    return fixer, {atom.name: atom for atom in residue.atoms()}


def _editedAlanineTemplate(fixer, edit):
    """Return the fixer's alanine template after edit() changes its positions.

    The positions are keyed by atom name.  Atoms that edit() removes are deleted from the template.
    """
    template = fixer.templates["ALA"]
    modeller = app.Modeller(template.topology, template.positions)
    atoms = list(modeller.topology.atoms())
    positions = {
        atom.name: position
        for atom, position in zip(
            atoms, modeller.positions.value_in_unit(unit.nanometer)
        )
    }
    edit(positions)
    modeller.delete([atom for atom in atoms if atom.name not in positions])
    names = [atom.name for atom in modeller.topology.atoms()]
    return Template(
        modeller.topology, [positions[name] for name in names] * unit.nanometer
    )


def _deleteCB(positions):
    del positions["CB"]


def _flatten(positions):
    for name, (x, y, _) in positions.items():
        positions[name] = Vec3(x, y, 0)


def _collapseCB(positions):
    positions["CB"] = positions["CA"]


def test_templateVolume():
    """The bundled L-alanine template has usable handedness.

    The input residue's own volume is 0.754.
    """
    fixer, _ = _alanine()
    assert templateVolume(fixer.templates["ALA"]) == pytest.approx(0.7136, abs=1e-4)


@pytest.mark.parametrize(
    "edit",
    [
        pytest.param(_deleteCB, id="incomplete"),
        pytest.param(_flatten, id="planar"),
        pytest.param(_collapseCB, id="zeroLengthBond"),
    ],
)
def test_templateWithoutUsableHandedness(edit):
    """A template without a complete, nonplanar center has no target volume."""
    fixer, _ = _alanine()
    assert templateVolume(_editedAlanineTemplate(fixer, edit)) is None


@pytest.mark.parametrize("newAtom", _quartet)
def test_findsCenterContainingAnyNewAtom(newAtom):
    """Marking any one quartet atom as new selects the center.

    Its target is the template's volume.
    """
    fixer, atoms = _alanine()
    (center,) = findChiralCenters([atoms[newAtom]], fixer._getTemplate)
    assert center.particles == tuple(atoms[name].index for name in _quartet)
    assert center.ideal == templateVolume(fixer.templates["ALA"])


def test_selectsEachResidueOnce():
    """A residue with several new quartet atoms gets one center."""
    fixer, atoms = _alanine()
    (center,) = findChiralCenters(
        [atoms[name] for name in _quartet], fixer._getTemplate
    )
    assert center.residue is atoms["CA"].residue


def test_ignoresCenterWithoutNewQuartetAtom():
    """A residue whose only new atom is outside the quartet gets no center."""
    fixer, atoms = _alanine()
    assert findChiralCenters([atoms["O"]], fixer._getTemplate) == []


def test_ignoresCenterWithoutTemplate():
    """A residue without a template gets no center."""
    _, atoms = _alanine()
    assert findChiralCenters([atoms["CB"]], {}.get) == []


def _system(fixer):
    """Return a System with one particle per atom of the fixer's topology."""
    system = mm.System()
    for _ in range(fixer.topology.getNumAtoms()):
        system.addParticle(1.0)
    return system


def test_addsForceOnlyForCenters():
    """No force is added without centers, and otherwise one force holds a bond per center."""
    fixer, atoms = _alanine()
    system = _system(fixer)
    assert addChiralityRestraints(system, [atoms["O"]], fixer._getTemplate) == []
    assert system.getNumForces() == 0
    (center,) = addChiralityRestraints(system, [atoms["CB"]], fixer._getTemplate)
    (force,) = system.getForces()
    assert force.getNumBonds() == 1
    particles, parameters = force.getBondParameters(0)
    assert tuple(particles) == center.particles
    assert tuple(parameters) == (center.ideal,)


def test_restraintEnergyMatchesVolume():
    """The OpenMM restraint evaluates the same normalized volume as normalizedChiralVolume()."""
    fixer, atoms = _alanine()
    system = _system(fixer)
    (center,) = addChiralityRestraints(system, [atoms["CB"]], fixer._getTemplate)

    # Generic coordinates make every coordinate difference in the restraint expression nonzero,
    # so a wrong term changes the energy.  The fixture itself lies largely in the z=0 plane.

    xyz = np.random.default_rng(0).normal(size=(fixer.topology.getNumAtoms(), 3))
    context = mm.Context(
        system, mm.VerletIntegrator(0.001), Platform.getPlatformByName("Reference")
    )
    context.setPositions([Vec3(*p) for p in xyz] * unit.nanometer)
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    volume = normalizedChiralVolume(xyz[list(center.particles)])
    assert energy.value_in_unit(unit.kilojoules_per_mole) == pytest.approx(
        restraintStrength * (volume / center.ideal - 1) ** 2
    )


def test_warnsAboutInvertedCenters():
    """A center whose final volume has the opposite sign to its template is reported by residue."""
    fixer, atoms = _alanine()
    centers = findChiralCenters([atoms["CB"]], fixer._getTemplate)
    mirrored = [
        Vec3(x, y, -z) for x, y, z in fixer.positions.value_in_unit(unit.nanometer)
    ] * unit.nanometer
    with pytest.warns(
        UserWarning, match="handedness does not match their templates: A:ALA2$"
    ):
        warnAboutChirality(centers, mirrored, stacklevel=1)


@pytest.mark.filterwarnings("error")
def test_noWarningForMatchingCenters():
    """A center that keeps its template's sign produces no warning."""
    fixer, atoms = _alanine()
    centers = findChiralCenters([atoms["CB"]], fixer._getTemplate)
    warnAboutChirality(centers, fixer.positions, stacklevel=1)


def test_warnsForZeroRebuiltChiralVolume():
    """Collinear retained N, CA, and C give zero volume for every position of the rebuilt CB.

    addMissingAtoms() reports that against its caller.
    """
    pdb = StringIO(
        """\
ATOM      1  N   ALA B  17A     -1.460   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA B  17A      0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA B  17A      1.520   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA B  17A      2.100   1.080   0.000  1.00  0.00           O
ATOM      5  OXT ALA B  17A      2.100  -1.080   0.000  1.00  0.00           O
END
"""
    )
    fixer = pdbfixer.PDBFixer(
        pdbfile=pdb, platform=Platform.getPlatformByName("Reference")
    )
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    with pytest.warns(
        UserWarning, match="handedness does not match their templates: B:ALA17A$"
    ) as record:
        fixer.addMissingAtoms(seed=1)
    assert [warning.filename for warning in record] == [__file__]
    atoms = {atom.name: atom.index for atom in fixer.topology.atoms()}
    xyz = np.asarray(fixer.positions.value_in_unit(unit.nanometer))
    assert normalizedChiralVolume(xyz[[atoms[name] for name in _quartet]]) == 0
