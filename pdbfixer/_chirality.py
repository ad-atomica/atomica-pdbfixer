"""Restrain rebuilt alpha-carbon stereocenters to the handedness of their residue templates.

addMissingAtoms() relaxes new atoms with the soft force field, which is invariant under inversion
of a tetrahedral carbon.  Without these restraints, minimization and clash-recovery dynamics can
turn a correctly placed residue into its mirror image.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.linalg as lin
import openmm as mm
from openmm import app, unit

if TYPE_CHECKING:
    from .pdbfixer import Template

# The atoms of the alpha-carbon stereocenter, in the order the functions below expect.  CA, the
# second atom, is the center.

chiralCenterAtoms = ("N", "CA", "C", "CB")

# A planar center costs one restraint strength in kJ/mol.  An exact mirror of the template costs
# four.  The planar penalty is large compared to kT at 300 K (2.5 kJ/mol) and the 43.9 kJ/mol
# improper torsion coefficient in soft.xml.  With a 5 fs step, the normalized restraint stayed
# stable at 10000 kJ/mol, and the unnormalized restraint became unstable above about 300 kJ/mol,
# as measured for the patch attached to https://github.com/openmm/pdbfixer/issues/352.

restraintStrength = 1000

# A template whose absolute normalized volume is below this has no usable handedness.  An
# undistorted tetrahedral center gives about 0.77.

minTemplateVolume = 0.1

# The same quantity as normalizedChiralVolume(), with a=p1-p2, b=p3-p2, and c=p4-p2.  It is
# written out in components because the equivalent product of sines of OpenMM's angle() and
# dihedral() has NaN forces when N, CA, and C are collinear.

restraintEnergy = (
    f"{restraintStrength}*(volume/ideal-1)^2;"
    "volume=(ax*(by*cz-bz*cy)+ay*(bz*cx-bx*cz)+az*(bx*cy-by*cx))"
    "/(distance(p1,p2)*distance(p3,p2)*distance(p4,p2));"
    "ax=x1-x2; ay=y1-y2; az=z1-z2;"
    "bx=x3-x2; by=y3-y2; bz=z3-z2;"
    "cx=x4-x2; cy=y4-y2; cz=z4-z2"
)


@dataclass(frozen=True)
class ChiralCenter:
    """A restrained stereocenter.  particles are in chiralCenterAtoms order."""

    residue: app.Residue
    particles: tuple[int, int, int, int]
    ideal: float

    @property
    def label(self) -> str:
        residue = self.residue
        return f"{residue.chain.id}:{residue.name}{residue.id}{residue.insertionCode.strip()}"

    def matchesTemplate(self, positions: unit.Quantity) -> bool:
        """Return whether the center's volume in positions has its template's sign.

        A zero volume, whose handedness is undefined, has neither sign.
        """
        volume = normalizedChiralVolume(_centerPoints(positions, self.particles))
        return volume * self.ideal > 0


def normalizedChiralVolume(points: np.ndarray) -> float:
    """Return the signed volume of the unit vectors from CA to N, C, and CB.

    points holds the N, CA, C, and CB positions in nm, as rows in that order.  The result is zero
    for a planar or degenerate center.
    """
    toN, toC, toCB = points[[0, 2, 3]] - points[1]
    scale = lin.norm(toN) * lin.norm(toC) * lin.norm(toCB)
    if scale == 0:
        return 0.0
    return float(np.dot(toN, np.cross(toC, toCB)) / scale)


def _centerPoints(positions: unit.Quantity, indices: Iterable[int]) -> np.ndarray:
    """Return the positions at indices in nm, as rows for normalizedChiralVolume()."""
    return np.array([positions[i].value_in_unit(unit.nanometer) for i in indices])


def _centerIndices(atoms: Iterable[app.Atom]) -> tuple[int, int, int, int] | None:
    """Return the indices of the N, CA, C, and CB atoms among atoms, or None if any is missing."""
    indices = {atom.name: atom.index for atom in atoms}
    if not all(name in indices for name in chiralCenterAtoms):
        return None
    n, ca, c, cb = (indices[name] for name in chiralCenterAtoms)
    return n, ca, c, cb


def templateVolume(template: Template | None) -> float | None:
    """Return a template's normalized chiral volume.

    Return None if there is no template or it has no usable handedness.
    """
    if template is None:
        return None
    indices = _centerIndices(template.topology.atoms())
    if indices is None:
        return None
    ideal = normalizedChiralVolume(_centerPoints(template.positions, indices))
    return ideal if abs(ideal) >= minTemplateVolume else None


def findChiralCenters(
    newAtoms: Iterable[app.Atom], getTemplate: Callable[[str], Template | None]
) -> list[ChiralCenter]:
    """Return each center that contains a new atom and can be restrained.

    A center can be restrained if its residue has all of chiralCenterAtoms and its template has
    usable handedness.  Each residue appears once, in the order its first new center atom appears.
    """
    centers = []
    for residue in dict.fromkeys(
        atom.residue for atom in newAtoms if atom.name in chiralCenterAtoms
    ):
        particles = _centerIndices(residue.atoms())
        ideal = templateVolume(getTemplate(residue.name))
        if particles is not None and ideal is not None:
            centers.append(ChiralCenter(residue, particles, ideal))
    return centers


def addChiralityRestraints(
    system: mm.System,
    newAtoms: Iterable[app.Atom],
    getTemplate: Callable[[str], Template | None],
) -> list[ChiralCenter]:
    """Restrain each center from findChiralCenters() to its template's volume.

    Return the restrained centers.
    """
    centers = findChiralCenters(newAtoms, getTemplate)
    if centers:
        force = mm.CustomCompoundBondForce(4, restraintEnergy)
        force.addPerBondParameter("ideal")
        for center in centers:
            force.addBond(center.particles, [center.ideal])
        system.addForce(force)
    return centers


def warnAboutChirality(
    centers: Iterable[ChiralCenter], positions: unit.Quantity, *, stacklevel: int
) -> None:
    """Warn, naming the residues, about centers whose positions don't match their templates.

    stacklevel has the same meaning as for warnings.warn(), counted from the caller of this
    function.
    """
    mismatched = [
        center.label for center in centers if not center.matchesTemplate(positions)
    ]
    if mismatched:
        warnings.warn(
            "Rebuilt centers whose handedness does not match their templates: "
            + ", ".join(mismatched),
            stacklevel=stacklevel + 1,
        )
