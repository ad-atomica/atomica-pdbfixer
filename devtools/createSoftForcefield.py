"""
createSoftForceField.py: Creates a force field XML file that is suitable
for removing clashes in very badly strained systems.

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2013 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

__author__ = "Peter Eastman"
__version__ = "1.0"

import openmm.app.element as elem
import openmm.app.forcefield as ff
from openmm import app

forcefield = app.ForceField("amber99sbildn.xml", "tip3p.xml")
bondK = 10000.0
angleK = 10.0

# Create the new force field file.

print("<ForceField>")

# Print the atom types, while identifying types and classes to omit.

print(" <AtomTypes>")
omitTypes = set()
omitClasses = set()
for atomType in forcefield._atomTypes:
    data = forcefield._atomTypes[atomType]
    if data.element is None or data.element == elem.hydrogen:
        omitTypes.add(atomType)
        omitClasses.add(data.atomClass)
    else:
        print(
            f'  <Type name="{atomType}" class="{data.atomClass}" element="{data.element.symbol}" mass="{data.mass:g}"/>'
        )
print(" </AtomTypes>")

# Print the residue templates.

print(" <Residues>")
skipTemplates = ("ASH", "CYM", "GLH", "HID", "HIE", "LYN")
for template in forcefield._templates.values():
    if template.name in skipTemplates or template.name[1:] in skipTemplates:
        continue
    print(f'  <Residue name="{template.name}">')
    atomIndex = {}
    for i, atom in enumerate(template.atoms):
        if atom.type not in omitTypes:
            print(f'   <Atom name="{atom.name}" type="{atom.type}"/>')
            atomIndex[i] = len(atomIndex)
    for a1, a2 in template.bonds:
        if a1 in atomIndex and a2 in atomIndex:
            print(f'   <Bond from="{atomIndex[a1]}" to="{atomIndex[a2]}"/>')
    for atom in template.externalBonds:
        if atom in atomIndex:
            print(f'   <ExternalBond from="{atomIndex[atom]}"/>')
    print("  </Residue>")
print(" </Residues>")

# Print the harmonic bonds.

print(" <HarmonicBondForce>")
bonds = next(f for f in forcefield._forces if isinstance(f, ff.HarmonicBondGenerator))
for i in range(len(bonds.types1)):
    type1 = next(iter(bonds.types1[i]))
    type2 = next(iter(bonds.types2[i]))
    if type1 not in omitTypes and type2 not in omitTypes:
        class1 = forcefield._atomTypes[type1].atomClass
        class2 = forcefield._atomTypes[type2].atomClass
        print(
            f'  <Bond class1="{class1}" class2="{class2}" length="{bonds.length[i]:g}" k="{bondK:g}"/>'
        )
print(" </HarmonicBondForce>")

# Print the harmonic angles.

print(" <HarmonicAngleForce>")
angles = next(f for f in forcefield._forces if isinstance(f, ff.HarmonicAngleGenerator))
for i in range(len(angles.types1)):
    type1 = next(iter(angles.types1[i]))
    type2 = next(iter(angles.types2[i]))
    type3 = next(iter(angles.types3[i]))
    if type1 not in omitTypes and type2 not in omitTypes and type3 not in omitTypes:
        class1 = forcefield._atomTypes[type1].atomClass
        class2 = forcefield._atomTypes[type2].atomClass
        class3 = forcefield._atomTypes[type3].atomClass
        print(
            f'  <Angle class1="{class1}" class2="{class2}" class3="{class3}" angle="{angles.angle[i]:g}" k="{angleK:g}"/>'
        )
print(" </HarmonicAngleForce>")

# Print the periodic torsions.

print(" <PeriodicTorsionForce>")
torsions = next(
    f for f in forcefield._forces if isinstance(f, ff.PeriodicTorsionGenerator)
)
wildcardType = forcefield._atomClasses[""]
for torsion in torsions.proper:
    type = [
        next(iter(torsion.types1)),
        next(iter(torsion.types2)),
        next(iter(torsion.types3)),
        next(iter(torsion.types4)),
    ]
    wildcard = [
        torsion.types1 == wildcardType,
        torsion.types2 == wildcardType,
        torsion.types3 == wildcardType,
        torsion.types4 == wildcardType,
    ]
    if all(type[i] not in omitTypes or wildcard[i] for i in range(4)):
        classes = tuple(
            "" if wildcard[i] else forcefield._atomTypes[type[i]].atomClass
            for i in range(4)
        )
        print(
            '  <Proper class1="{}" class2="{}" class3="{}" class4="{}"'.format(
                *classes
            ),
            end=" ",
        )
        for i in range(len(torsion.k)):
            print(
                f' periodicity{i + 1}="{torsion.periodicity[i]}" phase{i + 1}="{torsion.phase[i]:g}" k{i + 1}="{torsion.k[i]:g}"',
                end=" ",
            )
        print("/>")
for torsion in torsions.improper:
    type = [
        next(iter(torsion.types1)),
        next(iter(torsion.types2)),
        next(iter(torsion.types3)),
        next(iter(torsion.types4)),
    ]
    wildcard = [
        torsion.types1 == wildcardType,
        torsion.types2 == wildcardType,
        torsion.types3 == wildcardType,
        torsion.types4 == wildcardType,
    ]
    if all(type[i] not in omitTypes or wildcard[i] for i in range(4)):
        classes = tuple(
            "" if wildcard[i] else forcefield._atomTypes[type[i]].atomClass
            for i in range(4)
        )
        print(
            '  <Improper class1="{}" class2="{}" class3="{}" class4="{}"'.format(
                *classes
            ),
            end=" ",
        )
        for i in range(len(torsion.k)):
            print(
                f' periodicity{i + 1}="{torsion.periodicity[i]}" phase{i + 1}="{torsion.phase[i]:g}" k{i + 1}="{torsion.k[i]:g}"',
                end=" ",
            )
        print("/>")
print(" </PeriodicTorsionForce>")

# Print the script to add the soft-core nonbonded force.

print(" <Script>")
print("""import openmm as mm
nb = mm.CustomNonbondedForce('C/((r/0.2)^4+1)')
nb.addGlobalParameter('C', 1.0)
sys.addForce(nb)
for i in range(sys.getNumParticles()):
    nb.addParticle([])
exclusions = set()
for bond in data.bonds:
    exclusions.add((min(bond.atom1, bond.atom2), max(bond.atom1, bond.atom2)))
for angle in data.angles:
    exclusions.add((min(angle[0], angle[2]), max(angle[0], angle[2])))
for a1, a2 in exclusions:
    nb.addExclusion(a1, a2)""")
print(" </Script>")

print("</ForceField>")
