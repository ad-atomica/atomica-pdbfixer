from openmm.app.forcefield import ForceField
from openmm.app.topology import Residue

def matchResidueToTemplate(
    res: Residue,
    template: ForceField._TemplateData,
    bondedToAtom: list[set[int]],
    ignoreExternalBonds: bool = False,
    ignoreExtraParticles: bool = False,
) -> list[int] | None: ...
