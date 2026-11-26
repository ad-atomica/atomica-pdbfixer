"""The OpenMM unit API used by PDBFixer."""

from collections.abc import Sequence

from openmm.vec3 import Vec3

from .quantity import Quantity as Quantity
from .unit import Unit as _Unit
from .unit_math import dot as dot
from .unit_math import norm as norm
from .unit_math import sqrt as sqrt
from .unit_math import sum as sum

# OpenMM attaches these operators to Unit when it imports unit_operators.
class Unit(_Unit):
    def __rmul__(self, value: float | Sequence[float | Vec3], /) -> Quantity: ...
    def __rtruediv__(self, value: float, /) -> Quantity: ...

# OpenMM generates the SI-prefixed names at import time.
nanometer: Unit
nanometers: Unit
picosecond: Unit
picoseconds: Unit
femtosecond: Unit
femtoseconds: Unit
angstrom: Unit
angstroms: Unit
kelvin: Unit
molar: Unit
kilocalories_per_mole: Unit
kilojoules_per_mole: Unit
