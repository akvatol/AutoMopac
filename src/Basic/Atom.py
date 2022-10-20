from dataclasses import dataclass, field

import mpmath as mpm

mpm.mp.dps = 100

@dataclass(slots=True)
class Atom:
    '''Interface for atoms for easy manipulating with structures
    '''
    atom: str
    coordinates: mpm.matrix
    _fragment: int = field(repr=False, default=None)
    _helix_num: int = field(repr=False, default=None)

    def __post_init__(self):
        self.coordinates = mpm.matrix(self.coordinates)

    def __repr__(self):
        return f'{self.atom}\t{float(self.coordinates[0]):^20.12E} {float(self.coordinates[1]):^20.12E} {float(self.coordinates[2]):^20.12E}'

    @classmethod
    def from_string(cls, line: str):
        """Turn XYZ-type line into atoms

        Args:
            line (str): Name x y z - format string
        """
        atom, coords = line.split()[0], line.split()[1:]
        coords = mpm.matrix([mpm.mpf(i) for i in coords])
        return cls(atom=atom, coordinates=coords)