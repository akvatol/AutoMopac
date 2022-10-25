from attrs import define, field

import mpmath as mpm

mpm.mp.dps = 100

@define(slots=True)
class Atom:
    '''Interface for atoms for easy manipulating with structures
    '''
    atom: str
    coordinates: mpm.matrix = field(converter=mpm.matrix)
    asymmetric: bool = field(repr=True, default=True)
    _orbit: int = field(repr=False, default=None)

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

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other):
        # TODO: Docstring
        # TODO: Tests it
        if isinstance(other, Atom):
            return (self.atom == other.atom) and (self.coordinates == other.coordinates)
        else:
            raise TypeError('Atom object can be compared only with another Atom')