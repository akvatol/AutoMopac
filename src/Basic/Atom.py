from collections.abc import Iterable

import mpmath as mpm
mpm.mp.dps = 100

class Atom:
    '''Interface for atoms for easy manipulating with structures
    '''
    def __init__(self, atom:str, coords:Iterable):
        self.atom: str = atom
        self.coordinates: mpm.matrix = mpm.matrix([mpm.mpf(i) for i in coords])
        # вспомогательные атрибуты, могут понадобиться при работе с BigDFT
        self._fragment: int = None
        self._helix_num: int = None

    def __repr__(self):
        return f'{self.atom}\t{float(self.coordinates[0]):^20.12E} {float(self.coordinates[1]):^20.12E} {float(self.coordinates[2]):^20.12E}'

    @classmethod
    def from_string(cls, line):
        """Turn XYZ-type line into atoms

        Args:
            line (srt): _description_
        """
        atom, coords = line.split()[0], line.split()[1:]
        coords = tuple(mpm.mpf(i) for i in coords)
        return cls(atom=atom, coords=coords)