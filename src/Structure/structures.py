from abc import ABC, abstractmethod
from collections.abc import Iterable

import mpmath as mpm

from symmetry import LineGroup, PointGroup, SrewAxis, symmetry_elemet


class Abstract_Structure(ABC):
    '''Prototype for all structure types (Abstract factory)

    TODO:
        * Prototype for all types of structures (3D, 2D, 1D, 0D)
        * Contains information about structure
        * Interface to manipulate with structure (Supercell, monomer scrolling, N_merization)
    '''
    def __init__(self):
        pass

    @abstractmethod
    def monomer():
        pass

    @abstractmethod
    def symcell(self):
        pass

    @abstractmethod
    def xyz(self):
        pass

    @abstractmethod
    def cell(self):
        pass

    @abstractmethod
    def group(self):
        pass

class Atom:
    '''Interface for atoms for easy manipulating with structures
    '''
    def __init__(self, atom:str, coords:Iterable):
        self.atom: str = atom
        self.coordinates: tuple[float] = tuple(mpm.mpf(i) for i in coords)
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

class Structure1DKernel(Abstract_Structure):
    # Класс содержит всю логику работы с одномерными объектами
    def __init__(self):
        pass

class Structure1D(Structure1DKernel):
    # Интерфейс к ядру 
    def __init__(self, symmetry: LineGroup, symcell: tuple) -> None:
        pass

    @property
    def xyz(self):
        pass

    def monomer(self):
        pass

    def symcell(self):
        pass

    def group(self):
        pass

    def cell(self):
        pass

class Helix:
    def __init__(self, structure: Structure1D ,Q_min: float, Q_max: float, q_max: int = None, r_max: int = None):
        pass

class Molecule(Abstract_Structure):
    def __init__(self, symmetry: PointGroup, symcell: tuple) -> None:
        self.group = symmetry.group
        self.symcell = symcell 

    @property
    def xyz(self):
        pass

    def monomer(self):
        pass

    def symcell(self):
        pass

    def group(self):
        pass

    def cell(self):
        pass
