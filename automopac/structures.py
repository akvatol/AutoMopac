from abc import ABC, abstractmethod
from symmetry import SrewAxis, LineGroup, PointGroup

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
    def get_monomer():
        pass

    @abstractmethod
    def get_symcell(self):
        pass

    @abstractmethod
    def get_xyz(self):
        pass

    @abstractmethod
    def get_cell(self):
        pass

    @abstractmethod
    def get_group(self):
        pass

class Atom:
    '''Interface for atoms for easy manipulating with structures (???)
    TODO:
        * Понять, нужен ли он мне
    '''
    def __init__(self, atom, coords):
        pass

class Structure1D(Abstract_Structure):
    def __init__(self, symmetry: LineGroup, symcell: tuple) -> None:
        pass

    def get_xyz(self):
        pass

    def get_monomer(self):
        pass

    def get_symcell(self):
        pass

    def get_group(self):
        pass

    def get_cell():
        pass

class Molecule(Abstract_Structure):
    def __init__(self, symmetry: PointGroup, symcell: tuple) -> None:
        pass

    def get_xyz(self):
        pass

    def get_monomer(self):
        pass

    def get_symcell(self):
        pass

    def get_group(self):
        pass

    def get_cell():
        pass