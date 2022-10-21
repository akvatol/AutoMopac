from abc import ABC, abstractmethod

class StructureTemplate(ABC):
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

class GroupTemplate:
    def __init__(self):
        pass

    @abstractmethod
    def group(self):
        pass

    @abstractmethod
    def generators(self):
        pass

    @abstractmethod
    def apply(self, atom):
        pass

    @abstractmethod
    def get_orbit(self, atom):
        pass

    @abstractmethod
    def get_stabilizer(self, atom):
        pass