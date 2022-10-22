from abc import ABC, abstractmethod

class StructureTemplate(ABC):
    '''Prototype for all structure types.
        * Prototype for 1D structures
        * Interface to manipulate with structure
    '''
    def __init__(self):
        pass

    @abstractmethod
    def monomer(self):
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
    def reduce_monomer_symmetry(self):
        pass

    @abstractmethod
    def reduce_screw_axis(self):
        pass

    @abstractmethod
    def analyze_orbits(self):
        pass

class GroupTemplate:
    """Prototype for all group types.
    """
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

    @abstractmethod
    def copy(self):
        pass