from ..Basic.Templates import GroupTemplate
from src.Basic.symmetry_element import symmetry_elemet
from src.Symmetry.utilites import make_generators, make_group

class PointGroup(GroupTemplate):

    def __init__(self, n: int = 1, v:bool = False, h:bool = False, I: bool = False, U: bool = False, axis: str = 'x'):
        self.__n = n
        self.__v = v
        self.__h = h
        self.__I = I
        self.__U = U
        self.axis = axis


    @property
    def n(self):
        return self.__n

    @property
    def v(self):
        return self.__v

    @property
    def h(self):
        return self.__h

    @property
    def I(self):
        return self.__I

    @property
    def U(self):
        return self.__U

    @property
    def generators(self):
        return make_generators(n=self.n, v=self.v, h=self.h, I=self.I, U=self.U, axis=self.axis)


    @property
    def group(self):
        return make_group(self.generators)

    def apply(self, atoms: tuple) -> tuple:
        """Apply all symmetry elements for set of atoms.

        Args:
            atoms (tuple): set of atoms for transformation

        Returns:
            tuple: generated structue (N1, N1`, N1``, ..., N2, N2`, N2``, ...)
        """
        strucure = frozenset([SE.apply(atom) for SE in self.group for atom in atoms])
        return strucure

    def find_extra_generators(self):
        # TODO
        old_group = self.group
        
        pass

    def popgen(self):
        #TODO
        pass