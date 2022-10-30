from attrs import field, define
from src.Basic.Atom import Atom
from src.Basic.Templates import GroupTemplate
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis


@define(slots=True)
class LineGroupBase(GroupTemplate):
    PG: PointGroup = field(kw_only=True)
    SA: ScrewAxis = field(kw_only=True)

    __generators: dict = field(init=False, repr=False)
    __group: tuple = field(init=False, repr=False)

    @property
    def Q(self):
        # TODO:
        return self.SA.Q

    def __repr__(self):
        return f'{self.PG}\n{self.SA}\n'
    
    @__generators.default
    def _make_generators(self):
        # TODO:
        d = self.PG.generators
        d.update(self.SA.generators)
        return d

    def _make_group(self) -> tuple:
        # TODO: regularize orders of group elements
        # * Умножаю все элементы винтовой оси на все элементы точечной группы
        data = []
        for i in self.SA.group:
            for j in self.PG.group:
                new_symm = i*j
                if new_symm in data:
                    pass
                else:
                    data.append(new_symm)
        return tuple(data)

    @property
    def generators(self):
        return self.__generators

    @property
    def group(self):
        self.__group = self.__group or self._make_group()
        return self.__group

    @property
    def A(self):
        return self.SA.A

    @property
    def Q(self):
        return self.SA.Q

    @property
    def f(self):
        return self.SA.f

class LineGroup(LineGroupBase):

    def apply(self, atoms:list[Atom]):
        # TODO: Docstring
        # TODO: Test it
        monomer = self.PG.apply(atoms)
        structure = self.SA.apply(monomer)
        return structure

    def get_orbit(self, atom:Atom) -> dict:
        # TODO: Docstring
        # TODO: Test it

        orbit = []
        for num, elements in enumerate(self.group):

            # if num == 0, it is identity opperation, so it is original atom
            # if num != 0, it is false, because it change atom
            asymmetric = not num

            new_atom = elements.apply(atom)
            new_atom.asymmetric = asymmetric

            if new_atom in orbit:
                pass
            else:
                orbit.append(new_atom)

        return {atom:orbit}

    def to_dict(self):
        # TODO:
        return dict(PG=self.PG.to_dict(), SA=self.SA.to_dict())

    @classmethod
    def from_dict(cls, parameters):
        # TODO:
        return cls(SA=ScrewAxis.from_dict(parameters['SA']), PG=PointGroup.from_dict(parameters['PG']))

    def get_stabilizer(self, atom:Atom) -> PointGroup:
        # TODO: Make it
        # TODO: Test it
        # * Элементы винтовой оси не могут быть стабилизатором
        return self.PG.get_stabilizer(self, atom)

    def reduce(self, type:str):
        # TODO: Make it
        # TODO: Test it
        pass

    def reduce_monomer_symmetry(self, generator: str):
        # TODO: Make it
        # TODO: Test it
        pass

    def reduce_screw_axis(self, index: str):
        # TODO: Make it
        # TODO: Test it
        pass

    def reduce_glide_plane(self):
        # TODO: Make it
        # TODO: Test it
        pass