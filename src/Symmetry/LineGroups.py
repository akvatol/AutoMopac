from attrs import define, frozen
from src.Basic.Atom import Atom
from src.Basic.Templates import GroupTemplate
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis


@frozen(slots=True)
class LineGroupBase(GroupTemplate):
    PG: PointGroup = define(kw_only=True)
    SA: ScrewAxis = define(kw_only=True)

    __generators: dict = define(init=False, repr=False)
    __group: dict = define(init=False, repr=False)

    @property
    def Q(self):
        return self.SA.Q

    def __repr__(self):
        return f'{self.PG}\n{self.SA}\n'
    
    @__generators.default
    def _make_generators(self):
        d = self.PG.generators
        d.update(self.SA.generators)
        return d

    @__group.gefault
    def _make_group(self) -> frozenset:
        # TODO: regularize orders of group elements
        # * Умножаю все элементы винтовой оси на все элементы точечной группы
        data = []
        for i in self.SA.group:
            for j in self.PG.group:
                new_symm = i*j
                if new_symm in data:
                    data.append(new_symm)
                else:
                    pass
        return frozenset(data)

    @property
    def generators(self):
        return self.__generators

    @property
    def group(self):
        return self.__group


class LineGroup(LineGroupBase):

    def apply(self, atoms:list[Atom]):
        # TODO: Docstring
        # TODO: Test it
        structure = frozenset([SE.apply(atom) for SE in self.group for atom in atoms])
        return structure

    def get_orbit(self, atom:Atom) -> dict:
        # TODO: Docstring
        # TODO: Test it
        # * Орбита только точечаня группа
        orbit = {atom:[]}
        for elements in self.group:
            new_atom = elements.apply(atom)
            if new_atom in orbit[atom]:
                pass
            else:
                orbit[atom].append(new_atom)

        return orbit

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