from src.Basic.symmetry_element import symmetry_element
from ..Basic.Templates import GroupTemplate
from ..Basic.Atom import Atom
from src.Basic.utilites import make_generators, make_group, detect_group
from attr import frozen, field

@frozen(slots=True)
class PointGroup(GroupTemplate):
    n: int = field(default=1, kw_only=True)
    v: bool = field(default=False, kw_only=True)
    h: bool = field(default=False, kw_only=True)
    I: bool = field(default=False, kw_only=True)
    U: bool = field(default=False, kw_only=True)
    axis: str = field(default='x', kw_only=True)
    __generators: dict = field(init=False, repr=False)
    __group: frozenset[symmetry_element] = field(init=False, repr=False)

    # TODO: Validation for I, v, h, U

    @__generators.default
    def _make_generators(self):
        return make_generators(dict(n=self.n, v=self.v, h=self.h, I=self.I, U=self.U, axis=self.axis))

    @__group.default
    def _make_group(self):
        return make_group(self.generators)

    @property
    def generators(self):
        return self.__generators

    @property
    def group(self):
        return self.__group

    @classmethod
    def from_dict(cls, parameter: dict):
        return cls(
            n=parameter.get("n", 1),
            I=parameter.get("I", False),
            U=parameter.get("U", False),
            v=parameter.get("v", False),
            h=parameter.get("h", False),
        )

    def to_dict(self):
        return dict(n=self.n, I=self.I, U=self.U, v=self.v, h=self.h)

    def copy(self):
        return self.from_dict(self.to_dict())

    def get_orbit(self, atom: Atom) -> dict[Atom:list[Atom]]:
        """Возвращает словарь вида {Atom:[Atom, Atom1, Atom2 ..., AtomN]}, где Atom1, Atom2 ..., AtomN получены из Atom преобразованияями симметрии

        Returns:
            dict[Atom:list[Atom]]: Словарь содержащий все орбиты 
        """
        orbit = {atom:[]}
        for elements in self.group:
            new_atom = elements.apply(atom)
            if new_atom in orbit[atom]:
                pass
            else:
                orbit[atom].append(new_atom)

        return orbit

    def get_stabilizer(self, atom: Atom):
        # Определяем число элементов в локальной группе = len(self.group)/len(self.get_orbit(atom).get(atom)), если 1 то E, если N то исходная группа
        # Определяем те элементы группы которые не меняют положение атома и сохраняем их в список
        # Ищем в нём каждый элемент из пресета (xyz) и удоляем если находим:
        # Сперва ищем элементы с индексом 2
        # То что осталось - группа Cn -> len(...) = n
        stabilizer_index = len(self.group)/len(self.get_orbit(atom).get(atom))

        #TODO: refactore this
        if stabilizer_index == len(self.group):
        # Это действие - копирование. 
            subgroup = self.copy()
        elif stabilizer_index == 1:
            # C1 group
            subgroup = self.from_dict({})
        else:
            subgroup_elements = []
            for element in self.group:
                if element.apply(Atom) == Atom:
                    subgroup_elements.append(element)
                else:
                    pass
            subgroup_parameters = detect_group(subgroup_elements, self.axis)
            subgroup = self.from_dict(subgroup_parameters)
        return subgroup

    def apply(self, atoms: tuple) -> frozenset:
        """Apply all symmetry elements for set of atoms.

        Args:
            atoms (tuple): set of atoms for transformation

        Returns:
            tuple: generated structue (N1, N1`, N1``, ..., N2, N2`, N2``, ...)
        """
        strucure = frozenset([SE.apply(atom) for SE in self.group for atom in atoms])
        return strucure

    def popgen(self, gen_name: str):
        """Delete generator and return (deleted generator, new PointGroup)

        Args:
            gen_name (str): Name of generator for deleting
        """

        generator = self.generators.get(gen_name)

        if generator:
            if gen_name == 'n':
                new_gen_value = 1
            else:
                new_gen_value = False

            parameters = self.to_dict()
            parameters[gen_name] = new_gen_value

        return generator, self.from_dict(parameter=parameters)