from hashlib import new
from attrs import field, frozen

from src.Basic.Atom import Atom
from src.Basic.SymmetryElement import SymmetryElement
from src.Basic.Templates import GroupTemplate
from src.Basic.utilites import detect_group, make_generators, make_group


@frozen(slots=True)
class PointGroupBase(GroupTemplate):
    # TODO: Документация
    n: int = field(default=1, kw_only=True)
    v: bool = field(default=False)
    h: bool = field(default=False)
    I: bool = field(default=False)
    U: bool = field(default=False)
    axis: str = field(default='x', kw_only=True)
    __generators: dict = field(init=False, repr=False)
    # Оставил этот параметр т.к. не влият на производительность. Мы всегда используем всего 1 точечную группу
    __group: tuple[SymmetryElement] = field(init=False, repr=False)

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

class PointGroup(PointGroupBase):

    def get_orbit(self, atom: Atom) -> list[Atom]:
        # TODO: Docstring

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

        return orbit

    def get_stabilizer(self, atom: Atom):
        # Определяем число элементов в локальной группе = len(self.group)/len(self.get_orbit(atom).get(atom)), если 1 то E, если N то исходная группа
        # Определяем те элементы группы которые не меняют положение атома и сохраняем их в список
        # Ищем в нём каждый элемент из пресета (xyz) и удоляем если находим:
        # Сперва ищем элементы с индексом 2
        # То что осталось - группа Cn -> len(...) = n
        stabilizer_index = len(self.group)/len(self.get_orbit(atom).get(atom))

        #TODO: refactore this
        # TODO: Test it
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
            # TODO: Сделать чтобы работало
            subgroup_parameters = detect_group(subgroup_elements, self.axis)
            subgroup = self.from_dict(subgroup_parameters)
        return subgroup

    def apply(self, atoms: list[Atom]) -> tuple[Atom]:
        """Apply all symmetry elements for set of atoms.

        Args:
            atoms (tuple): set of atoms for transformation

        Returns:
            tuple: generated structue (N1, N1`, N1``, ..., N2, N2`, N2``, ...)
        """
        # TODO: Test it
        structure = []

        for _orbit, atom in enumerate(atoms):
            for _order, SE in enumerate(self.group):
                new_atom = SE.apply(atom)

                # Если атом в списке переходим к следующему шагу
                # Если нет - продолжаем
                if new_atom in structure:
                    continue

                # * Присваиваем номер орбиты и определяем ассиметричен ли атом
                # * Если атом получен единичным элементом, то он ассиметричен
                new_atom.asymmetric = not _order
                new_atom._orbit = _orbit + 1
                structure.append(new_atom)

        return tuple(structure)

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

    def copy(self):
        # TODO: Docstring
        return self.from_dict(self.to_dict())

    def to_dict(self):
        # TODO: Docstring
        return dict(n=self.n, I=self.I, U=self.U, v=self.v, h=self.h, axis=self.axis)

    @classmethod
    def from_dict(cls, parameter: dict):
        # TODO: Docstring
        return cls(
            n=parameter.get("n", 1),
            I=parameter.get("I", False),
            U=parameter.get("U", False),
            v=parameter.get("v", False),
            h=parameter.get("h", False),
            axis=parameter.get("axis", "x"),
        )
