from attrs import field, frozen
from src.Basic.symmetry_element import symmetry_element

from src.Basic.Templates import GroupTemplate
from src.Basic.Atom import Atom
from src.Basic.utilites import _make_srew_axis, _positive_validator, find_rp
import mpmath as mpm

mpm.mp.dps = 100


@frozen(slots=True)
class ScrewAxisBase(GroupTemplate):
    # TODO: Docstring
    q: int = field(kw_only=True, default=1, validator=_positive_validator)
    p: int = field(kw_only=True, default=1, validator=_positive_validator)
    A: int = field(kw_only=True, validator=_positive_validator)
    axis: str = field(kw_only=True, default='x')
    r: int = field(init=False)
    Q: int = field(init=False)
    
    __generators: dict = field(init=False, repr=False)
    __group: frozenset = field(init=False, repr=False)

    @p.validator
    def _p_validation(instance, attribute, value):
        if instance.p > instance.q:
            raise ValueError('p value cannot be greater than q')

    @r.default
    def _find_r(self):
        r = find_rp(self.q, self.p)
        if r:
            return r
        else:
            raise ValueError('q and p must be coprime numbers')

    @Q.default
    def _set_Q(self) -> float:
        return self.q/self.r

    @__generators.default
    def _srew_axis_generator(self) -> symmetry_element:
        return _make_srew_axis(q=self.q, p=self.p, A=self.A, axis=self.axis)

    @__generators.validator
    def _screw_axis_validation(instance, attribute, value):
        if value.get('q').order != instance.q:
            raise ValueError(f'Generated screw axis does not constent\n q = {instance.q}, screw axis order = {value.order}')
    
    @__group.default
    def _make_screw_axis_group(self) -> tuple:
        return self.__generators['q'].get_all_powers()

    @property
    def group(self):
        return self.__group

    @property
    def generators(self):
        return self.__generators


class ScrewAxis(ScrewAxisBase):

    def reduce(self, index: int):
        # TODO: Docstring
        # TODO: Test it
        if index < self.Q:
            # TODO: Перевести
            raise ValueError("Группа не может быть меньше чем порядок оси Q")
        elif index == self.Q:
            SA = ScrewAxis(q=1, p=1, A=self.A)
        else:
            SA = self._reduce_screw_axis(index)
            SA = ScrewAxis(q=SA[0], p=SA[1], A=self.A, )
        return SA

    def _reduce_screw_axis(self, index: int ) -> list[int]:
        # TODO: Docstring
        if self.q % index == 0:
            new_q = self.q / index
        else:
            new_q = self.q * index
            new_p = find_rp(self.q, self.r)
        return new_q, new_p

    def apply(self, atoms:list[Atom]) -> tuple:
        # TODO: Docstring
        # TODO: Checck if Q = 1 SA works fine
        # TODO: Test it
        # * Screw-axis always change the atom
        structure = []
        for _order, SE in enumerate(self.group):
            for atom in atoms:
                new_atom = SE.apply(atom)

                # На всякий случай:( Проверка на дурака
                if new_atom in structure:
                    continue

                if _order == 0:
                    asymmetric = atom.asymmetric
                else:
                    asymmetric = False

                new_atom = atom._orbit
                new_atom.asymmetric = asymmetric

                structure.append(new_atom)

        return tuple(structure)

    def get_stabilizer(self, atom:Atom):
        # TODO: Docstring
        # TODO: Test it
        # * Stabilizer of screw axis is L1 Group
        # ? Do i need that?
        return ScrewAxis(q=1, p=1, A=self.A)

    def get_orbit(self, atom: Atom):
        """Возвращает словарь вида {Atom:[Atom, Atom1, Atom2 ..., AtomN]}, где Atom1, Atom2 ..., AtomN получены из Atom преобразованияями симметрии

        Returns:
            dict[Atom:list[Atom]]: Словарь содержащий все орбиты 
        """
        # TODO: Test it
        orbit = {atom:[]}
        for elements in self.group:
            new_atom = elements.apply(atom)
            if new_atom in orbit[atom]:
                pass
            else:
                orbit[atom].append(new_atom)

        return orbit

    def to_dict(self):
        # TODO: Docstring
        # TODO: Test it
        return {'q': self.q, 'p':self.p, 'A':self.A}

    @classmethod
    def from_dict(cls, parameters:dict[str, int]):
        # TODO: Docstring
        # TODO: Test it
        return cls(q=parameters.get('q', 1), p=parameters.get('p', 1), A=parameters.get('A'), axis=parameters.get('axis', 'x'))

    def copy(self):
        # TODO: Docstring
        # TODO: Test it
        return self.from_dict(self.to_dict())