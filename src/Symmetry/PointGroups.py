from ..Basic.Templates import GroupTemplate
from src.Symmetry.utilites import make_generators, make_group
from dataclasses import dataclass, field

@dataclass(slots=True)
class PointGroup(GroupTemplate):
    n: int = 1
    v: bool = False
    h: bool = False
    I: bool = False
    U: bool = False
    axis: str = "x"
    generators: dict = field(init=False, repr=False)
    group: frozenset = field(init=False, repr=False)

    def __post_init__(self):
        self.generators = make_generators(dict(n=self.n, v=self.v, h=self.h, I=self.I, U=self.U, axis=self.axis))
        self.group = make_group(self.generators)

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

    def apply(self, atoms: tuple) -> tuple:
        """Apply all symmetry elements for set of atoms.

        Args:
            atoms (tuple): set of atoms for transformation

        Returns:
            tuple: generated structue (N1, N1`, N1``, ..., N2, N2`, N2``, ...)
        """
        strucure = frozenset([SE.apply(atom) for SE in self.group for atom in atoms])
        return strucure

    def find_extra_generators(self) -> list:
        extra_generators = []
        old_generators = self.generators

        for generator in old_generators:

            if old_generators.get(generator):

                newer_generators = old_generators.copy()
                del newer_generators[generator]
                newer_group = make_group(newer_generators)

                if newer_group == self.group:
                    extra_generators.append(generator)

        return extra_generators

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