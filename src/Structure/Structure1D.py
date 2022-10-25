import json

from attrs import define, field
from src.Basic.Atom import Atom
from src.Symmetry.LineGroups import LineGroup
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis


@define(kw_only=True)
class Structure1DBase:
    # TODO: Docstring
    # TODO: Test it

    LG: LineGroup
    symcell: tuple[Atom] = field(converter=tuple)

    status: str = field(default=None)

    properties: dict[str, float] = field(repr=False, factory=dict)
    containers: dict[str, tuple[float]] = field(repr=False, factory=dict)

    monomer_C1: tuple[Atom] = field(repr=False, default=None)
    structure_noPG: tuple[Atom] = field(repr=False, default=None)
    structure_C1: tuple[Atom] = field(repr=False, default=None)

    @property
    def A(self):
        return self.LG.A

    @property
    def monomer(self):
        return self.LG.PG.apply(self.symcell)

    @property
    def structure(self):
        return self.LG.apply(self.symcell)

    @property
    def monomer_atoms(self):
        return [self.structure.index(atom) for atom in self.monomer]

    @property
    def symcell_atoms(self):
        return [self.structure.index(atom) for atom in self.symcell]


class Structure1D(Structure1DBase):

    def get_orbits(self) -> dict[int, Atom]:
        # TODO: Docstring
        # TODO: Test it

        orbits = {}
        for atom in self.structure:
            if orbits.get(atom._orbit):
                orbits[atom._orbit].append(atom)
            else:
                orbits[atom._orbit] = []
                orbits[atom._orbit].append(atom)

        return orbits

    def to_xyz(self, path):
        # TODO:
        structure = self.structure
        with open(path, 'w') as wf:
            wf.write(f'{len(structure)}\n\n')
            for atom in structure:
                wf.write(str(atom) + '\n')

            x = self.A if self.LG.SA.axis == 'x' else 0
            y = self.A if self.LG.SA.axis == 'y' else 0
            z = self.A if self.LG.SA.axis == 'z' else 0

            wf.write(f'tv\t{x:^20.12E} {y:^20.12E} {z:^20.12E}')

    def to_dict(self):
        # TODO:
        d = dict(LG=self.LG.to_dict(),
                 symcell=[str(atom) for atom in self.symcell],
                 status=self.status,
                 properties=self.properties,
                 containers=self.containers,
                 monomer_C1=self.monomer_C1,
                 structure_noPG=self.structure_noPG,
                 structure_C1=self.structure_C1)
        return d

    @classmethod
    def from_dict(cls, parameters:dict):
        # TODO:
        _ = cls(LG=LineGroup.from_dict(parameters.get('LG')),
                symcell=[Atom.from_string(atom) for atom in parameters.get('symcell', [None]) if atom],
                status=parameters.get('status'),
                properties=parameters.get('properties'),
                containers=parameters.get('containers'),
                monomer_C1=parameters.get('monomer_C1'),
                structure_noPG=parameters.get('structure_noPG'),
                structure_C1=parameters.get('structure_C1'))
        return _

    def to_json(self, path):
        # TODO:
        with open(path, 'w') as wf:
            json.dump(fp=wf, obj=self.to_dict())

    @classmethod
    def from_json(cls, path):
        # TODO:
        with open(path, 'r') as fr:
            d = json.load(fp=fr)
        return cls.from_dict(d)

    def to_f34(self, path):
        # TODO:
        pass
