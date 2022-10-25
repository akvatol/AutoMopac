import pytest

from src.Structure.Structure1D import Structure1D

from src.Symmetry.ScrewAxis import ScrewAxis
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.LineGroups import LineGroup

from src.Basic.Atom import Atom

PGs = (PointGroup(), PointGroup(n=2), PointGroup(n=3), PointGroup(n=3, v=True))
SAs = (ScrewAxis(q=2, p=1, A=2), ScrewAxis(q=3, p=1, A=3), ScrewAxis(q=4, p=1, A=4))
LGs = (LineGroup(PG=PGs[0], SA=SAs[0]), LineGroup(PG=PGs[1], SA=SAs[1]), LineGroup(PG=PGs[2], SA=SAs[2]), LineGroup(PG=PGs[3], SA=SAs[2]))

@pytest.mark.parametrize('struct, expected', [(Structure1D(LG=LGs[0], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0')]),
                                              (Structure1D(LG=LGs[1], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0')]),
                                              (Structure1D(LG=LGs[2], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0')]),
                                              (Structure1D(LG=LGs[3], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0')])])
def test_structure_monomer(struct, expected):
    assert list(struct.monomer) == expected


@pytest.mark.parametrize('struct, expected', [(Structure1D(LG=LGs[0], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0'), Atom.from_string('C 1 0 0')]),
                                              (Structure1D(LG=LGs[1], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0'), Atom.from_string('C 1 0 0'), Atom.from_string('C 2 0 0')]),
                                              (Structure1D(LG=LGs[2], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0'), Atom.from_string('C 1 0 0'), Atom.from_string('C 2 0 0'), Atom.from_string('C 3 0 0')]),
                                              (Structure1D(LG=LGs[3], symcell=[Atom.from_string('C 0 0 0')]), [Atom.from_string('C 0 0 0'), Atom.from_string('C 1 0 0'), Atom.from_string('C 2 0 0'), Atom.from_string('C 3 0 0')])])
def test_structure(struct, expected):
    assert list(struct.structure) == expected