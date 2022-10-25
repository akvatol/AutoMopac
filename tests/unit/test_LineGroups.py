import pytest
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis
from src.Symmetry.LineGroups import LineGroup

PGs = (PointGroup(), PointGroup(n=2), PointGroup(n=3), PointGroup(n=3, v=True))
SAs = (ScrewAxis(q=2, p=1, A=2), ScrewAxis(q=3, p=1, A=3), ScrewAxis(q=4, p=1, A=4))
LGs = (LineGroup(PG=PGs[0], SA=SAs[0]), LineGroup(PG=PGs[1], SA=SAs[1]), LineGroup(PG=PGs[2], SA=SAs[2]), LineGroup(PG=PGs[3], SA=SAs[2]))

@pytest.mark.parametrize('PG, SA, expected', [(PointGroup(n=1), ScrewAxis(q=2, p=1, A=2), 2),
                                              (PointGroup(n=1), ScrewAxis(q=3, p=1, A=3), 3),
                                              (PointGroup(n=2), ScrewAxis(q=2, p=1, A=2), 2),
                                              (PointGroup(n=3), ScrewAxis(q=3, p=1, A=3), 3)])
def test_LG_Q(PG, SA, expected):
    assert LineGroup(PG=PG, SA=SA).Q == expected

@pytest.mark.parametrize('LG, expected', [(LGs[0], dict(PG=dict(n=1, v=False, h=False, I=False, U=False, axis='x'), SA=dict(q=2, p=1, A=2, axis='x'))),
                                          (LGs[1], dict(PG=dict(n=2, v=False, h=False, I=False, U=False, axis='x'), SA=dict(q=3, p=1, A=3, axis='x'))),
                                          (LGs[2], dict(PG=dict(n=3, v=False, h=False, I=False, U=False, axis='x'), SA=dict(q=4, p=1, A=4, axis='x'))),
                                          (LGs[3], dict(PG=dict(n=3, v=True, h=False, I=False, U=False, axis='x'), SA=dict(q=4, p=1, A=4, axis='x')))])
def test_LG_to_dict(LG, expected):
    assert LG.to_dict() == expected