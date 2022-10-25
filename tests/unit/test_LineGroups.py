import pytest
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis
from src.Symmetry.LineGroups import LineGroup


@pytest.mark.parametrize('PG, SA, expected', [(PointGroup(n=1), ScrewAxis(q=2, p=1, A=2), 2),
                                              (PointGroup(n=1), ScrewAxis(q=3, p=1, A=3), 3),
                                              (PointGroup(n=2), ScrewAxis(q=2, p=1, A=2), 2),
                                              (PointGroup(n=3), ScrewAxis(q=3, p=1, A=3), 3)])
def test_LG_Q(PG, SA, expected):
    assert LineGroup(PG=PG, SA=SA).Q == expected