import pytest
from src.Basic.Atom import Atom
from src.Basic.symmetry_element import symmetry_element
from src.Symmetry.ScrewAxis import ScrewAxis
import mpmath as mpm

mpm.mp.dps = 100

SAs = (ScrewAxis(q=1, p=1, A=10),
       ScrewAxis(q=2, p=1, A=10),
       ScrewAxis(q=3, p=1, A=12),
       ScrewAxis(q=17, p=6, A=0.100754840323E+02),
       ScrewAxis(q=28, p=17, A=0.167892743723E+02))

@pytest.mark.parametrize("SA, expected", [(dict(q=1, p=1, A=10), SAs[0]),
                                          (dict(q=2, p=1, A=10), SAs[1]),
                                          (dict(q=3, p=1, A=12), SAs[2]),
                                          (dict(q=17, p=6, A=0.100754840323E+02), SAs[3]),
                                          (dict(q=28, p=17, A=0.167892743723E+02), SAs[4])])
def test_from_dict(SA, expected):
    assert ScrewAxis.from_dict(SA) == expected

@pytest.mark.parametrize("SA,expected", [(SAs[0], dict(q=1, p=1, A=10)),
                                         (SAs[1], dict(q=2, p=1, A=10)),
                                         (SAs[2], dict(q=3, p=1, A=12)),
                                         (SAs[3], dict(q=17, p=6, A=0.100754840323E+02)),
                                         (SAs[4], dict(q=28, p=17, A=0.167892743723E+02))])
def test_to_dict(SA, expected):
    assert SA.to_dict() == expected

@pytest.mark.parametrize("SA,expected", [(SAs[0], 1),
                                         (SAs[1], 2),
                                         (SAs[2], 3),
                                         (SAs[3], 17),
                                         (SAs[4], 28)])
def test_groups_len(SA, expected):
    assert len(SA.group) == expected

@pytest.mark.parametrize("SA,expected", [(SAs[1], 2),
                                         (SAs[1], 2)])
def test_SA(SA, expected):
    assert SA.generators['q'] == symmetry_element(rotation=mpm.matrix([[1, 0, 0],[0, -1, 0],[0, 0, -1]]), translation=mpm.matrix([5, 0, 0]), translation_vector=mpm.matrix([10, 0, 0]))

@pytest.mark.parametrize("SA,expected", [(SAs[0], 0),
                                         (SAs[1], 5),
                                         (SAs[2], 4),
                                         (SAs[3], mpm.fdiv(mpm.fmul(0.100754840323E+02, 6), 17)),
                                         (SAs[4], mpm.fdiv(mpm.fmul(0.167892743723E+02, 17), 28))])
def test_generators_translation(SA, expected):
    assert mpm.almosteq(SA.generators['q'].translation[0] - expected, 0, abs_eps=1e-15)