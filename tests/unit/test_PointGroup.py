from tokenize import group
import mpmath as mpm
import pytest
from src.Basic.symmetry_element import symmetry_element
from src.Symmetry.PointGroups import PointGroup
from src.Basic.utilites import SymmetryElements

make_cn = SymmetryElements.make_cn_x

@pytest.mark.parametrize(
    "group, expected",
    [
        (PointGroup(), 1),
        (PointGroup(n=2), 2),
        (PointGroup(n=3), 3),
        (PointGroup(n=1, I=True, h=True), 4),
        (PointGroup(n=3, v=True), 6),
        (PointGroup(n=3, U=True), 6),
        (PointGroup(n=3, U=True, h=True), 12)
    ],
)
def test_group_len(group, expected):
    assert len(group.group) == expected

@pytest.mark.parametrize(
    "group, expected",
    [
        (PointGroup(), {'n': make_cn(1)}),
        (PointGroup(n=2), {'n': make_cn(2)}),
        (PointGroup(n=3), {'n': make_cn(3)}),
        (PointGroup(n=3, v=True),
         {'n': make_cn(3),
          'v':symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))}),
        (PointGroup(n=3, U=True, h=True),
         {'n': make_cn(3),
          'h': symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])),
          'U':symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]))}),
    ],
)
def test_generators(group, expected):
    assert group.generators == expected


@pytest.mark.parametrize(
    "group, expected",
    [
        (PointGroup(), {'n': 1, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=2), {'n': 2, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3), {'n': 3, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3, v=True), {'n': 3, 'v': True, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3, U=True), {'n': 3, 'v': False, 'h': False, 'U': True, 'I':False}),
        (PointGroup(n=3, U=True, h=True), {'n': 3, 'v': False, 'h': True, 'U': True, 'I':False}),
    ],
)
def test_to_dict(group, expected):
    assert group.to_dict() == expected

@pytest.mark.parametrize(
    "expected, group",
    [
        (PointGroup(), {'n': 1, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=2), {'n': 2, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3), {'n': 3, 'v': False, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3, v=True), {'n': 3, 'v': True, 'h': False, 'U': False, 'I':False}),
        (PointGroup(n=3, U=True), {'n': 3, 'v': False, 'h': False, 'U': True, 'I':False}),
        (PointGroup(n=3, U=True, h=True), {'n': 3, 'v': False, 'h': True, 'U': True, 'I':False}),
    ],
)
def test_from_dict(group, expected):
    assert PointGroup.from_dict(group) == expected

def test_popgen():
    assert PointGroup(n=3, v=True).popgen('n') == (make_cn(3), PointGroup(v=True))
    assert PointGroup(n=3, U=True).popgen('U') == (symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])), PointGroup(n=3))
    assert PointGroup(n=4, v=True).popgen('v') == (symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]])), PointGroup(n=4))
    assert PointGroup(n=5, h=True).popgen('h') == (symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])), PointGroup(n=5))
    assert PointGroup(n=6, I=True).popgen('I') == (symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])), PointGroup(n=6))
    assert PointGroup(n=3, U=True, h=True).popgen('h') == (symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])), PointGroup(n=3, U=True))

@pytest.mark.parametrize(
    "group, expected",
    [
        (PointGroup(), []),
        (PointGroup(n=2), []),
        (PointGroup(n=3, v=True), []),
        (PointGroup(n=6, I=True, h=True), ['h', 'I']),
        (PointGroup(n=6, U=True, I=False, h=False, v=True), []),
        (PointGroup(n=6, U=True, h=True, I=False, v=False), []),
        (PointGroup(n=6, U=True, h=True, v=True, I=False), ['U', 'h', 'v']),
        (PointGroup(n=6, I=True, h=True, U=True, v=True), ['h', 'I', 'U', 'v']),
        (PointGroup(n=1, I=True, h=True, U=True, v=True), ['n', 'I', 'U', 'v']),
    ],
)
@pytest.mark.skip(reason='There is a bug')
def test_find_extra_generators(group, expected):
    assert set(group.find_extra_generators()) == set(expected)