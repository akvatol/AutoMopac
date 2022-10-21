import pytest
from src.Basic.Atom import Atom
from src.Basic.symmetry_element import symmetry_element
from src.Symmetry.PointGroups import PointGroup


@pytest.mark.parametrize('PG, atom, expected', [
    (PointGroup(), A:= Atom.from_string('C 0 0 0'), {A: [A]}),
    (PointGroup(n=12), A:= Atom.from_string('C 2 0 0'), {A: [A]}),
    (PointGroup(n=12, axis='z'), A:= Atom.from_string('C 0 0 2'), {A: [A]}),
    (PointGroup(n=12, axis='y'), A:= Atom.from_string('C 0 2 0'), {A: [A]}),
    (PointGroup(n=2), A:= Atom.from_string('C 2 0 1'), {A: [A, Atom.from_string('C 2 0 -1')]}),
    (PointGroup(I=True), A:= Atom.from_string('C 2 0 1'), {A: [A, Atom.from_string('C -2 0 -1')]}),
    (PointGroup(n=2, I=True), A:= Atom.from_string('C 2 0 1'), {A: [A, Atom.from_string('C -2 0 -1'), Atom.from_string('C 2 0 -1'), Atom.from_string('C -2 0 1')]}),
    (PointGroup(n=2, I=True, h=True), A:= Atom.from_string('C 2 0 1'), {A: [A, Atom.from_string('C -2 0 -1'), Atom.from_string('C 2 0 -1'), Atom.from_string('C -2 0 1')]}),
    ])
def test_get_orbit(PG, atom, expected):
    assert len(PG.get_orbit(atom).get(atom)) == len(PG.get_orbit(atom).get(atom))
    assert set(PG.get_orbit(atom).get(atom)) == set(expected.get(atom))

