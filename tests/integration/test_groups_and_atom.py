import pytest
from src.Basic.Atom import Atom
from src.Basic.symmetry_element import symmetry_element
from src.Symmetry.PointGroups import PointGroup
from src.Symmetry.ScrewAxis import ScrewAxis


SAs = (ScrewAxis(q=1, p=1, A=10),
       ScrewAxis(q=2, p=1, A=10),
       ScrewAxis(q=3, p=1, A=12),
       ScrewAxis(q=17, p=6, A=0.100754840323E+02),
       ScrewAxis(q=28, p=17, A=0.167892743723E+02))

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
def test_PG_get_orbit(PG, atom, expected):
    assert set(PG.get_orbit(atom)) == set(expected.get(atom))

@pytest.mark.parametrize("SA,atom,expected", [
    (SAs[0], A:= Atom.from_string('C 2 0 1'), {A:[A]}),
    (SAs[1], A:= Atom.from_string('C 2 0 1'), {A:[A, Atom.from_string('C 7 0 -1')]}),
    (SAs[2], A:= Atom.from_string('C 2 0 0'), {A:[A, Atom.from_string('C 6 0 0'), Atom.from_string('C 10 0 0')]}),
])
def test_SA_get_orbit(SA, atom, expected):
    assert set(SA.get_orbit(atom)) == set(expected.get(atom))