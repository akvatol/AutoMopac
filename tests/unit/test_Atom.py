import pytest
from src.Basic.Atom import Atom
import mpmath as mpm
mpm.mp.dps = 100

@pytest.mark.parametrize("string, expected", [
    ('C 1 2 3', Atom(atom='C', coordinates=[1, 2, 3])),
    ('1 1 2 3', Atom(atom='1', coordinates=[1, 2, 3])),
    ('256 1 2 3', Atom(atom='256', coordinates=[1, 2, 3]))
    ])
def test_atom_from_strig(string, expected):
    assert Atom.from_string(string) == expected

@pytest.mark.parametrize("atom1, atom2, expected", [(Atom.from_string('C 0 0 0'), Atom.from_string('C 0 0 0'), True),
                                                    (Atom.from_string('C 0 0 1'), Atom.from_string('C 0 0 0'), False),
                                                    (Atom.from_string('N 0 0 0'), Atom.from_string('C 0 0 0'), False),
                                                    (Atom.from_string('12 0 0 0'), Atom.from_string('12 0 0 0'), True),
                                                    (Atom.from_string('C -10 0 10'), Atom.from_string('C 10 0 -10'), False)])
def test_atom_eq(atom1, atom2, expected):
    value = atom1 == atom2
    assert value == expected