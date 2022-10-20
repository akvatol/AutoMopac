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