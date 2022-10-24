import pytest
import src.Basic.utilites as ut

@pytest.mark.parametrize('q, p, r', [(1, 1, 1), (2, 1, 1), (6, 3, None), (6, 2, None)])
def test_find_rp(q, p, r):
    assert ut.find_rp(q, p) == r