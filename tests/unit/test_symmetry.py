import pytest
from src.Symmetry.symmetry_element import symmetry_element
from src.Symmetry.utilites import point_group_symbol_parser
from src.Symmetry.utilites import make_cn
import mpmath as mpm
mpm.mp.dps = 100

I = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix(3, 1))

def test_symmetry_mul():
    assert make_cn(180)*make_cn(180) == make_cn(360)
    assert make_cn(90)*make_cn(180) == make_cn(270)
    assert  make_cn(90)*make_cn(120) != make_cn(220)

def test_symmetry_pow():
    assert make_cn(90)**3 == make_cn(270)
    assert make_cn(90)**4 == make_cn(0)
    assert make_cn(90)**2 == make_cn(90)*make_cn(90)
    assert make_cn(40)**0 == make_cn(360)
    assert make_cn(40)**3 == make_cn(60)**2
    assert make_cn(40)**3 != make_cn(40)**2

def test_find_order():
    assert I.find_order() == 2
    assert (make_cn(90)**3).find_order() == 4
    assert (make_cn(90)**2).find_order() == 2
    assert make_cn(160).find_order() == make_cn(40).find_order()
    assert make_cn(120).find_order() == 3
    for i in range(1, 50 + 1):
        assert make_cn(360/i).find_order() == i

def test_if_in():
    assert make_cn(180) in [make_cn(90), make_cn(90)**2]
    assert make_cn(0) in [make_cn(90)**i for i in range(4)]
    assert make_cn(160) in [make_cn(40)**i for i in range(16)]
    
def test_raises():
    with pytest.raises(ValueError):
        make_cn(90)*2

def test_point_group_pareser():
    assert point_group_symbol_parser('C1') == ('C', '1', '')
    assert point_group_symbol_parser('C12') == ('C', '12', '')
    assert point_group_symbol_parser('S2') == ('S', '2', '')
    assert point_group_symbol_parser('D12d') == ('D', '12', 'D')
    assert point_group_symbol_parser('X12s') == ('', '12', '')
    assert point_group_symbol_parser('D12dv') == ('D', '12', 'D')
    assert point_group_symbol_parser('DAG12WAS') == ('', '12', '')

def test_translation_eq():
    SE = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix([1, 0, 0]))
    SE2 = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix([2.5, 0, 0]))
    assert not SE.translation_eq(I)
    assert SE.translation_eq(SE, translation_vector=mpm.matrix([1, 0, 0]))
    assert SE.translation_eq(SE)
    assert not SE.translation_eq(symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix([2, 0, 0])), translation_vector=mpm.matrix([2, 0, 0]))
    assert not SE.translation_eq(symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix([4, 0, 0])), translation_vector=mpm.matrix([2, 0, 0]))
    assert SE2.translation_eq(symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix([5, 0, 0])), translation_vector=mpm.matrix([2.5, 0, 0]))