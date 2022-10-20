import pytest
from src.Basic.symmetry_element import symmetry_element
from src.Symmetry.utilites import point_group_symbol_parser
from src.Symmetry.utilites import make_cn_z as make_cn
import mpmath as mpm
mpm.mp.dps = 100

I = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), translation=mpm.matrix(3, 1))

def test_symmetry_element_mul():
    assert make_cn(2)*make_cn(2) == make_cn(1)
    assert  make_cn(4)*make_cn(3) != make_cn(4)
    assert make_cn(4)*make_cn(2)*make_cn(4) == make_cn(1)

def test_symmetry_element_pow():
    assert make_cn(4)**4 == make_cn(1)
    assert make_cn(4)**2 == make_cn(2)
    assert make_cn(9)**0 == make_cn(1)
    assert make_cn(9)**3 == make_cn(6)**2
    assert make_cn(9)**3 != make_cn(9)**2
    assert make_cn(180)**2 == make_cn(90)
    assert make_cn(360)**2 == make_cn(180)

def test__symmetry_element_order():
    assert I.order == 2
    assert (make_cn(4)**3).order == 4
    assert (make_cn(4)**2).order == 2
    assert make_cn(9).order == (make_cn(18)*make_cn(18)).order
    assert make_cn(3).order == 3
    for i in range(1, 50 + 1):
        assert make_cn(i).order == i

def test_if_symmetry_element_in():
    assert make_cn(2) in [make_cn(4), make_cn(4)**2]
    assert make_cn(1) in [make_cn(4)**i for i in range(4)]
    assert make_cn(9) in [make_cn(9)**i for i in range(16)]
    assert make_cn(90) in make_cn(180).get_all_powers()
    
def test_symmetry_element__raises():
    with pytest.raises(ValueError):
        make_cn(4)*2

    with pytest.raises(ValueError):
        symmetry_element(rotation=mpm.matrix(4))

    with pytest.raises(ValueError):
        symmetry_element(rotation=mpm.matrix(3), translation_vector=mpm.matrix([1, 2, 3, 4]))

    with pytest.raises(ValueError):
        symmetry_element(rotation=mpm.matrix(3), translation=mpm.matrix([1, 2, 3, 4]))

@pytest.mark.skip(reason='Function does not needed anymore')
def test_point_group_pareser():
    assert point_group_symbol_parser('C1') == ('C', '1', '')
    assert point_group_symbol_parser('C12') == ('C', '12', '')
    assert point_group_symbol_parser('S2') == ('S', '2', '')
    assert point_group_symbol_parser('D12d') == ('D', '12', 'D')
    assert point_group_symbol_parser('X12s') == ('', '12', '')
    assert point_group_symbol_parser('D12dv') == ('D', '12', 'D')
    assert point_group_symbol_parser('DAG12WAS') == ('', '12', '')

@pytest.mark.parametrize(
    'SE1,SE2,expected',
    [
        (make_cn(2), make_cn(3), True),
        (symmetry_element(mpm.matrix(3), mpm.matrix([7, 0, 0]), mpm.matrix([8, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([7, 0, 0]), mpm.matrix([5, 0, 0])),
         False),
        (symmetry_element(mpm.matrix(3), mpm.matrix([7, 0, 0]), mpm.matrix([8, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([7, 0, 0]), mpm.matrix([8, 0, 0])),
         True),
        (symmetry_element(mpm.matrix(3), mpm.matrix([8, 0, 0]), mpm.matrix([2, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([0, 0, 0]), mpm.matrix([5, 0, 0])),
         False),
        (symmetry_element(mpm.matrix(3), mpm.matrix([0, 0, 0]), mpm.matrix([0, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([0, 0, 0]), mpm.matrix([0, 0, 0])),
         True),
        (symmetry_element(mpm.matrix(3), mpm.matrix([8, 0, 0]), mpm.matrix([4, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([0, 0, 0]), mpm.matrix([4, 0, 0])),
         True),
        (symmetry_element(mpm.matrix(3), mpm.matrix([6, 0, 0]), mpm.matrix([4, 0, 0])),
         symmetry_element(mpm.matrix(3), mpm.matrix([2, 0, 0]), mpm.matrix([4, 0, 0])),
         True),
     ]
)
def test_translation_eq(SE1, SE2, expected):
    assert SE1.translation_eq(SE2) == expected


@pytest.mark.parametrize('tp,tv,expected', [((0, 0, 0), (0, 0, 0), (0, 0, 0)),
                                            ((1, 0, 0), (0, 0, 0), (0, 0 , 0)),
                                            ((1, 0, 0), (1, 0, 0), (0, 0, 0)),
                                            ((3, 0, 0), (1, 0, 0), (0, 0, 0)),
                                            ((1.5, 0, 0), (1, 0, 0), (.5, 0, 0)),
                                            ((1, 0, 0), (2, 0, 0), (1, 0, 0))])
def test_reduce_translation_part(tp, tv, expected):
    tp = mpm.matrix(tp)
    tv = mpm.matrix(tv)
    expected = mpm.matrix(expected)
    assert symmetry_element.reduce_translation_part(tp, tv)