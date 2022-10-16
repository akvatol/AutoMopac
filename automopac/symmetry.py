import re
import warnings
from abc import ABC, abstractmethod
from collections import namedtuple
from collections.abc import Iterable
from itertools import count

import mpmath as mpm
import numpy as np

from automopac.symmetry_element import symmetry_elemet

'''Модуль содержит в себе классы симметрии молекул и полимеров. 
Данные классы никак не взаимодействют со структурами, хранят группы и 
обеспечивают нормальное взаимодейтсвие групп между собой.
'''

#TODO: Добавить счетчик в метакласс

# TODO: Убрать отсюда
# I = symmetry_elemet(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]))
# sigma_xz = symmetry_elemet(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
# sigma_yz = symmetry_elemet(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))
# sigma_yx = symmetry_elemet(rotation=mpm.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
# C2x = symmetry_elemet(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
# C2y = symmetry_elemet(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]))
# C2z = symmetry_elemet(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]))

# http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Matrix_Symm_Op.htm
# x_set = {'sigma_v':sigma_xz, 'I':I, 'sigma_h':sigma_yz, 'C2_':C2y}
# y_set = {'sigma_v':sigma_yz, 'I':I, 'wsigma_h':sigma_xz, 'C2_':C2x}
# z_set = {'sigma_v':sigma_xz, 'I':I, 'sigma_h':sigma_yx, 'C2_':C2y}

# TODO: Доделать генераторы y установку




def all_LG1F_generator(q_max: int, n: int, Q_interval: tuple) -> tuple:
    """Generate all LG with q <= q_max, n = n, and Q_interval[0] <= Q <= Q_interval[1].

    Args:
        q_max (int): Maximal values of q_tilde (q_tilde = q/n)
        n (int): Index n of Cn group of monomer
        Q_interval (tuple): Q_interval - tuple with minimal and maximal Q

    Returns:
        tuple: list of sorted line groups (srew-axis if n = 1)
    """

    LineGroup1 = namedtuple("LineGroup1", "Q q p r n")

    # (r*p)%q = 1 original formula from Damnjanovich
    q_max = q_max + 1

    q = np.arange(0, q_max)
    r = np.arange(0, q_max - 1)
    # r = p
    data = []

    rp_mod_q = np.mod.outer(np.multiply.outer(r, r), q)
    indexes = np.nonzero(np.where(rp_mod_q == 1, rp_mod_q, 0))

    # q here - q_tilde, p - p_tilde
    for r, p, q in (zip(indexes[0], indexes[1], indexes[2])):
        if q < r or q < p:
            continue
        Q = q/r
        if Q >= Q_interval[0] and Q <= Q_interval[1]:
            data.append(LineGroup1(Q, q*n, p*n, r, n))

    return tuple(sorted(data, key=lambda x: x.Q))


def find_rp(q: int, p: int) -> int:
    """Return r or p depends on given arguments: if f(q, p) --> r else f(q, r) --> p;

    Args:
        q (int): q - is a q from line group theory (order of a screw-axis)
        p (int): in this case it may be r or p from LGT

    Returns:
        int: r or p depends on input
    """


    for r in range(1, 100):
        mod_q = r * p % q
        if mod_q == 1:
            break
    return int(r)


def point_group_symbol_parser(symbol: str) -> tuple:
    """Takes Schönflies symbol and return it`s readable representation

    Args:
        symbol (str): Schönflies symbol of point group

    Returns:
        tuple: Schönflies symbol representation
    """    

    regex = r'([CSD]?)(\d+)([VDH]?)|(C[IS])'
    symbol = symbol.upper()
    search = re.search(regex, symbol)
    if search:
        if search.group(4):
            return search.group(4), '', ''
        else:
            s1 = search.group(1)
            s2 = search.group(2)
            s3 = search.group(3)
            return s1, s2, s3


class AbstractGroup(ABC):
    '''
    Group - contain multiplication of all generators and generators
    '''

    def __init__(self):
        self._ids = count(0)
        pass

    @abstractmethod
    def get_generators(self):
        pass

    @abstractmethod
    def group(self):
        '''Return container with all group elements;
        '''
        pass

    @abstractmethod
    def check_group():
        '''Check if group is really a group'''
        pass



class PointGroup(AbstractGroup):
    '''Symmetry groups without translational part
    '''

    def __init__(self, generator, symbol:str=None):

        self.symbol = symbol

        if isinstance(generator, symmetry_elemet):
            self.generators = tuple([generator])
        elif isinstance(generator, Iterable) and all(isinstance(i, symmetry_elemet) for i in generator):
            self.generators = tuple(generator)
        else:
            raise ValueError(f'generator must be symmetry element or collections of symmetryelements {self.symbol}')

        if any(i.translation != mpm.matrix(3,1) for i in self.generators):
            raise ValueError(f'Elements of point group cannot contain translation part {self.symbol}')

    def get_generators(self):
        return self.generators
    
    def __generate_group(self):
        generators_elements = [i.get_all_powers() for i in self.generators]
        generators_elements.sort(key=len, reverse=True)
        group = []

        for num, subgroup in enumerate(generators_elements):
            if num == 0:
                group = frozenset(subgroup)
            else:
                group = frozenset((i*j for i in group for j in subgroup))

        return group

    @property
    def group(self):
        return self.__generate_group()

    def __eq__(self, other):
        if isinstance(other, PointGroup):
            return self.group == other.group

    def check_group(self):
        pass

    @classmethod
    def from_string(symbol: str = 'C1', axis: str = 'x'):
        # TODO Доделать этот метод инициализации, он удобен
        """Produce group from Schönflies symbol. Groups T, Th, Td, O, Oh, I, Ih not available."""
        # Parsing 1 letter, number, 3rd lettr
        # 1 letter - C, S, D
        # number - any
        # 3rd letter - None, v, h, d
        # Cn - Cn , Sn = I*Cn/2, Cv = Cn*sigma_v, Ch=Cn*sigma_h, Dn=Cn*C2_, Dnh = Dn*sigma_h, Dnd = Dn*sigma_v
        if symbol in ['T', 'Th', 'Td', 'O', 'Oh', 'I', 'Ih']:
            raise ValueError(
                'Groups T, Th, Td, O, Oh, I, Ih not available for generating from string')
        if axis not in 'xyz':
            raise ValueError('Axis value can be x, y or z only')
        # axis_dict = {'x':x_set, 'y':y_set, 'z':z_set}
        pass


class ScrewAxis(AbstractGroup):

    def __init__(self, q: int, p: int = None, r: int = None, A: int = 1, axis: str = 'x'):

        self.q = mpm.mpf(q)
        self.axis = axis
        self.A = mpm.mpf(A)
        self.angle = mpm.mpf(2*mpm.pi/q)

        if r and p:
            if find_rp(q, p) == find_rp(q, r):
                self.r = mpm.mpf(r)
                self.p = mpm.mpf(p)
            else:
                raise ValueError('p and r values does not match.')
        elif p or r:
            self.p = mpm.mpf(p) or mpm.mpf(find_rp(q, r))
            self.r = mpm.mpf(r) or mpm.mpf(find_rp(q, r))
        else:
            raise ValueError(
                'p or r value shoud be passed for screw-axis initialization.')

        self.Q = q/r

        # TODO: Убрать когда сделаю y установку
        if axis == 'y':
            raise ValueError('Y axis is not implemented yet')

        self.screw_generators = {'x': symmetry_elemet(rotation=mpm.matrix([
                                                        [1, 0, 0],
                                                        [0., mpm.cos(self.angle), -mpm.sin(self.angle)],
                                                        [0., mpm.sin(self.angle), mpm.cos(self.angle)]]),
                                    translation=mpm.matrix([A/q, 0, 0])),
                'y': '',
                'z': symmetry_elemet(rotation=mpm.matrix([[mpm.cos(self.angle), mpm.sin(self.angle), 0.],
                                                        [mpm.sin(self.angle), mpm.cos(self.angle), 0.],
                                                        [0, 0, 1]]),
                                    translation=mpm.matrix([A*p/q, 0, 0]))}

        self.generators = tuple([self.screw_generators[self.axis]])

    @property
    def group(self):
        return frozenset((self.generators[0] for i in range(self.q)))

    def get_generators(self):
        return self.generators

    def check_group(self):
        pass

class LineGroup(AbstractGroup):
    def __init__(self, screwaxis: ScrewAxis, pointgroup: PointGroup,  name: str = None):
        self.name = name
        self.generators = screwaxis.get_generators() + pointgroup.get_generators()

    def check_group(self):
        pass

    def get_generators(self):
        return self.generators

    def __generate_group(self):
        generators_elements = [i.get_all_powers() for i in self.generators]
        generators_elements.sort(key=len, reverse=True)
        group = []

        for num, subgroup in enumerate(generators_elements):
            if num == 0:
                group = frozenset(subgroup)
            else:
                group = frozenset((i*j for i in group for j in subgroup))

        return group

    @property
    def group(self):
        return self.__generate_group()