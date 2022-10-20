import re
import warnings
from abc import ABC, abstractmethod
from collections import namedtuple
from collections.abc import Iterable
from itertools import count

import mpmath as mpm
import numpy as np

from src.Basic.symmetry_element import symmetry_elemet
from src.Symmetry.utilites import find_rp

# TODO: перевевсти
'''Модуль содержит в себе классы симметрии молекул и полимеров. 
Данные классы никак не взаимодействют со структурами, хранят группы и 
обеспечивают нормальное взаимодейтсвие групп между собой.
'''

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

        # G = x1*H + x2*H + ... + xn*H
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


class ScrewAxis(AbstractGroup):

    def __init__(self, q: int, p: int = 1, r: int = 1, A: int = 1, axis: str = 'x'):

        self.q = mpm.mpf(q)
        self.axis = axis
        self.angle = mpm.mpf(2*mpm.pi/q)

        # TODO: encapsulate it and make property (with setter and getter)
        self.A = mpm.mpf(A)

        if r and p:
            _r = find_rp(q, p)
            _p = find_rp(q, r)
            if (_r == r) and (_p == p):
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
    
    

