import re
from collections import namedtuple

import mpmath as mpm
import numpy as np
from src.Basic.Atom import Atom
from src.Basic.symmetry_element import symmetry_element
from dataclasses import dataclass

mpm.mp.mpds = 100

class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]
@dataclass(slots=True, frozen=True)
class SymmetryElements(metaclass=Singleton):
    # ? Возможно стоит перенести этот класс в отдельный файл
    I = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]))
    sigma_xz = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
    sigma_yz = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    sigma_yx = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
    C2x = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
    C2y = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]))
    C2z = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]))

    @staticmethod
    def make_cn_z(n: int) -> symmetry_element:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            symmetry_element: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [[mpm.cos(angle), -mpm.sin(angle), 0.0], 
            [mpm.sin(angle), mpm.cos(angle), 0.0],
            [0, 0, 1]],
        )
        trans = mpm.matrix(3, 1)
        return symmetry_element(rotation=rot, translation=trans)

    @staticmethod
    def make_cn_x(n: int) -> symmetry_element:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            symmetry_element: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [
                [1, 0, 0],
                [0.0, mpm.cos(angle), -mpm.sin(angle)],
                [0.0, mpm.sin(angle), mpm.cos(angle)],
            ]
        )
        trans = mpm.matrix(3, 1)
        return symmetry_element(rotation=rot, translation=trans)

    @staticmethod
    def make_cn_y(n: int) -> symmetry_element:
        """Support function for make Cn symmetry element

        Args:
            n (int): axis order

        Returns:
            symmetry_element: Cn group generator
        """
        angle = mpm.radians(360 / n)
        rot = mpm.matrix(
            [
                [mpm.cos(angle), 0, -mpm.sin(angle)],
                [0, 1, 0],
                [mpm.sin(angle), 0.0, mpm.cos(angle)],
            ]
        )
        trans = mpm.matrix(3, 1)
        return symmetry_element(rotation=rot, translation=trans)

# TODO: refactor this
# http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Matrix_Symm_Op.htm
PRESETS = {
    "x": {
        "v": SymmetryElements.sigma_xz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_yz,
        "U": SymmetryElements.C2y,
    },
    "y": {
        "v": SymmetryElements.sigma_yz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_xz,
        "U": SymmetryElements.C2z,
    },
    "z": {
        "v": SymmetryElements.sigma_xz,
        "I": SymmetryElements.I,
        "h": SymmetryElements.sigma_yx,
        "U": SymmetryElements.C2x
    },
}

def all_LG1F_generator(q_max: int, n: int, Q_interval: tuple) -> tuple:
    # TODO: Refactor this, delete namedtuples
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
    for r, p, q in zip(indexes[0], indexes[1], indexes[2]):
        if q < r or q < p:
            continue
        Q = q / r
        if Q >= Q_interval[0] and Q <= Q_interval[1]:
            data.append(LineGroup1(Q, q * n, p * n, r, n))

    return tuple(sorted(data, key=lambda x: x.Q))


def find_rp(q: int, p: int) -> int:
    """Return r or p depends on given arguments: if f(q, p) --> r else f(q, r) --> p;

    Args:
        q (int): q - is a q from line group theory (order of a screw-axis)
        p (int): in this case it may be r or p from LGT

    Returns:
        int: r or p depends on input
    """

    if (q == 1) and (p == 1):
        return 1

    for r in range(1, 100):
        mod_q = r * p % q
        if mod_q == 1:
            return int(r)


def point_group_symbol_parser(symbol: str) -> tuple:
    """Takes Schönflies symbol and return it`s readable representation

    Args:
        symbol (str): Schönflies symbol of point group

    Returns:
        tuple: Schönflies symbol representation
    """

    regex = r"([CSD]?)(\d+)([VDH]?)|(C[IS])"
    symbol = symbol.upper()
    search = re.search(regex, symbol)
    if search:
        if search.group(4):
            return search.group(4), "", ""
        else:
            s1 = search.group(1)
            s2 = search.group(2)
            s3 = search.group(3)
            return s1, s2, s3

def make_generators(parameters: dict) -> dict[str, symmetry_element]:
    # TODO: DOCSTRING
    """_summary_

    Returns:
        _type_: _description_
    """

    axis = parameters.get('axis', 'x')
    preset = PRESETS.get(axis)

    if axis == 'x':
        make_cn = SymmetryElements.make_cn_x(n=parameters.get('n', 1))
    elif axis == 'z':
        make_cn = SymmetryElements.make_cn_z(n=parameters.get('n', 1))
    else:
        make_cn = SymmetryElements.make_cn_y(n=parameters.get('n', 1))

    preset.update(dict(n=make_cn))

    generators_ = {name: preset.get(name) for name in parameters if parameters.get(name) and name != 'axis'}

    return generators_

def make_group(generators: dict) -> tuple:
    # TODO Docstring
    # TODO: regularize orders of group elements
    """_summary_

    Args:
        generators (frozenset): _description_

    Returns:
        frozenset: _description_
    """
    generators_elements = [generators[element].get_all_powers() for element in generators]
    generators_elements.sort(key=len, reverse=True)
    group = list(generators_elements[0])

    # G = x1*H + x2*H + ... + xn*H

    # optimize it
    for subgroup in generators_elements[1:]:
        for s_element in subgroup:
            if s_element in group:
                continue
            else:
                for element in group:
                    new_element = element*s_element
                    if new_element in group:
                        continue
                    else:
                        group.append(new_element)

    return tuple(group)

def detect_group(group: frozenset, axis: str) -> dict:
    # TODO: make it
    preset = PRESETS.get(axis)
    group_parameters = {}
    for element in preset:
        if element == 'n':
            pass

def _positive_validator(instance, attribute, value):
    # TODO; Docstring
    # Validator for point and line groups parameters
    if value <= 0:
        raise ValueError('Value vector must be positive')

def _make_srew_axis(q:int, p:int, A:float, axis:str = 'x') -> dict[str, symmetry_element]:
    # TODO: Docstring
    rotation_ = make_generators({'n':q, 'axis':axis})
    rotation = rotation_['n'].rotation
    translation_part = mpm.fmul(A, (mpm.fdiv(p, q)))
    screw_generators = {
        'x':symmetry_element(rotation=rotation, translation=mpm.matrix([translation_part, 0, 0]), translation_vector=mpm.matrix([A, 0, 0]),),
        'y':symmetry_element(rotation=rotation, translation=mpm.matrix([0, translation_part, 0]), translation_vector=mpm.matrix([0, A, 0]),),
        'z':symmetry_element(rotation=rotation, translation=mpm.matrix([0, 0, translation_part]), translation_vector=mpm.matrix([0, 0, A]),),
        }

    return {'q':screw_generators[axis]}