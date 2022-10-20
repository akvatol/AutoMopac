import re
from collections import namedtuple

import mpmath as mpm
import numpy as np
from src.Basic.Atom import Atom
from src.Basic.symmetry_element import symmetry_element

mpm.mp.mpds = 100

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


# TODO: Make it better 
I = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]))
sigma_xz = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
sigma_yz = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))
sigma_yx = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
C2x = symmetry_element(rotation=mpm.matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
C2y = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]))
C2z = symmetry_element(rotation=mpm.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]))

# TODO: This code smells
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
            [mpm.cos(angle), -mpm.sin(angle), 0.0],
            [0, 1, 0],
            [0.0, mpm.sin(angle), mpm.cos(angle)],
        ]
    )
    trans = mpm.matrix(3, 1)
    return symmetry_element(rotation=rot, translation=trans)

def make_generators(parameters: dict) -> dict[str:symmetry_element]:
    """_summary_

    Returns:
        _type_: _description_
    """

    # TODO: refactor this
    # http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Matrix_Symm_Op.htm
    presets = {
        "x": {
            "v": sigma_xz,
            "I": I,
            "h": sigma_yz,
            "U": C2y,
            "n": make_cn_x,
        },
        "y": {
            "v": sigma_yz,
            "I": I,
            "h": sigma_xz,
            "U": C2z,
            "n": make_cn_y,
        },
        "z": {
            "v": sigma_xz,
            "I": I,
            "h": sigma_yx,
            "U": C2x,
            "n": make_cn_z,
        },
    }

    preset = presets.get("axis", "x")

    generators = {name:preset[name] for name in parameters if parameters.get(name)}

    return generators

def make_group(generators: frozenset) -> frozenset:
    # TODO Написать
    """_summary_

    Args:
        generators (frozenset): _description_

    Returns:
        frozenset: _description_
    """
    generators_elements = [generators[element].get_all_powers() for element in generators]
    generators_elements.sort(key=len, reverse=True)
    group = []

    # G = x1*H + x2*H + ... + xn*H
    for num, subgroup in enumerate(generators_elements):
        if num == 0:
            group = frozenset(subgroup)
        else:
            group = frozenset((i*j for i in group for j in subgroup))

    return group
