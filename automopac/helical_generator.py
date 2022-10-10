import numpy as np
from collections import namedtuple

from abc import ABC, abstractmethod

def all_LG1F_generator(q_max:int, n:int , Q_interval:tuple):
    ''' Generate all LG with q <= q_max, n = n, and Q_interval[0] <= Q <= Q_interval[1].

    Args:
        q_tilde_max: int
            Maximal values of q_tilde (q_tilde = q/n)
        n: int
            Index n of Cn group of monomer
        Q_interval: list
            Q_interval - tuple with minimal and maximal Q

    Returns:
        tuple with sorted line groups (first family)
    '''
    
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
        if Q >=Q_interval[0] and Q<=Q_interval[1]:
            data.append(LineGroup1(Q, q*n, p*n, r, n))

    return tuple(sorted(data, key=lambda x: x.Q))
