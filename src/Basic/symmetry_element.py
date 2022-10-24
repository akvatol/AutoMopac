"""Basic part of all symmetry groups.
"""

import warnings

import mpmath as mpm
from src.Basic.Atom import Atom

mpm.mp.mpds = 100
# Note that symmetry_element objects very sensitive to parameter mpm.mp.dps

# TODO: refactore this
class symmetry_element:
    """Matrix representation of symmetry element. Contain both translational and rotational part of symmetry element."""

    def __init__(self, rotation: mpm.matrix, translation: mpm.matrix = None, translation_vector: mpm.matrix = None):

        """Matrix representation of symmetry element

        Args:
            rotation (mpm.matrix): 3x3 matrix representation of rotational part of symmetry element.
            translation (mpm.matrix): translation part of symmetry element.
            translation_vector (mpm.matrix): translation vector of whole system. Needs to identify true (E|A) element in LineGroups.
        """

        if rotation.rows == rotation.cols and rotation.cols == 3:
            self.rotation = rotation
        else:
            raise ValueError(
                f"Rotational part of symmetry element shoud be mpmath.matrix 3x3 size instead {rotation}"
            )

        self.translation_vector = translation_vector or mpm.matrix(3, 1)

        if self.translation_vector.cols == 1 and self.translation_vector.rows == 3:
            pass
        else:
            raise ValueError(
                "Translational vector shoud be mpmath.matrix 3x1 size"
            )

        self.translation = translation or mpm.matrix(3, 1)
        if self.translation.cols == 1 and self.translation.rows == 3:
            pass
        else:
            raise ValueError(
                "Translational part of symmetry element shoud be mpmath.matrix 3x1 size"
            )

        self.translation = self.reduce_translation_part(self.translation, self.translation_vector)

    @property
    def order(self, limit=1000, tol=1e-15):
        """Return q if given, else try to find element order."""

        for i in range(1, limit + 1):
            delta = self.rotation**i - mpm.diag([1, 1, 1])
            if all(mpm.almosteq(delta[i, i], 0, abs_eps=tol) for i in range(2 + 1)):
                return i

        warnings.warn("Warning, cannot find element order, try to increase limit")
        return float("inf")

    def apply(self, atom: Atom) -> Atom:
        """Apply symmetry element to given atom."""
        return Atom(atom.atom, coordinates=mpm.chop((self.rotation * atom.coordinates) + self.translation))

    def rotation_eq(self, other):
        """Check if rotational part of two symmetry elements are equal."""

        delta_rot = self.rotation - other.rotation
        # check if all elements in delta_rot is almost equal to 0
        return all(
            mpm.almosteq(delta_rot[i, j], 0, abs_eps=10 ** -(mpm.mp.dps - 2))
            for i in range(2 + 1)
            for j in range(2 + 1))

    def translation_eq(self, other):
        """Check if translation part of two symmetry elements differ by translation vector"""
        if self.translation_vector == other.translation_vector:
            tp1 = self.reduce_translation_part(self.translation, self.translation_vector)
            tp2 = self.reduce_translation_part(other.translation, self.translation_vector)
            eq = tp1 == tp2
        else:
            eq = False
        return eq

    def get_all_powers(self) -> frozenset:
        n = self.order
        if n == float("inf"):
            raise ValueError(f"Order of element {self} to big or cannot be defined")
        return frozenset(self**i for i in range(1, n + 1))

    def __eq__(self, other):
        """Check if rotational part of symmetry elements is exacly the same. **May give wrong answer if elements contain translation**!!!"""
        if isinstance(other, symmetry_element):
            delta_rot = self.rotation - other.rotation
            rot_eq = self.rotation_eq(other)
            return rot_eq

    def __repr__(self):
        rot_ = "{} {} {}\n{} {} {}\n{} {} {}\n".format(
            *[mpm.nstr(mpm.chop(i, tol=1e-15)) for i in self.rotation]
        )
        rot_trans_ = rot_ + "{} {} {}".format(
            *[mpm.nstr(mpm.chop(i, tol=1e-15)) for i in self.translation]
        )
        return rot_trans_

    @staticmethod
    def reduce_translation_part(tp: mpm.matrix, tv: mpm.matrix):
        '''Уменьшает трансляционную часть на n*A, - где A - вектор трансляции сиситемы. Если вектор трансляции раввено 0, зануляет трансляционную часть.'''
        new_tp = mpm.matrix(3, 1)
        for i in range(3):
            if mpm.almosteq(tv[i], 0, abs_eps=1e-8):
                new_tp[i] = 0
            else:
                new_tp[i] = mpm.fsub(tp[i], mpm.fmul(tv[i], int(mpm.fdiv(tp[i], tv[i]))))

        return new_tp

    def __mul__(self, other):
        """Produce new symmetry element"""
        if isinstance(other, symmetry_element):
            mpm.mp.dps += 10
            new_rotation = self.rotation * other.rotation
            new_translation = self.rotation * other.translation + self.translation
            self.reduce_translation_part(new_translation, self.translation_vector)
            mpm.mp.dps -= 10
        else:
            raise ValueError(
                "Symmetry element can be multiplied only for another symmetry element"
            )
        return symmetry_element(rotation=new_rotation, translation=new_translation)

    def __hash__(self):
        return hash(str(self))

    def __pow__(self, other:int):
        """Return power of curent symmentry element.
        """

        new_rotation = self.rotation**other
        new_translation = self.translation * other
        return symmetry_element(rotation=new_rotation, translation=new_translation, translation_vector=self.translation_vector)


def main():
    SA = symmetry_element(rotation=mpm.matrix([[1, 0, 0],[0, -1, 0],[0, 0, -1]]), translation=mpm.matrix([5, 0, 0]), translation_vector=mpm.matrix([10, 0, 0]))
    print(SA)
    print(SA.get_all_powers())

if __name__ == '__main__':
    main()

