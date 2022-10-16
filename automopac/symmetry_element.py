import warnings

import mpmath as mpm
mpm.mp.mpds = 100
# Note that symmetry_element objects very sensitive to parameter mpm.mp.dps

class symmetry_elemet:
    '''Matrix representation of symmetry element. Contain both translational and rotational part of symmetry element.
    '''
    __slots__ = ['rotation', 'translation']

    def __init__(self, rotation: mpm.matrix, translation: mpm.matrix = None):
        
        """Matrix representation of symmetry element

        Args:
            rotation (mpm.matrix): 3x3 matrix representation of rotational part of symmetry element
            translation (mpm.matrix): translation part of symmetry element
        """

        if rotation.rows == rotation.cols and rotation.cols == 3:
            self.rotation = rotation
        else:
            raise ValueError(
                f'Rotational part of symmetry element shoud be mpmath.matrix 3x3 size instead {rotation}')

        self.translation = translation or mpm.matrix(3, 1)
        if translation.cols == 1 and translation.rows == 3:
            pass
        else:
            raise ValueError(
                'Translational part of symmetry element shoud be mpmath.matrix 3x1 size')

    def transform_atom():
        pass

    def rotation_eq(self, other):
        """_summary_

        Args:
            other (symmetry_element): other symmetry element

        Returns:
            bool: shows if rotation part of two symmetry elements equal
        """
        delta_rot = self.rotation - other.rotation
        # check if all elements in delta_rot is almost equal to 0
        return all(mpm.almosteq(delta_rot[i, j], 0, abs_eps=10**-(mpm.mp.dps - 2)) for i in range(2 + 1) for j in range(2 + 1))

    def translation_eq(self, other, translation_vector=mpm.matrix(3, 1)):
        '''Check if translation part of two symmetry elements differ by translation vector
        '''
        # T = n*A + x_0 -> 0 = n*A + x_0 - T -> 0 = n*A + delta
        delta_trans = self.translation - other.translation
        eq_list = []
        for i, j in zip(translation_vector, delta_trans):
            
            # Если трансляционная часть отличается на вектор трансляции -> True
            try:
                # Если findroot не может найти корни, то выдаёт ошибку
                cond_ = mpm.almosteq(mpm.findroot(lambda n: n*i - j, 1), int(mpm.findroot(lambda n: n*i - j, 1)))
            except:
                cond_ = False

            # Если трансляционная часть равна -> True
            if abs(i) == abs(j):
                eq_list.append(True)
            elif cond_:
                eq_list.append(True)
            else:
                eq_list.append(False)
            
        return True if all(eq_list) else False

    def find_order(self, limit=1000, tol=1e-15, q:int=None):
        '''Return q if given, else try to find element order.
        '''
        
        # foolproof, infinite element can be here
        if q:
            return q

        for i in range(1, limit + 1):
            delta = self.rotation**i - mpm.diag([1, 1, 1])
            if all(mpm.almosteq(delta[i, i], 0, abs_eps=tol) for i in range(2 + 1)):
                return i

        warnings.warn(
            'Warning, cannot find element order, try to increase limit')
        return float('inf')

    def get_all_powers(self):
        n = self.find_order()
        if n == float('inf'):
            raise ValueError(f'Order of element {self} to big or cannot be defined')
        return tuple(self ** i for i in range(1, n + 1))

    def __eq__(self, other):
        '''Check if rotational part of symmetry elements is exacly the same. **May give wrong answer if elements contain translation**!!!
        '''
        if isinstance(other, symmetry_elemet):
            delta_rot = self.rotation - other.rotation
            rot_eq = all(mpm.almosteq(delta_rot[i, j], 0, abs_eps=10**-(
                mpm.mp.dps - 2)) for i in range(2 + 1) for j in range(2 + 1))
            return rot_eq

    def __repr__(self):
        rot_ = '{} {} {}\n{} {} {}\n{} {} {}\n'.format(
            *[mpm.nstr(mpm.chop(i, tol=1e-30)) for i in self.rotation])
        rot_trans_ = rot_ + \
            '{} {} {}'.format(*[mpm.nstr(mpm.chop(i, tol=1e-30))
                              for i in self.translation])
        return rot_trans_

    def __mul__(self, other):
        '''Produce new symmetry element
        '''
        if isinstance(other, symmetry_elemet):
            mpm.mp.dps += 10
            new_rotation = self.rotation*other.rotation
            new_translation = self.rotation*other.translation + self.translation
            mpm.mp.dps -= 10
        else:
            raise ValueError(
                'Symmetry element can be multiplied only for another symmetry element')
        return symmetry_elemet(rotation=new_rotation, translation=new_translation)

    def __hash__(self):
        return hash(str(self))

    def __pow__(self, other):
        new_rotation = self.rotation**other
        new_translation = self.translation*other
        return symmetry_elemet(rotation=new_rotation, translation=new_translation)

