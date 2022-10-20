from abc import ABC, abstractmethod

from src.Symmetry import LineGroup, PointGroup, SrewAxis, symmetry_elemet




class Structure1DKernel(Abstract_Structure):
    # Класс содержит всю логику работы с одномерными объектами
    def __init__(self):
        pass

class Structure1D(Structure1DKernel):
    # Интерфейс к ядру 
    def __init__(self, symmetry: LineGroup, symcell: tuple) -> None:
        pass

    @property
    def xyz(self):
        pass

    def monomer(self):
        pass

    def symcell(self):
        pass

    def group(self):
        pass

    def cell(self):
        pass

class Helix:
    def __init__(self, structure: Structure1D ,Q_min: float, Q_max: float, q_max: int = None, r_max: int = None):
        pass