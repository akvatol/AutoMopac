from abc import ABC, abstractmethod

class ParserTemplate(ABC):

    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ParserTemplate, cls).__new__(cls)
        return cls.instance

    @abstractmethod
    def read_xyz_line(self):
        pass

    @abstractmethod
    def read_optimization(self):
        pass

    @abstractmethod
    def read_charge(self):
        pass

    @abstractmethod
    def read_method(self):
        pass

    @abstractmethod
    def read_homo_lumo(self):
        pass

    @abstractmethod
    def read_molecular_weight(self):
        pass

    @abstractmethod
    def run(self):
        pass