from abc import ABC, abstractmethod

class CalculatorTemplate(ABC):

# TODO: Разобраться почему не работает
#    def __new__(cls):
#        if not hasattr(cls, 'instance'):
#             cls.instance = super(CalculatorTemplate, cls).__new__(cls)
#        return cls.instance

    @abstractmethod
    def _make_template(self):
        pass

    @abstractmethod
    def run_task(self):
        pass

    @abstractmethod
    def read_output(self):
        pass

    @abstractmethod
    def run(self):
        pass