from abc import ABC, abstractmethod

class CalculatorTemplate(ABC):

    def write_input(self, path:str, template:str):
        with open(path, 'w') as fw:
            fw.write(template)

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