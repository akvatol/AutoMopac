import os
from pathlib import Path
import re

import mpmath as mpm
from src.QChem.Parsers.MopacParser import MopacParser
from src.QChem.Parsers.Template import ParserTemplate
from src.QChem.TaskRunner.TaskRunner import Process, execute_process
from src.Structure.Structure1D import Structure1D

from .Template import CalculatorTemplate


class MopacCalculator(CalculatorTemplate):
    
    def __init__(self, compound:Structure1D, project_dir:Path, parser: ParserTemplate):
        self.compound = compound
        self.root_dir = project_dir
        self.parser = parser

    def _make_template(self):
        template =['AUX LARGE CHARGE=0 SINGLET GEO-OK PM6-D3 HESS=10 NOREOR PDBOUT PRTXYZ  ITRY=500  NOSYM', '\n']
        for atom in self.compound:
            atom_ = atom.split()
            _ = f'{atom_[0]} {atom_[1]} * {atom_[2]} * {atom_[3]} *\n'
            template.append(_)
        A = self.compound.A
        tv = f'tv {str(A) + " *" if self.compound.SA.axis == "x" else 0} {str(A) + " *" if self.compound.SA.axis == "y" else 0} {str(A) + " *" if self.compound.SA.axis == "y" else 0}'
        template.append(tv)
        return ''.join(template)

    def save_input(self, template:str):
        angle = mpm.fdiv(360, self.compound.LG.SA.Q)
        angle_q_p = f'{angle}_{self.compound.LG.SA.q}_{self.compound.LG.SA.p}'
        # This calc dir
        new_path = self.root_dir / mpm.nstr(angle_q_p, n=8)

        # Создали папку
        if new_path.is_dir():
            pass
        else:
            os.mkdir(new_path)

        # Cохранили папку
        with open(new_path / 'input.mop', 'w') as fw:
            fw.write(template)

        return new_path

    def run_task(self, inputdir:Path, exec:str = 'mopac'):
        process = Process(args=[inputdir / 'input.mop'], executable=exec)
        return execute_process(process)

    def read_output(self, inputdir:Path):
        return self.parser.run(inputdir / 'input.out')

    def run(self, exec:str = 'mopac'):
        # Make intput
        calc_dir = self.save_input(self._make_template())
        # Run it
        self.run_task(inputdir=calc_dir, exec=exec)
        # Parse input
        data_scalars, data_arrays = self.read_output(calc_dir)
        return data_scalars, data_arrays
