import os
from pathlib import Path
import pysnooper as pnp

import mpmath as mpm
from src.QChem.Parsers.Template import ParserTemplate
import subprocess
from src.Structure.Structure1D import Structure1D

from .Template import CalculatorTemplate


class MopacCalculator(CalculatorTemplate):


    def __init__(self, project_dir:Path, parser: ParserTemplate):

        # TODO: Refactore this
        # Delete __init___ put everything in functions

        self.root_dir = project_dir
        self.parser = parser

    def _make_template(self, compound:Structure1D) -> str:
        template =['AUX LARGE CHARGE=0 SINGLET GEO-OK NOREOR HESS=0 PDBOUT PRTXYZ ITRY=500 AM1 NOSYM RECALC=10\n', '\n\n']
        for atom in compound.structure:
            atom_ = str(atom).split()
            _ = f'{atom_[0]} {atom_[1]} 1 {atom_[2]} 1 {atom_[3]} 1\n'
            template.append(_)
        A = compound.A
        tv = f'tv {str(A) + " 1" if compound.LG.SA.axis == "x" else "0.0 0"} {str(A) + " 1" if compound.LG.SA.axis == "y" else "0.0 0"} {str(A) + " 1" if compound.LG.SA.axis == "y" else "0.0 0"}\n'
        template.append(tv)
        return ''.join(template)

        # TODO: Убрать compound из аргументов
    def save_input(self, template:str, compound:Structure1D) -> Path:
        angle = mpm.nstr(mpm.fdiv(360, compound.LG.SA.Q), n=8)
        angle_q_p = f'{angle}_{compound.LG.SA.q}_{compound.LG.SA.p}'
        # This calc dir
        new_path = self.root_dir / angle_q_p

        # Создали папку
        if new_path.is_dir():
            pass
        else:
            os.mkdir(new_path)

        # Cохранили папку
        with open(new_path / 'input.mop', 'w') as fw:
            fw.write(template)

        return new_path # calc_dir

    def run_task(self, calc_dir:Path, exec:str = 'mopac'):
        subprocess.call([exec, calc_dir/'input.mop'])

    def read_output(self, calc_dir:Path):
        return self.parser.run(calc_dir / 'input.out')

    def run(self, compound:Structure1D, exec:str = 'mopac'):
        # Make intput
        compound = compound
        calc_dir = self.save_input(self._make_template(compound=compound), compound=compound)
        print(calc_dir)
        # Run it
        self.run_task(calc_dir=calc_dir, exec=exec)
        # Parse input
        data = self.read_output(calc_dir)
        return data
