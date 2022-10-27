import re

from src.Basic.Atom import Atom

from .Template import ParserTemplate

class MopacParser(ParserTemplate):
    regex = {
        'method':r'\s+([\w\-\d]+)\s+[\w\s]+',
        'HoFKJ':r'\s+FINAL HEAT OF FORMATION\s+=\s+[-\d.]+\s+KCAL\/MOL\s+=\s+([-\d.]+)\s+KJ\/MOL',
        'IonizationPotentialEV':r'\s+[\D\s=]+([\d.]+)\s+[\D]+',
        'HomoLumo':r'\s+[\D\s=]\s+([-\d.]+)\s+([-\d.]+)',
        'MolecularWeightPoingGroup':r'\s+[\D\s=]+([\d.]+)\s+[\D\s]+:\s+(\w+\d+)',
        'xyz':r'\s*[\d]+\s+([\w]+)\([\w\d\s]+\)\s+([-\d\.]+)\s+\*?\s+([-\d\.]+)\s+\*?\s+([-\d\.]+)',
        'charge_mull':r'\s*([\d]+)\s+([\w]+)[\(\)\d\w\s]+\s+([-\d\.]+)\s+[\s\d\.]+'
    }

    def read_xyz_line(self, line):
        data = re.search(self.regex['xyz'], line)
        if data:
            atom, x, y, z = data.group(1), data.group(2), data.group(3), data.group(4)
            return Atom.from_string(f'{atom} {x} {y} {z}')

    def read_HoF(self, line):
        data = re.search(self.regex['HoFKJ'], line)
        value = None
        if data:
            value = data.group(1)
        return value

    def read_optimization(self):
        pass

    def read_charge(self, line):
        value = None
        data = re.search(self.regex['charge_mull'], line)
        if data:
            value = data.group(3)
        return value

    def read_method(self, line):
        value = None
        data = re.search(self.regex['method'], line)
        if data:
            value = data.group(1)
        return value

    def read_homo_lumo(self, line):
        homo, lumo = None, None
        data = re.search(self.regex['HomoLumo'], line)
        if data:
            homo, lumo = data.group(1), data.group(2)
        return homo, lumo

    def read_molecular_weight(self, line):
        weight, group = None, None
        data = re.search(self.regex['MolecularWeightPoingGroup'], line)
        if data:
            weight, group = data.group(1), data.group(2)
        return weight, group

    def read_IP(self, line):
        value = None
        data = re.search(self.regex['IonizationPotentialEV'], line)
        if data:
            value = data.group(1)
        return value

    def run(self, path):
        data = {}
        data_containers = {'structure':[], 'charges':[]}
        flag = 'method'
        with open(path, 'r') as fr:
            for line in fr:

                if 'Geometry optimization using' in line:
                    flag = 'cycles'
                if 'NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS' in line:
                    flag = 'charge'
                if 'SCF FIELD WAS ACHIEVED' in line:
                    flag = 'prop'
                if (flag == 'prop') and ('CHEMICAL' in line):
                    flag = 'structure'
                if 'ATOMIC ORBITAL ELECTRON POPULATIONS' in line:
                    flag = 'Orb'

                if flag == 'method':
                    if 'CALCULATION RESULTS' in line:
                        method = self.read_method(line)
                        data['method'] = method

                if flag == 'cycles':
                    if 'CYCLE:' in line:
                        continue

                if flag == 'structure':
                    atom = self.read_xyz_line(line)
                    if atom:
                        data_containers['structure'].append(atom)

                if flag == 'prop':
                    if 'FINAL HEAT OF FORMATION' in line:
                        value = self.read_HoF(line)
                        if value:
                            data['HEAT OF FORMATION'] = float(value)
                    if 'IONIZATION POTENTIAL' in line:
                        value = self.read_IP(line)
                        if value:
                            data['IONIZATION POTENTIAL'] = float(value)
                    if 'HOMO LUMO ENERGIES' in line:
                        homo, lumo = self.read_homo_lumo(line)
                        if homo and lumo:
                            data['homo'] = float(homo)
                            data['lumo'] = float(lumo)
                    if 'MOLECULAR WEIGHT' in line:
                        weight, group = self.read_molecular_weight(line)
                        if weight and group:
                            data['WEIGHT'] = float(weight)
                            data['Mopac Group'] = group
                if flag == 'charge':
                    charge = self.read_charge(line)
                    if charge:
                        data_containers['charges'].append(float(charge))
        return data, data_containers
