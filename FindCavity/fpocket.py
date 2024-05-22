from . import overlap
import os
import pandas as pd
import re

class Fpocket:
    def __init__(self, protein_folder, fpocket_folder):
        self.pockets = self.get_pockets(protein_folder, fpocket_folder)

    @staticmethod
    def find_cavity(ligands, fpocket_cavities):
        """
        checkes for every ligand the overlap with the cavities to return the right cavity for the ligand.
        :param ligands: list of pandas dataframes with ligand coordinates
        :param fpocket_cavities: folder with all output cavities from fpocket
        :return: corresponding cavity to ligand
        """
        cavities = os.listdir(fpocket_cavities)
        cavities = [f for f in cavities if f.endswith('.pdb')]
        results = dict()
        count=1
        for i in ligands:
            results[count] = {}
            for cavity in cavities:
                cavity_coordinates = Fpocket.load_fpocket_pdb_to_xyz(os.path.join(fpocket_cavities, cavity))
                ligand_coordinates = i[['x','y','z']].astype(float)
                results[count][cavity] = overlap.total_overlap(ligand_coordinates, cavity_coordinates)
            count += 1
        return Fpocket.best(results)

    @staticmethod
    def best(dict):
        """

        :param dict: nested dictionary where every key has a dictionary
        :return: the key in the nested dictionary with highest value
        """
        if len(dict)==0:
            return None
        temp = dict[1]
        for key1, d in dict.items():
            for key2, value in d.items():
                if value > temp[key2]:
                    temp[key2] = value
        winner = max(temp, key=lambda k: temp[k])
        if float(temp[winner]) > 0.5:
            return winner
        else:
            return None

    @staticmethod
    def check_line(line):
        if re.match(r'.*\d{2}-\d{2}.*', line):
            line = re.sub(r'(.*\d{2})(-\d{2}.*)', r'\1 \2', line)
        if re.match(r'.*\.\d{2}\d*\..*', line):
            line = re.sub(r'(.*\.\d{2})(\d*\..*)', r'\1 \2', line)
        if re.match(r'(ATOM|HETATM)\d+.*', line):
            line = re.sub(r'(ATOM|HETATM)(\d+.*)', r'\1 \2', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d[A-Z]{4}\s[A-Z].*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d)([A-Z]{4}\s[A-Z].*)', r'\1 \3', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d?\s+[A-Z]{3,4}\s[A-Z]{1}\d+\s+.*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d?\s+[A-Z]{3,4}\s[A-Z]{1})(\d+\s*.*)', r'\1 \3', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]+\d*\s+[A-Z]{3}\s+[A-Z]\d+.*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]+\d*\s+[A-Z]{3}\s+[A-Z])(\d+.*)', r'\1 \3', line)
        return line

    @staticmethod
    def get_ligands(pdb_file):
        """
        takes ligands in the pdb file
        :param pdb_file: the pdb file with enzyme and ligand
        :return: list with pandas dataframe with all ligands not separated (every element is part of other hetero/homo dimer)
        """
        results = []
        df = pd.DataFrame(columns=["atom", "AA", "chain", "number", "x", "y", "z"])
        with open(pdb_file, 'r') as file:
            lines = file.readlines()
        for line in lines:
            if line.startswith('HETATM') and not 'HOH' in line:
                line = Fpocket.check_line(line)
                df.loc[len(df)] = line.split()[2:9]
            elif line.startswith("ENDMDL"):
                results.append(df.copy())
                df = pd.DataFrame(columns=["atom", "AA", "chain", "number", "x", "y", "z"])
            else:
                continue
        results.append(df)
        split_results = []
        for i in results:
            temp = Fpocket.split_ligands(i)
            for i in temp:
                if len(i) > 5:
                    split_results.append(i)
        return split_results

    @staticmethod
    def load_fpocket_pdb_to_xyz(file_path):
        """

        :param file_path: fpocket pocket output file
        :return: pandas dataframe with xyz coordinates of cavity
        """
        df = pd.DataFrame(columns=['x', 'y', 'z'])
        with open(file_path, "r") as file:
            lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith('HETATM'):
                line = Fpocket.check_line(line)
                df.loc[len(df)] = list(filter(lambda x: len(x) > 0, re.split(r'[\n\t\s]+', line)))[6:9]
        df = df.apply(pd.to_numeric)
        return df

    @staticmethod
    def split_ligands(df):
        """

        :param df: df with all ligands
        :return: list with df's of separated ligands
        """
        ligands = []
        for i in df["chain"].unique():
            ligands.append(df[df['chain'] == i])
        return ligands

    @staticmethod
    def get_pockets(protein_folder, fpocket_folder):
        list_proteins = os.listdir(protein_folder)
        list_fpocket = os.listdir(fpocket_folder)
        results = {}
        for file in list_proteins:
            if os.path.isfile(os.path.join(protein_folder,file)):
                print(file)
                correspond = [x for x in list_fpocket if file.replace('.pdb', '_out') in x][0]
                protein_file = os.path.join(protein_folder, file)
                pocket_folder = os.path.join(fpocket_folder, correspond, "pockets")
                ligands = Fpocket.get_ligands(protein_file)
                pocket = Fpocket.find_cavity(ligands, pocket_folder)
                results[file] = pocket
        df = pd.DataFrame(list(results.items()), columns=["id", "pocket"])
        return df