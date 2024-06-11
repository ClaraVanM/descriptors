import os
import pandas as pd
import re
import numpy as np


class Structures:
    AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    def __init__(self, pdbfile):
        self.name = pdbfile.split('/')[-1].split('.')[0]
        self.structures = self.getStructure(pdbfile)
        self.protein = self.extract_protein()


    @staticmethod
    def getStructure(pdbfile):
        """

        :param pdbfile: pdb structure file
        :return: dataframe with structure coordinates, amino acids and their order, atom and type of structure
        """
        assert os.path.isfile(pdbfile), "No valid file."
        df = pd.DataFrame(columns=['type', 'atom', 'AA', 'AA_number', 'x', 'y', 'z'])
        with open(pdbfile, 'r') as f:
            line = f.readline()
            while line:
                if line.startswith('ATOM') or line.startswith('HETATM') and not "HOH" in line and not 'CA' in line and not "FE" in line and not 'SO4' in line and not 'ZN' in line and not "5FC" in line:
                    line = Structures.check_pdbline(line)
                    df.loc[len(df)] = [line.split()[i] for i in [0,2, 3, 5, 6, 7, 8]]
                line = f.readline()
        df[['x', 'y', 'z']] = df[['x', 'y', 'z']].astype(float)
        df[["x", "y", "z"]] = Structures.outliers(df[["x", "y", "z"]])
        df.dropna(inplace=True)
        df['AA'] = df['AA'].apply(Structures.check_aa)
        df = Structures.check_atoms(df)
        return df

    @staticmethod
    def check_pdbline(line):
        """
        checks structure of the line and corrects if necessary
        :param line: a line of a pdb file
        :return:corrected line
        """
        if re.match(r'.*\d{2}-\d{2}.*', line):
            line = re.sub(r'(.*\d{2})(-\d{2}.*)', r'\1 \2', line)
        if re.match(r'.*\.\d{2}\d*\..*', line):
            line = re.sub(r'(.*\.\d{2})(\d*\..*)', r'\1 \2', line)
        if re.match(r'HETATM\d+.*', line):
            line = re.sub(r'^(HETATM)(\d+.*)', r'\1 \2', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d[A-Z]{4}\s[A-Z].*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d)([A-Z]{4}\s[A-Z].*)', r'\1 \3', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d?\s+[A-Z]{3,4}\s[A-Z]{1}\d+\s+.*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]{2}\d?\s+[A-Z]{3,4}\s[A-Z]{1})(\d+\s*.*)', r'\1 \3', line)
        if re.match(r'^(ATOM|HETATM)\s+\d+\s+[A-Z]+\d*\s+[A-Z]{3}\s+[A-Z]\d+.*', line):
            line = re.sub(r'(^(ATOM|HETATM)\s+\d+\s+[A-Z]+\d*\s+[A-Z]{3}\s+[A-Z])(\d+.*)', r'\1 \3', line)
        return line

    @staticmethod
    def check_atoms(df):
        """
        checks if line that is noted to be part of the structure is not a ligand
        :param df: structure df
        :return: corrected df
        """
        for index, row in df.iterrows():
            if row['type'] == 'ATOM':
                if row['AA'] not in Structures.AA_symb:
                    df.loc[index,'type'] = 'HETATM'
        return df

    @staticmethod
    def outliers(cavity):
        """
        detects if some points are far away and remove them
        :param cavity: cavity coordinates
        :return: cavity coordinates without outliers
        """
        q1 = np.percentile(cavity, 25, axis=0)
        q3 = np.percentile(cavity, 75, axis=0)
        iqr = q3 - q1
        lower_bound = q1 - 2 * iqr
        upper_bound = q3 + 2 * iqr
        outliers = np.any((cavity < lower_bound) | (cavity > upper_bound), axis=1)
        cleaned_cavity = cavity[~outliers]
        return cleaned_cavity

    def extract_protein(self):
        """

        :return: protein without ligand
        """
        return self.structures[self.structures['type']=='ATOM']

    @staticmethod
    def check_aa(aa):
        """
        checks if amino acid consists out of 3 characters and removes chain identifier if not
        :param aa: amino acid
        :return: corrected amino acid
        """
        if len(aa) > 3:
            return aa[1:]
        else:
            return aa