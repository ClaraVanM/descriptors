import pymolPy3
import os
import pandas as pd
import re
import numpy as np

AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def get_cavity_atoms(protein_file, fpocket_out, pocket_file):
    """
    Selects all atoms/AA around the predicted cavity and puts them in new file
    :param protein_file: full protein.pdb file
    :param fpocket_out: output folder of fpocket after cavity prediction
    :return: None
    """
    pm = pymolPy3.pymolPy3(0)
    pm(f"load {protein_file}")
    pm(f"load {fpocket_out}/pockets/{pocket_file}")
    pm(f"select cavity, {pocket_file.split('.')[0]}")
    protein_name = protein_file.split("/")[-1].split(".")[0]
    pm(f"select prot, {protein_name}")
    pm(f"select neighborhood, prot near_to 10 of cavity")
    pm(f"save {protein_name}_neighbor.pdb, neighborhood")
    return protein_name+'_neighbor.pdb'


def load_pdb(file):
    """
    load a pdb file and output a pandas dataframe with x, y, z coordinates and additional information from pdb file
    :param file: path to file
    :return: pandas dataframe with coordinates of cavity and ligand
    """
    assert os.path.isfile(file), "No valid file."
    df_cavity = pd.DataFrame(columns=['atom', 'AA','AA_number', 'x', 'y', 'z'])
    df_ligand = pd.DataFrame(columns=['atom', 'AA','AA_number', 'x','y','z'])
    with open(file, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('ATOM'):
                # sometimes indents in pdb files are wronge so next 3 if statements correct them if they occure
                if re.match(r'.*\d{2}-\d{2}.*', line):
                    line = re.sub(r'(.*\d{2})(-\d{2}.*)', r'\1 \2', line)
                if re.match(r'.*\.\d{2}\d*\..*', line):
                    line = re.sub(r'(.*\.\d{2})(\d*\..*)', r'\1 \2', line)
                if re.match(r'^(ATOM|HETATM)\s*\d*\s*[A-Z]{2}\d[A-Z]{4}\s[A-Z].*', line):
                    line = re.sub(r'(^(ATOM|HETATM)\s*\d*\s*[A-Z]{2}\d)([A-Z]{4}\s[A-Z].*)', r'\1 \3', line)
                if re.match(r'^(ATOM|HETATM)\s*\d*\s*[A-Z]{2}\d?\s*[A-Z]{3,4}\s[A-Z]{1}\d*\s*.*', line):
                    line = re.sub(r'(^(ATOM|HETATM)\s*\d*\s*[A-Z]{2}\d?\s*[A-Z]{3,4}\s[A-Z]{1})(\d*\s*.*)', r'\1 \3', line)
                df_cavity.loc[len(df_cavity)] = [line.split()[i] for i in [2,3,5,6,7,8]]
            if line.startswith("HETATM") and not "HOH" in line and not 'CA' in line and not "FE" in line and not 'SO4'in line:
                if re.match(r'.*\d{2}-\d{2}.*', line):
                    line = re.sub(r'(.*\d{2})(-\d{2}.*)', r'\1 \2', line)
                if re.match(r'.*\.\d{2}\d*\..*', line):
                    line = re.sub(r'(.*\.\d{2})(\d*\..*)', r'\1 \2', line)
                if re.match(r'HETATM\d*.*', line):
                    line = re.sub(r'(HETATM)(\d*.*)', r'\1 \2', line)
                df_ligand.loc[len(df_ligand)] = [line.split()[i] for i in [2,3,5,6,7,8]]
            line = f.readline()
    df_cavity[['x', 'y', 'z']] = df_cavity[['x', 'y', 'z']].astype(float)
    df_ligand[['x', 'y', 'z']] = df_ligand[['x', 'y', 'z']].astype(float)
    df_cavity[["x", "y", "z"]] = outliers(df_cavity[["x", "y", "z"]])
    df_cavity.dropna()
    #sometimes ligand under atom instead of hetatm fix here
    ligand_rows = df_cavity[~df_cavity['AA'].isin(AA_symb.keys())]
    df_ligand = pd.concat([ligand_rows, df_ligand])
    df_cavity = df_cavity[df_cavity['AA'].isin(AA_symb.keys())]
    return df_cavity, df_ligand


def get_sequence(cavity):
    """

    :param pdbfile: pdb file of structure
    :return: string with sequence of structure
    """
    cavity['AA_number'] = cavity['AA_number'].astype(float)
    cavity = cavity.sort_values(by=['AA_number'], ascending=True)
    cavity = cavity.drop_duplicates(subset=["AA_number"])
    sequence = ''.join(cavity["AA"].map(AA_symb))
    return sequence


def outliers(cavity):
    """
    detects if some points are far away and remove them
    :param cavity: cavity coordinates
    :return: cavity coordinates without outliers
    """
    q1 = np.percentile(cavity, 25, axis=0)
    q3 = np.percentile(cavity, 75, axis=0)
    iqr = q3-q1
    lower_bound = q1 - 5 * iqr
    upper_bound = q3 +5*iqr
    outliers = np.any((cavity < lower_bound) | (cavity > upper_bound), axis=1)
    cleaned_cavity = cavity[~outliers]
    return cleaned_cavity