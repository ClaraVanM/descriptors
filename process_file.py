import pymolPy3
import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio import SeqIO


def get_cavity_atoms(protein_file, fpocket_out):
    """
    Selects all atoms/AA around the predicted cavity and puts them in new file
    :param protein_file: full protein.pdb file
    :param fpocket_out: output folder of fpocket after cavity prediction
    :return: None
    """
    pm = pymolPy3.pymolPy3(0)
    pm(f"load {protein_file}")
    pm(f"load {fpocket_out}/pockets/pocket1_atm.pdb")
    pm(f"select cavity, pocket1_atm")
    protein_name = protein_file.split("/")[-1].split(".")[0]
    pm(f"select prot, {protein_name}")
    pm(f"select neighborhood, prot near_to 10 of cavity")
    print(protein_name+'_neighbor.pdb')
    pm(f"save {protein_name}_neighbor.pdb, neighborhood")
    return protein_name+'_neighbor.pdb'


def load_pdb(file):
    """
    load a pdb file and output a pandas dataframe with x, y, z coordinates and additional information from pdb file
    :param file: path to file
    :return: pandas dataframe with coordinates
    """
    print(file)
    assert os.path.isfile(file), "No valid file."
    df = pd.DataFrame(columns=['atom', 'AA','AA_number', 'x', 'y', 'z'])
    with open(file, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('ATOM'):
                df.loc[len(df)] = [line.split()[i] for i in [2,3,5,6,7,8]]
            line = f.readline()
    df[['x', 'y', 'z']] = df[['x', 'y', 'z']].astype(float)
    return df


def get_sequence(pdbfile):
    record = list(SeqIO.parse(pdbfile, "pdb-atom"))
    return record[0].seq