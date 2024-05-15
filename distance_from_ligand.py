import numpy as np
from scipy.spatial.distance import euclidean
import pandas as pd
from scipy.spatial.distance import cdist
properties = {"hydrophobicity":{"ARG":1, "LYS":1, "GLU":1, "ASP":1, "GLN":1, "ASN":1,"GLY":2,"ALA":2,"SER":2,"THR":2,"PRO":2,"HIS":2,"TYR":2,"CYS":3,"LEU":3,"VAL":3,"ILE":3,"MET":3,"PHE":3,"TRP":3},
              "normalized_vdw":{"GLY":1,"CYS":1,"ALA":1,"SER":1,"THR":1,"PRO":1,"ASP":1,"ASN":2,"VAL":2,"GLU":2,"GLN":2,"ILE":2,"LEU":2,"MET":3,"HIS":3,"LYS":3,"PHE":3,"ARG":3,"TYR":3,"TRP":3},
              "polarity":{"LEU":1,"ILE":1,"PHE":1,"TRP":1,"CYS":1,"MET":1,"VAL":1,"TYR":1,"ALA":2,"PRO":2,"THR":2,"GLY":2,"SER":2,"LYS":3,"ASN":3,"HIS":3,"ARG":3,"GLN":3,"GLU":3,"ASP":3},
              "charge":{"LYS":1,"ARG":1,"ALA":2,"ASN":2,"CYS":2,"GLN":2,"GLY":2,"HIS":2,"ILE":2,"LEU":2,"MET":2,"PHE":2,"PRO":2,"SER":2,"THR":2,"TRP":2,"TYR":2,"VAL":2,"ASP":3,"GLU":3},
              "secondary_struct":{"GLU":1,"ALA":1,"LEU":1,"MET":1,"GLN":1,"LYS":1,"ARG":1,"HIS":1,"VAL":2,"ILE":2,"TYR":2,"CYS":2,"TRP":2,"PHE":2,"THR":2,"GLY":3,"ASN":3,"PRO":3,"SER":3,"ASP":3},
              "solvent_accessibility":{"ALA":1,"LEU":1,"PHE":1,"CYS":1,"GLY":1,"ILE":1,"VAL":1,"TRP":1,"ARG":2,"LYS":2,"GLN":2,"GLU":2,"ASN":2,"ASP":2,"MET":3,"PRO":3,"SER":3,"THR":3,"HIS":3,"TYR":3},
              "polarizability":{"GLY":1,"ALA":1,"SER":1,"ASP":1,"THR":1,"CYS":2,"PRO":2,"ASN":2,"VAL":2,"GLU":2,"GLN":2,"ILE":2,"LEU":2,"LYS":3,"MET":3,"HIS":3,"PHE":3,"ARG":3,"TYR":3,"TRP":3}}
AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

##problem: if an amino acid is not at that distance, an NaN value will be returned. this should be 0

def ctd_comp(cavity):
    cavity =  cavity.copy()
    ctd = {}
    for prop in properties.keys():
        cavity.loc[:,"prop"] = [properties[prop][aa] for aa in cavity['AA']]
        for i in cavity["group"].unique():
            df_temp = cavity[cavity["group"] == i]
            prop_count = df_temp['prop'].value_counts().to_dict()
            ctd[prop + str(i)] = prop_count
    return ctd


def COG(data):
    """

    :param data: processed pdb file
    :return: center of the xyz coordinates
    """
    return data[["x","y","z"]].mean()


def dist_from_ligand(cavity, COG_ligand):
    """

    :param cavity: cavity dataframe
    :param COG_ligand: center of ligand
    :return: cavity with added distance from ligand for each point
    """
    cavity['dist_lig'] = [euclidean(point, COG_ligand.values) for point in cavity[['x','y','z']].values]
    return cavity


def sequence(cavity):
    """

    :param cavity: cavity dataframe with distance from ligand
    :return: sequence of cavity based on distance from ligand
    """
    sorted_cavity = cavity.sort_values(by='dist_lig')
    sequence = ''.join(sorted_cavity['AA'].map(AA_symb))
    return sequence


def get_sequence(cavity, ligand):
    COG_ligand = COG(ligand)
    # drop dublicate AA number by taking center of gravity of every AA
    cog = cavity.groupby('AA_number')[['x', 'y', 'z']].mean()
    cavity = pd.merge(cavity, cog, on='AA_number', suffixes=('', 'center'))
    cavity = cavity.drop(['atom', 'x','y','z', 'index'], axis=1)
    cavity = cavity.rename(columns={'xcenter':"x", 'ycenter':"y", 'zcenter':'z'})
    cavity = dist_from_ligand(cavity, COG_ligand)
    seq = sequence(cavity)
    return seq, cavity


def divide_cavity(cavity):
    deepness = cavity['dist_lig'].max() - cavity['dist_lig'].min()
    # make 5 intervals
    jumps = deepness / 5
    cavity['group'] = -1
    for i in range(5):
        selection_index = cavity[
            (cavity["dist_lig"] >= cavity['dist_lig'].min() + jumps * (i)) & (
                    cavity["dist_lig"] <= cavity['dist_lig'].min() + jumps * (i + 1))].index
        cavity.loc[selection_index, 'group'] = i
    return cavity


def AA_per_buriedness(df):
    """

    :param df: df with buriedness level column and distances from axis.
    :return: dictionary with frequency of amino acids per buriedness level
    """
    total_count = {}
    for i in df["group"].unique():
        df_temp = df[df["group"]==i]
        aa_count = df_temp['AA'].value_counts()
        total_count[i] = {'group_' + aa+str(i):aa_count.get(aa,0) for aa in AA_symb.keys()}
    return total_count


