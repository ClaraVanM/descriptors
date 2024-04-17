import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from shape import COG
from autocorrelation import get_normalized_props
from skspatial.objects import Sphere

AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# do the autocorrelation, but now per buriedness, to have more effect of
#structure in the sequence measure

def get_sequence_per_buriedness(dataframe):
    """

    :param dataframe: pandas dataframe with coordinates, AA, AA_numbers and buriedness
    :return: list with 5 elements, each element is a string of the sequence of that buriedness level
    """
    sequences = []
    df = dataframe.drop_duplicates(subset=['AA_number'])
    for i in range(5):
        subset = df[df['buriedness']==i]
        sequence = subset["AA"].str.cat()
        sequences.append(sequence)
    return sequences

#first idea is to put all coordinates on a circle and then calculate autocorrelation and other distributions
#but then starting piont needs to be aligned for all structures in order to be comparable
#sequence order is already aligned. so we work with this even tho they are then randomly oriented in space.
df = pd.read_csv('dataset')

#give every coordinate of AA a value from AA index and project these values on a cone/ ball
# then see per buriedness if there is a concentrated spot for that property?
def add_properties(dataframe, properties = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]):
    """

    :param dataframe: pandas dataframe with atoms, AA, coordinates
    :param properties: which properties do we want to add
    :return: dataframe with AA and their corresponding values of the properties
    """
    props = get_normalized_props(properties)
    # choose the closest point of aa to axis as representative for aa
    df_index= dataframe.groupby('AA')['dist_from_axis'].idxmin()
    df = dataframe.loc[df_index]
    df['AA_symb'] = df['AA'].map(AA_symb)
    for i, values in props.items():
        print(i, values)
        df[i] = df['AA_symb'].map(values)
    return df



