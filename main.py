import pandas as pd
import os

from sequence import Sequence
from Getdata.Structures import Structures
from Getdata.Cavity import Cavity
from shape.Shape import Shape
from FindCavity.fpocket import Fpocket
from Distance.Distance import Distance


def get_results(protein_file, fpocket, pocket):
    """

    :param protein_file: pdb file of protein
    :param fpocket: fpocket output of protein
    :return: calculates all descriptors and puts them in dictionary
    """
    results = dict()
    #import structures
    protein = Structures(protein_file)
    cavity = Cavity(protein_file, fpocket, pocket)
    assert not cavity.ligand.empty
    results['name'] = protein.name
    #get shape descriptors
    shape = Shape(cavity.cavity, cavity.ligand)
    results.update(shape.getDescriptors())
    #get sequence descriptors
    seq = Sequence.Sequence(cavity.cavity)
    results.update(seq.sequence_descriptors())
    # get distance descriptors
    d = Distance(cavity.cavity, cavity.ligand)
    results.update(d.getDescriptors())
    return results


def main():
    # proces pdb files
    b = True
    protein_folder = "/home/r0934354/Downloads/not_3.2.1/structures"
    fpocket_folder = "/home/r0934354/Downloads/not_3.2.1/fpocket"
    correspond = pd.read_csv('ids_with_pockets_not_3.2.1.csv', index_col=0)
    fpocket_list = os.listdir(fpocket_folder)
    for file in os.listdir(protein_folder):
        print(file)
        if pd.isna(list(correspond.loc[correspond['id'] == file, ['pocket']]['pocket'])[0]):
            continue

        else:
            try:
                pocket = list(correspond.loc[correspond['id'] == file, ['pocket']]['pocket'])[0]
                print(pocket)
                if os.path.isfile(os.path.join(protein_folder,file)):
                    name = file.replace('.pdb','_out')
                    fpocket = [x for x in fpocket_list if name in x][0]
                    if b:
                        df = pd.DataFrame(get_results(os.path.join(protein_folder, file), os.path.join(fpocket_folder, fpocket), pocket), index=[0])
                        b = False
                    else:
                        df.loc[len(df)] = get_results(os.path.join(protein_folder, file), os.path.join(fpocket_folder, fpocket), pocket)
            except:
                df.loc[len(df)] = [file] + [-2] * (len(df.columns)-1)
    return df




if __name__ == "__main__":
    """pockets = Fpocket("/home/r0934354/Downloads/not_3.2.1/structures", "/home/r0934354/Downloads/not_3.2.1/fpocket")
    pockets.pockets.to_csv("ids_with_pockets_not_3.2.1.csv")"""
    '''df = main()
    df.to_csv('out.csv')'''
    results = get_results("C:/Users/32496/Desktop/not_3.2.1/structures\\8PPZ.pdb", "C:/Users/32496/Desktop/not_3.2.1/fpocket\\8PPZ_out", "pocket1_atm.pdb")
