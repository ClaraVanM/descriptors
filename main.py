import pandas as pd
import sys
import os

import process_file
import shape
import depth_comp
import composition
import CTD
import pseAAC
import amPseAAC
import autocorrelation
import sequence_order


def get_results(protein_file, fpocket):
    """

    :param protein_file: pdb file of protein
    :param fpocket: fpocket output of protein
    :return: calculates all descriptors and puts them in dictionary
    """
    cavity_selection = process_file.get_cavity_atoms(protein_file, fpocket)
    cavity = process_file.load_pdb(cavity_selection)
    protein = process_file.load_pdb(protein_file)

    # center data
    cavity[["x", "y", "z"]] = cavity[["x", "y", "z"]] - shape.COG(protein)
    protein[["x", "y", "z"]] = protein[["x", "y", "z"]] - shape.COG(protein)

    results = dict()
    results['name'] = protein_file.split('/')[-1].split('.')[0]
    # get axis and do projection of cavity on this axis, used later for buriedness
    axis = shape.find_cavity_axis(cavity)
    cavity_pr = shape.projection(cavity, axis)

    # first descriptor depth, df is cavity df with residues devided into buriedness groups
    ################################################################################################################## depth = float
    df, depth = shape.add_buriedness(cavity, cavity_pr, axis)
    results['depth'] = depth

    # second descriptor, narrowness for every buriedness level
    ################################################################################################################### l_nar is list with floats
    l_nar = shape.list_narrowness(df, axis, shape.COG(cavity))
    for i in range(len(l_nar) - 1):
        results['narrowness_' + str(i)] = l_nar[i]

    # add distances from redidues to axis to df
    df = shape.residue_dist_from_axis(df, axis)

    # 3th and 4th descriptors, AA frequency per buriedness level and exposed residues
    ################################################################################################################# 3 is nested dictionary, 4 is dictionary
    descriptor3 = depth_comp.AA_per_buriedness(df)
    for values in descriptor3.values():
        results.update(values)
    results.update(depth_comp.exposed_aa(df, l_nar))

    # get sequence with x's where gap in cavity sequence
    sequence = process_file.get_sequence(cavity_selection)

    # 5th, 6th, 7th and 8th descriptor
    ################################################################################################################# 5, 6, 7, 8 dictionaries
    results.update(composition.aa_composition(sequence))
    results.update(composition.dipeptide_composition(sequence))
    results.update(composition.tripeptide_composition(sequence))
    results.update(composition.conjoint_triad(sequence))

    # 9th, 10th, 11th descriptor
    ##################################################################################################################### 9,10,11 nested dictionaries
    CTD_comp = CTD.ctd_composition(sequence)
    CTD_trans = CTD.ctd_transition(sequence)
    CTD_distr = CTD.ctd_distribution(sequence)
    for values in CTD_comp.values():
        results.update(values)
    for values in CTD_trans.values():
        results.update(values)
    for values in CTD_distr.values():
        results.update(values)

    # 12th descriptor
    ################################################################################################# dictionary
    results.update(pseAAC.pseaac(sequence))

    # 13th descriptor
    ################################################################################################ dictionary
    results.update(amPseAAC.amp_pse_AAC(sequence))

    # 14th, 15th and 16th descriptor
    ##################################################################################################### nested dictionary
    moreau_broto = autocorrelation.autocorrelation(sequence, autocorrelation.moreaubroto_ac)
    for values in moreau_broto.values():
        results.update(values)
    moran = autocorrelation.autocorrelation(sequence, autocorrelation.moran_ac)
    for values in moran.values():
        results.update(values)
    geary = autocorrelation.autocorrelation(sequence, autocorrelation.geary_ac)
    for values in geary.values():
        results.update(values)

    # 17th descriptor
    ##################################################################################################### dictionary
    results.update(sequence_order.tau_qsoc(sequence))
    return results

def main():
    # proces pdb files
    b = True
    protein_folder = sys.argv[1]   # folder with all protein .pdb files
    fpocket_folder = sys.argv[2]        # folder with fpocket output for every protein
    fpocket_list = os.listdir(fpocket_folder)
    for file in os.listdir(protein_folder):
        if os.path.isfile(os.path.join(protein_folder,file)):
            name = file.replace('.pdb','_out')
            fpocket = [x for x in fpocket_list if name in x][0]
            if b:
                df = pd.DataFrame(get_results(os.path.join(protein_folder, file), os.path.join(fpocket_folder, fpocket)), index=[0])
                b = False
            else:
                df.loc[len(df)] = get_results(os.path.join(protein_folder, file), os.path.join(fpocket_folder, fpocket))
    return df


if __name__ == "__main__":
    df = main()
    df.to_csv('out.csv')



#ipv AA per buriedness, AA verdelen in eigenschappen zoals hydro, arom,... === descriptor 3
#descriptor 4 is C, N, O genoeg of OE, NE,...