import pandas as pd

import process_file
import shape
import depth_comp
import composition
import CTD
import pseAAC
import amPseAAC
import autocorrelation
import sequence_order


if __name__ == "__main__":
    #proces pdb files
    process_file.get_cavity_atoms("1B3X.pdb", "1B3X_out")
    cavity = process_file.load_pdb("1B3X_neighbor.pdb")
    protein = process_file.load_pdb("1GNY.pdb")

    #center data
    cavity[["x","y","z"]] = cavity[["x","y","z"]] - shape.COG(protein)
    protein[["x","y","z"]] = protein[["x","y","z"]] - shape.COG(protein)

    #get axis and do projection of cavity on this axis, used later for buriedness
    axis = shape.find_cavity_axis(cavity)
    cavity_pr = shape.projection(cavity, axis)

    #first descriptor depth, df is cavity df with residues devided into buriedness groups
    ################################################################################################################## depth = float
    df, depth = shape.add_buriedness(cavity,cavity_pr, axis)

    #second descriptor, narrowness for every buriedness level
    ################################################################################################################### l_nar is list with floats
    l_nar = shape.list_narrowness(df, axis, shape.COG(cavity))

    #add distances from redidues to axis to df
    df = shape.residue_dist_from_axis(df, axis)

    #3th and 4th descriptors, AA frequency per buriedness level and exposed residues
    ################################################################################################################# 3 is nested dictionary, 4 is dictionary
    descriptor3 = depth_comp.AA_per_buriedness(df)
    descriptor4 = depth_comp.exposed_aa(df, l_nar)

    #get sequence with x's where gap in cavity sequence
    sequence = process_file.get_sequence("1B3X_neighbor.pdb")

    #5th, 6th, 7th and 8th descriptor
    ################################################################################################################# 5, 6, 7, 8 dictionaries
    AAC = composition.aa_composition(sequence)
    di_AAC = composition.dipeptide_composition(sequence)
    tri_AAC = composition.tripeptide_composition(sequence)
    triads_AAC = composition.conjoint_triad(sequence)

    #9th, 10th, 11th descriptor
    ##################################################################################################################### 9,10,11 nested dictionaries
    CTD_comp = CTD.ctd_composition(sequence)
    CTD_trans = CTD.ctd_transition(sequence)
    CTD_distr = CTD.ctd_distribution(sequence)

    #12th descriptor
    ################################################################################################# dictionary
    pse_AAC = pseAAC.pseaac(sequence)

    #13th descriptor
    ################################################################################################ dictionary
    am_pse_AAC = amPseAAC.amp_pse_AAC(sequence)

    #14th, 15th and 16th descriptor
    ##################################################################################################### nested dictionary
    moreau_broto = autocorrelation.autocorrelation(sequence, autocorrelation.moreaubroto_ac)
    moran = autocorrelation.autocorrelation(sequence, autocorrelation.moran_ac)
    geary = autocorrelation.autocorrelation(sequence, autocorrelation.geary_ac)

    #17th descriptor
    ##################################################################################################### dictionary
    qsoc = sequence_order.tau_qsoc(sequence)

    #create df
    df = pd.DataFrame(columns=['name', "depth"])
    df['name'] = "1B3X"
    df["depth"] = depth
    print(df)
