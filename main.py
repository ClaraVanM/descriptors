import pandas as pd
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
import distance_from_ligand

def check_aa(aa):
    if len(aa)>3:
        return aa[1:]
    else:
        return aa


def get_structures(protein_file, fpocket, pocket):
    """

    :param protein_file: pdb file of protein
    :param fpocket: fpocket output folder
    :param pocket: the corresponding output pocket for the binding cavity
    :return: cavity, ligand and name
    """
    cavity_selection = process_file.get_cavity_atoms(protein_file, fpocket, pocket)
    cavity, ligand = process_file.load_pdb(cavity_selection)
    protein, t = process_file.load_pdb(protein_file)
    # center data
    cavity[["x", "y", "z"]] = cavity[["x", "y", "z"]] - shape.COG(protein)
    protein[["x", "y", "z"]] = protein[["x", "y", "z"]] - shape.COG(protein)
    name = protein_file.split('/')[-1].split('.')[0]
    cavity['AA'] = cavity['AA'].apply(check_aa)
    return cavity, ligand, name, cavity_selection


def structure_descriptors(cavity, ligand):
    """

    :param cavity: cavity dataset
    :param ligand: ligand dataset
    :return: descriptors based on structure, depth, narrowness, AA comp per depth and exposed atoms per depth
    """
    axis = shape.find_cavity_axis(cavity, ligand)
    df = shape.residue_dist_from_axis(cavity, axis)
    cavity_pr = shape.projection(cavity, axis)
    df, depth = shape.add_buriedness(cavity, cavity_pr, axis)
    l_nar = shape.list_narrowness(df, axis, shape.COG(cavity))

    AA_comp = depth_comp.AA_per_buriedness(df)
    exposed = depth_comp.exposed_aa(df, l_nar)
    return depth, l_nar, AA_comp, exposed


def sequence_descriptors(cavity_name):
    sequence = process_file.get_sequence(cavity_name)
    aa_comp = composition.aa_composition(sequence)
    dipep_comp = composition.dipeptide_composition(sequence)
    tripep_comp = composition.tripeptide_composition(sequence)
    triad = composition.conjoint_triad(sequence)
    CTD_comp = CTD.ctd_composition(sequence)
    CTD_trans = CTD.ctd_transition(sequence)
    CTD_distr = CTD.ctd_distribution(sequence)
    pseaac = pseAAC.pseaac(sequence)
    ampseaac = amPseAAC.amp_pse_AAC(sequence)
    moreau_broto = autocorrelation.autocorrelation(sequence, autocorrelation.moreaubroto_ac)
    moran = autocorrelation.autocorrelation(sequence, autocorrelation.moran_ac)
    geary = autocorrelation.autocorrelation(sequence, autocorrelation.geary_ac)
    qsoc = sequence_order.tau_qsoc(sequence)
    return aa_comp, dipep_comp, tripep_comp, triad, CTD_comp, CTD_trans, CTD_distr, pseaac, ampseaac, moreau_broto, moran, geary, qsoc


def seq_struc_descriptors(cavity, ligand):
    sequence, cavity = distance_from_ligand.get_sequence(cavity, ligand)
    pseaac = pseAAC.pseaac(sequence)
    moran = autocorrelation.autocorrelation(sequence, autocorrelation.moran_ac)
    qsoc = sequence_order.tau_qsoc(sequence)
    cavity = distance_from_ligand.divide_cavity(cavity)
    AA_groups = distance_from_ligand.AA_per_buriedness(cavity)
    ctd_groups = distance_from_ligand.ctd_comp(cavity)
    return pseaac, moran, qsoc, AA_groups, ctd_groups


def get_results(protein_file, fpocket, pocket):
    """

    :param protein_file: pdb file of protein
    :param fpocket: fpocket output of protein
    :return: calculates all descriptors and puts them in dictionary
    """
    results = dict()
    cavity, ligand, name, cavity_name = get_structures(protein_file, fpocket, pocket)
    assert not ligand.empty
    results['name'] = name

    depth, l_nar, AA_comp, exposed = structure_descriptors(cavity.copy(), ligand)
    results['depth'] = depth
    for i in range(len(l_nar)):
        results['narrowness_' + str(i)] = l_nar[i]
    for values in AA_comp.values():
        results.update(values)
    results.update(exposed)

    aa_comp, dipep_comp, tripep_comp, triad, CTD_comp, CTD_trans, CTD_distr, pseaac, ampseaac, moreaubroto, moran, geary, qsoc = sequence_descriptors(cavity.copy())
    results.update(aa_comp)
    results.update(dipep_comp)
    results.update(tripep_comp)
    results.update(triad)
    for values in CTD_comp.values():
        results.update(values)
    for values in CTD_trans.values():
        results.update(values)
    for values in CTD_distr.values():
        results.update(values)
    results.update(pseaac)
    results.update(ampseaac)
    for values in moreaubroto.values():
        results.update(values)
    for values in moran.values():
        results.update(values)
    for values in geary.values():
        results.update(values)
    results.update(qsoc)

    pseaac_ss, moran_ss, qsoc_ss, AA_groups, ctd_groups = seq_struc_descriptors(cavity.copy(), ligand)
    results.update(pseaac_ss)
    results.update(qsoc)
    for values in moran_ss.values():
        results.update(values)
    for values in AA_groups.values():
        results.update(values)
    for values in ctd_groups.values():
        results.update(values)
    return results


def main():
    # proces pdb files
    b = True
    protein_folder = "/home/r0934354/Downloads/not_3.2.1/structures"
    fpocket_folder = "/home/r0934354/Downloads/not_3.2.1/fpocket"
    correspond = pd.read_csv('ids_with_pockets_not_3.2.1', index_col=0)
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
    df = main()
    df.to_csv('out.csv')
    """results = get_results("/home/r0934354/Downloads/EC3.2.1/structures/5D5A.pdb", "/home/r0934354/Downloads/EC3.2.1/fpocket/5D5A_out", "pocket22_atm.pdb")
    print(results)"""
