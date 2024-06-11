from sequence.autocorrelation import get_normalized_props
from sequence.composition import aa_composition

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


def am_correlation(Ri, Rj):
    """

    :param Ri: ith amino acid
    :param Rj: jth amino acid
    :return: hydrophobicity and hydrophilicity correlation
    """
    props = get_normalized_props(['PRAM900101', 'GRAR740102'])
    theta1 = props['PRAM900101'][Ri] * props['PRAM900101'][Rj]
    theta2 = props['GRAR740102'][Ri] * props['GRAR740102'][Rj]
    return theta1, theta2


def am_sequence_order_cor_factor(sequence, k=1):
    """

    :param sequence: protein sequence
    :param k: window
    :return: sequence order correlation factor
    """
    phoob = 0
    philic = 0
    count = 0
    for i in range(len(sequence)-k):
        if not 'X' in (sequence[i], sequence[i+k]) :
            count +=1
            theta1, theta2 = am_correlation(sequence[i], sequence[i+k])
            phoob += theta1
            philic += theta2
    if count != 0:
        phoob = phoob/count
        philic = philic/count
    else:
        phoob = 0
        philic = 0
    return phoob, philic


def amp_pse_AAC(sequence, l=30, weight = 0.5):
    """

    :param sequence: protein sequence
    :param l: max window
    :param weight: weight of sequence order effect
    :return: amphiphilic pseudo amino acid count
    """
    aa_comp = aa_composition(sequence)
    sum_thetas = 0
    #in formula 2L, because phoob and philic, but in code at same time, so L.
    for i in range(l):
        sum_thetas += sum(am_sequence_order_cor_factor(sequence, i+1))
    devider = sum(aa_comp.values()) + weight*sum_thetas
    #first 20 components
    amp_pse_aac = {'ampPSE_'+aa: aa_comp['AAC_'+aa] / devider for aa in amino_acids}
    #other components
    numerator = [am_sequence_order_cor_factor(sequence, x + 1) for x in range(21, 21 + l)]
    for i in range(l):
        amp_pse_aac['amPse_'+ str(i)+".1"] = numerator[i][0] / devider
        amp_pse_aac['amPse_' + str(i)+".2"] = numerator[i][1] / devider
    return amp_pse_aac
