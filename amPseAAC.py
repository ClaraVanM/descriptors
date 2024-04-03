from autocorrelation import get_normalized_props
from composition import aa_composition

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


def am_correlation(Ri, Rj):
    props = get_normalized_props(['PRAM900101', 'GRAR740102'])
    theta1 = props[0][Ri] * props[0][Rj]
    theta2 = props[1][Ri] * props[1][Rj]
    return theta1, theta2


def am_sequence_order_cor_factor(sequence, k=1):
    phoob = 0
    philic = 0
    for i in range(len(sequence)-k):
        theta1, theta2 = am_correlation(sequence[i], sequence[i+k])
        phoob += theta1
        philic += theta2
    phoob = phoob/(len(sequence)-k)
    philic = philic/(len(sequence)-k)
    return phoob, philic


def amp_pse_AAC(sequence, l=30, weight = 0.5):
    sequence = sequence.upper()
    aa_comp = aa_composition(sequence)
    sum_thetas = 0
    for i in range(2*l):
        sum_thetas += sum(am_sequence_order_cor_factor(sequence, i+1))
    devider = sum(aa_comp.values()) + weight*sum_thetas
    #first 20 components
    amp_pse_aac = {aa: aa_comp[aa] / devider for aa in amino_acids}
    #other components
    numerator = [am_sequence_order_cor_factor(sequence, x + 1) for x in range(20, 2*l)]
    for i in range(20, 2*l):
        amp_pse_aac[i+1] = (weight * numerator[i]) / devider
    return amp_pse_aac
