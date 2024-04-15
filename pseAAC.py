from autocorrelation import get_normalized_props
from composition import aa_composition

# list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]


def pseaac(sequence, l=30, weight=0.05, properties=['PRAM900101', 'GRAR740102', 'PONJ960101']):
    """

    :param sequence: protein sequence
    :param l: max window, integer < length of sequence
    :param weight: weight of sequence order effect
    :param properties: identifiers properties from aaindex
    :return: pseudo amino acid count
    """
    props = get_normalized_props(properties)
    sum_thetas = 0
    for i in range(l):
        sum_thetas += sequence_order_correlation_factor(sequence, props, i+1)
    aa_comp = aa_composition(sequence)
    devider = sum(aa_comp.values()) + weight*sum_thetas
    # first 20 components
    pseudo_aac = {aa:aa_comp[aa]/devider for aa in amino_acids}
    #other components
    numerator = [sequence_order_correlation_factor(sequence, props, x+1) for x in range (l-20)]
    for i in range(l-20):
        pseudo_aac[i] = (weight * numerator[i])/devider
    return pseudo_aac


def correlation(Ri, Rj, props):
    """

    :param Ri: ith amino acid
    :param Rj: jth amino acid
    :param props: normalized property values from aaindex
    :return: amino acid correlation factor
    """
    theta = 0
    for i in props.values():
        theta += (i[Ri] - i[Rj])**2
    return theta/len(props)


# k determines the order of the factor
def sequence_order_correlation_factor(sequence, props, k=1):
    """

    :param sequence: protein sequence
    :param props: normalized property values from aaindex
    :param k: window
    :return: sequence order-correlated factor theta
    """
    cor = float(0)
    count = 0
    for i in range(len(sequence)-k):
        if not 'X' in (sequence[i], sequence[i+k]):
            count += 1
            cor += correlation(sequence[i], sequence[i+k], props)
    if count == 0:
        cor = 0
    else:
        cor = cor/count
    return cor
