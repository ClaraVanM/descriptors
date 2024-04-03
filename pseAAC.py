from autocorrelation import get_normalized_props
from composition import aa_composition

# list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


def pseaac(sequence, l=30, weight=0.05, properties=['PRAM900101', 'GRAR740102', 'PONJ960101']):
    #l determines the number of thetas, every theta increases the calculating window
    sequence = sequence.upper()
    if l < 0 or l > len(sequence):
        l = 30
    props = get_normalized_props(properties)
    sum_thetas = 0
    for i in range(l):
        sum_thetas += sequence_order_correlation_factor(sequence, i+1, props)
    aa_comp = aa_composition(sequence)
    devider = sum(aa_comp.values()) + weight*sum_thetas
    # first 20 components
    pseudo_aac = {aa:aa_comp[aa]/devider for aa in amino_acids}
    #other components
    numerator = [sequence_order_correlation_factor(sequence, x+1, props) for x in range (l-20)]
    for i in range(l-20):
        pseudo_aac[i] = (weight * numerator[i])/devider
    return pseudo_aac


def correlation(Ri, Rj, props):
    theta = 0
    for i in props:
        theta += (i[Ri] - i[Rj])**2
    return theta/len(props)


# k determines the order of the factor
def sequence_order_correlation_factor(sequence, props, k=1):
    cor = []
    for i in range(len(sequence)-k):
        cor.append(correlation(sequence[i], sequence[i+k], props))
    return cor