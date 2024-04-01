from aaindex import aaindex1
from autocorrelation import get_normalized_props

def am_correlation(Ri, Rj):
    props = get_normalized_props(['PRAM900101', 'GRAR740102'])
    theta1 = props[0][Ri] * props[0][Rj]
    theta2 = props[1][Ri] * props[1][Rj]
    return theta1, theta2

def am_sequence_order_cor_factor(sequence, k=1):
    phoob = []
    philic = []
    for i in range(len(sequence)-k):
        theta1, theta2 = am_correlation(sequence[i], sequence[i+k])
        phoob.append(theta1)
        philic.append(theta2)
    cor_phoob = sum(phoob)/(len(sequence)-k)
    cor_philic = sum(philic)/(len(sequence)-k)
    return cor_phoob, cor_philic


def amp_pse_AAC(sequence, l=30, weight = 0.5):
    return None