from itertools import product
import math
# <editor-fold desc="AA grouped per property">
properties = {"hydrophobicity":{"R":1, "K":1, "E":1, "D":1, "Q":1, "N":1,"G":2,"A":2,"S":2,"T":2,"P":2,"H":2,"Y":2,"C":3,"L":3,"V":3,"I":3,"M":3,"F":3,"W":3},
              "normalized_vdw":{"G":1,"A":1,"S":1,"T":1,"P":1,"D":1,"N":2,"V":2,"E":2,"Q":2,"I":2,"L":2,"M":3,"H":3,"K":3,"F":3,"R":3,"Y":3,"W":3},
              "polarity":{"L":1,"I":1,"F":1,"W":1,"C":1,"M":1,"V":1,"Y":1,"A":2,"P":2,"T":2,"G":2,"S":2,"K":3,"N":3,"H":3,"R":3,"Q":3,"E":3,"D":3},
              "charge":{"K":1,"R":1,"A":2,"N":2,"C":2,"Q":2,"G":2,"H":2,"I":2,"L":2,"M":2,"F":2,"P":2,"S":2,"T":2,"W":2,"Y":2,"V":2,"D":3,"E":3},
              "secondary_struct":{"E":1,"A":1,"L":1,"M":1,"Q":1,"K":1,"R":1,"H":1,"V":2,"I":2,"Y":2,"C":2,"W":2,"F":2,"T":2,"G":3,"N":3,"P":3,"S":3,"D":3},
              "solvent_accessibility":{"A":1,"L":1,"F":1,"C":1,"G":1,"I":1,"V":1,"W":1,"R":2,"K":2,"Q":2,"E":2,"N":2,"D":2,"M":3,"P":3,"S":3,"T":3,"H":3,"Y":3},
              "polarizability":{"G":1,"A":1,"S":1,"D":1,"T":1,"C":2,"P":2,"N":2,"V":2,"E":2,"Q":2,"I":2,"L":2,"K":3,"M":3,"H":3,"F":3,"R":3,"Y":3,"W":3}}
# hydrophobicity: '1' -> Polar; '2' -> Neutral, '3'
# normalized_vdw: '1' -> (0-2.78); '2' -> (2.95-4.0), '3' -> (4.03-8.08)
# polarity:'1' -> (4.9-6.2); '2' -> (8.0-9.2), '3' -> (10.4-13.0)
# charge: '1' -> Positive; '2' -> Neutral, '3' -> Negative
# secondary structure: '1' -> Helix; '2' -> Strand, '3' -> coil
# solvent_accessibility'1' -> Buried; '2' -> Exposed, '3' -> Intermediate
# polarizability: '1' -> (0-0.108); '2' -> (0.128-0.186), '3' -> (0.219-0.409)
# </editor-fold>


def str_to_num(sequence, prop):
    sequence = sequence.upper()
    seq=""
    for i in sequence:
        seq += str(properties[prop][i])
    return seq


def ctd_composition(sequence):
    sequence = sequence.upper()
    comp = {}
    count = {}
    for prop, values in properties.items():
        converted_sequence = str_to_num(sequence, prop)
        for i in range(1,4):
            count[i] = converted_sequence.count(str(i))
        norm_count = {i:count/len(sequence) for i, count in count.items()}
        comp[prop] = norm_count
    return comp


def ctd_transition(sequence):
    sequence = sequence.upper()
    total_trans = {}
    for prop, values in properties.items():
        trans_values = {''.join(map(str, triad)): 0 for triad in product(range(1, 4), repeat=2)}
        convert = str_to_num(sequence, prop)
        for i in range(len(convert)-1):
            trans_values[convert[i:i+2]] += 1
        total_trans[prop] = trans_values
    return total_trans

#should devide by length to normalize
def ctd_distribution(sequence):
    """
    for every property there are 3 categories. this function calculates where the first, the first 25%, 50%, 75% and 100% of
    this category are in the protein sequence
    :param sequence:
    :return:
    """
    sequence = sequence.upper()
    total_distr = {}
    for prop, values in properties.items():
        convert = str_to_num(sequence, prop)
        distr = {}
        for i in set(values.values()):
            count = convert.count(str(i))
            occurences = [i for i, number in enumerate(convert)]
            distr["0.01"] = occurences[0]
            for j in [0.25, 0.50, 0.75, 1]:
                distr[str(j)] = occurences[math.ceil(count*j)-1]
        total_distr[prop] = distr
    return distr


print(ctd_distribution('AAA'))