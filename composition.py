from itertools import product


# list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
groups = {"A":1,"G":1,"V":1,"I":2,"L":2,"F":2,"P":2,"Y":3,"M":3,"T":3,"S":3,"H":4,"N":4,"Q":4,"W":4, "R":5,"K":5, "D":6,"E":6, "C":7}


def aa_composition(sequence):
    sequence = sequence.upper()
    aa = {a:(sequence.count(a)/len(sequence)) for a in amino_acids}
    return aa


def dipeptide_composition(sequence):
    #drawback: does not count overlaps, 'AAA' is 1 count for 'AA'
    sequence = sequence.upper()
    dipeptide_composition= {i+j : sequence.count(i+j)/(len(sequence)/2) for i in amino_acids for j in amino_acids}
    return dipeptide_composition


def tripeptide_composition(sequence):
    sequence = sequence.upper()
    tripep_comp = {i+j+h : sequence.count(i+j+h)/(len(sequence)/3) for i in amino_acids for j in amino_acids for h in amino_acids}
    return tripep_comp


def conjoint_triad(sequence):
    sequence = sequence.upper()
    for i in sequence:
        sequence = sequence.replace(i, str(groups[i]))
    possible_triads = [''.join(map(str,triad)) for triad in product(range(1,8), repeat=3)]
    conjoint = {triad:0 for triad in possible_triads}
    for i in range(len(sequence)-2):
        conjoint[sequence[i:i+3]] +=1
    #normalize
    norm_conjoint = {triad:count/(len(sequence)-2) for triad, count in conjoint.items()}
    return norm_conjoint

