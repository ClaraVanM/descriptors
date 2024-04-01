from aaindex import aaindex1

#list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

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


print(tripeptide_composition('aaa'))
