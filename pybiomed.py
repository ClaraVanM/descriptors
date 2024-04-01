import cavity

from PyBioMed.PyProtein import PyProtein
from Bio.PDB import PDBParser
from Bio import SeqIO
from PyBioMed.PyProtein import AAComposition, CTD


def get_sequence_from_pdb(pdb_file):
    sequence = ""



#following code gives me: WRGNASGSTSHSGIXXXXXXXXFXXXGDGVGAVFDIXXXXXXXXXXXXXXXXXXXXXXXXXXXXQIFAQLKEDWSKGEWDCXXXXXXXXXXXXXXXXXXXXEXXXKFXXXXRDVQVGIQAK
"""cavity_atoms = cavity.get_cavity_atoms("1GNY.pdb", "1GNY_out")
for record in SeqIO.parse("1GNY_cavity.pdb", "pdb-atom"):
    print(record.seq)"""


#following code, gives me WRGNASGSTSHSGIFGDGVGAVFDIQIFAQLKEDWSKGEWDCEKFRDVQVGIQAK
p = PDBParser()
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
struct = p.get_structure("1GNY3", "1GNY_cavity.pdb")
for model in struct:
    for chain in model:
        seq=[]
        for residue in chain:
            seq.append(d3to1[residue.resname])
       #print(''.join(seq))


pr = PyProtein.PyProtein(seq)
aac = AAComposition.CalculateAAComposition(seq)
ctd = CTD.CalculateC(seq)
print(pr.getDPComp())
#???faults in pybiomed package code

