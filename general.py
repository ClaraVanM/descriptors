import numpy as np
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser

s1 = get_structure("1B3X.pdb")
s2 = get_structure("1GNY.pdb")
coords1, seq1 = get_residue_data(next(s1.get_chains()))
coords2, seq2 = get_residue_data(next(s2.get_chains()))
res = tm_align(coords1, coords2, seq1, seq2)

p = PDBParser()
struct = p.get_structure("1GNY3", "1GNY.pdb")
sr = ShrakeRupley()
sr.compute(struct, level="S")
print(struct)