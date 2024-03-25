import pymolPy3
import os
from PyBioMed.PyProtein import PyProtein

pm = pymolPy3.pymolPy3(0)
pm(f"load 1B3X.pdb")
pm(f"load 1B3X_out/pockets/pocket1_atm.pdb")
pm(f"select cavity, pocket1_atm")
pm(f"select prot, 1B3X")
pm(f"select neighborhood, prot near_to 5 of cavity")
pm(f"save neighbor.pdb, neighborhood")

pr = PyProtein.PyProtein("AADDDQDA")
#print(pr.GetDPComp())