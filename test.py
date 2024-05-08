import pandas as pd
import distance_from_ligand
import autocorrelation
import process_file
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv("descriptors_ec3.2.1.csv")

cavity_file = "1O8S_neighbor.pdb"
cavity, ligand = process_file.load_pdb(cavity_file)


