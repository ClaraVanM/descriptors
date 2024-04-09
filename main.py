import pandas as pd

import process_file
import shape
import depth_comp

# <editor-fold desc="workflow">

"""
1) select neigborhood
2) calculate descriptors
        composition
        CTD
        ammPseAAC
        pseAAC
        autocorrelation
        sequence order
        buriedness
        deepness
        exposed residues
"""
# </editor-fold>

process_file.get_cavity_atoms("1B3X.pdb", "1B3X_out")
cavity = process_file.load_pdb("1B3X_neighbor.pdb")
protein = process_file.load_pdb("1GNY.pdb")
#center data
cavity[["x","y","z"]] = cavity[["x","y","z"]] - shape.COG(protein)
protein[["x","y","z"]] = protein[["x","y","z"]] - shape.COG(protein)

axis = shape.find_cavity_axis(cavity)
cavity_pr = shape.projection(cavity, axis)
#depth is first descriptor
df, depth = shape.add_buriedness(cavity,cavity_pr, axis)
narrow = shape.narrowness(df[df['buriedness'] ==1][['x','y','z']], axis, shape.COG(cavity))
#second descriptor
l_nar = shape.list_narrowness(df, axis, shape.COG(cavity))
#add distances of redidues to axis to df
df = shape.residue_dist_from_axis(df, axis)

#next use df to collect composition of protein based on buriedness, so per buriedness say how many of each AA?
df.to_csv('dataset')
l_nar = pd.DataFrame(l_nar, columns=['Values'])
l_nar.to_csv('list')
depth_comp.AA_per_buriedness(df)
