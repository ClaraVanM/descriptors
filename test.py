import pandas as pd
import process_file
import shape
import depth_comp
import matplotlib
import matplotlib.pyplot as plt
from skspatial.objects import Sphere
from scipy.spatial.distance import cdist
import numpy as np
from skspatial.objects import Line


def outliers(cavity):
    q1 = np.percentile(cavity, 25, axis=0)
    q3 = np.percentile(cavity, 75, axis=0)
    iqr = q3-q1
    lower_bound = q1 - 5 * iqr
    upper_bound = q3 +5*iqr
    outliers = np.any((cavity < lower_bound) | (cavity > upper_bound), axis=1)
    cleaned_cavity = cavity[~outliers]
    return cleaned_cavity

cavity_file = "1O8S_neighbor.pdb"
cavity, ligand = process_file.load_pdb(cavity_file)
cavity[["x","y","z"]] = outliers(cavity[["x","y","z"]])
cavity.dropna(inplace=True)
cavity['index'] == cavity.index
axis = shape.find_cavity_axis(cavity, ligand)
projection = shape.projection(cavity, axis)
cavity, depth = shape.add_buriedness(cavity, projection, axis)
"""narrowl = shape.list_narrowness(cavity, axis, shape.COG(cavity))
print(narrowl)"""



matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(cavity['x'], cavity['y'], cavity['z'],s=1)
ax.scatter(ligand['x'], ligand['y'], ligand['z'], s=1)
ax.scatter(projection['x'], projection['y'], projection['z'])
plt.show()


"""Auto correlation not within sequence but ranking amino acids based on distance from cog ligand

1) buriedness
2) distance from ligand"""