import os
import pandas as pd
import numpy as np
from skspatial.objects import Line
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_pdb(file):
    assert os.path.isfile(file), "No valid file."
    df = pd.DataFrame(columns=['x', 'y', 'z'])
    with open(file, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('ATOM'):
                df.loc[len(df)] = line.split()[6:9]
            line = f.readline()
    df[['x', 'y', 'z']] = df[['x','y', 'z']].astype(float)
    return df


def COG(data):
    return data.mean()


def projection(cavity, protein):
    #v1 = vector pointing from cog prot to cog cavity
    v1 = COG(cavity) - COG(protein)
    line1 = Line(point=COG(protein), direction=v1)
    #project
    projection_cavity = pd.DataFrame(columns=['x', 'y', 'z'])
    for point in cavity.to_numpy():
        projection_cavity.loc[len(projection_cavity)] = line1.project_point(point)
    return projection_cavity


cavity = load_pdb("neighbor.pdb")
protein = load_pdb("1B3X.pdb")
#center
mean = COG(protein)
protein = protein - mean
centered_cavity = cavity - mean
projected_cavity = projection(centered_cavity, protein)
projected_cavity["distance"] = np.linalg.norm(projected_cavity - mean, axis=1)
projected_cavity["index"] = projected_cavity.index
projected_cavity = projected_cavity.sort_values(by=['distance'], ascending=False)
cavity_opening_indexes = projected_cavity["index"].head(10).tolist()
cavity_opening = centered_cavity.iloc[cavity_opening_indexes]
#next step is projecting them on the plane
print(cavity_opening)

#https://scikit-spatial.readthedocs.io/en/stable/objects/line.html
"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(protein['x'], protein['y'], protein['z'], c='y', s=0.5)
ax.scatter(centered_cavity['x'], centered_cavity['y'], centered_cavity['z'], c='r')
ax.scatter(projected_cavity['x'], projected_cavity['y'], projected_cavity['z'],c='g')
plt.show()"""


#1) make line from center protein to center cavity
#2) find furthest point of cavity following that line
#3) take all other points with same distance +- x
#4) project those points on 2D plane
#5) look at that projection and find a good 2D curve to fit and take standard dev or diameter
