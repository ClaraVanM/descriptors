import os

import matplotlib
import pandas as pd
import numpy as np
from skspatial.objects import Line
from skspatial.objects import Plane
from skspatial.objects import Sphere
from scipy.spatial.distance import cdist

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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


cavity = load_pdb("cavity.pdb")
protein = load_pdb("1B3X.pdb")
#center data
mean = COG(protein)
protein = protein - mean
centered_cavity = cavity - mean
#project cavity on line
projected_cavity = projection(centered_cavity, protein)
#calculate distance from cog prot
projected_cavity["distance"] = np.linalg.norm(projected_cavity - mean, axis=1)
projected_cavity["index"] = projected_cavity.index
#sort distances to get highest distance first
projected_cavity = projected_cavity.sort_values(by=['distance'], ascending=True)
#deepness of cavity (max-min distance)
deepness = projected_cavity['distance'].max() - projected_cavity['distance'].min()
#make 10 intervals
buriedness = list()
jumps = deepness/5
for i in range(5):
    selection_index = projected_cavity[(projected_cavity["distance"] >= projected_cavity['distance'].min()+jumps*(i)) & (projected_cavity["distance"] <= projected_cavity['distance'].min()+jumps*(i+1))]['index']
    selection_cavity = centered_cavity.iloc[selection_index]
    buriedness.append(selection_cavity)
#project points on plane
#make the plane
plane1 = Plane(point = COG(protein), normal= (COG(centered_cavity) - COG(protein)))
#project the points
plane_projection = list()
for i in buriedness:
    plane_pr = pd.DataFrame(columns=["x", "y", "z"])
    for point in i.to_numpy():
        plane_pr.loc[len(plane_pr)] = plane1.project_point(point)
    plane_projection.append(plane_pr)


#find the most appropriate line that represents the middle of the cavity.
# go from center of cavity and make a line with every point, the two neighboring lines with highest angle between them represent the opening
#can be done in 3D circle,
sphere1 = Sphere(point=COG(centered_cavity), radius=15)
projection_sphere = pd.DataFrame(columns=["x", "y", "z"])
for point in centered_cavity.to_numpy():
    pr = sphere1.project_point(point)
    projection_sphere.loc[len(projection_sphere)] = pr
grid_sphere = Sphere(point=COG(centered_cavity), radius=15).to_points(n_angles = 20).unique()
#compute pairwise distance
#output is for every entry of grid_sphere an array of distances to every projection_sphere point
distances = cdist(grid_sphere, projection_sphere)
#filter distances with threshold and sum remaining number
distances = np.sum(distances < 4, axis=1)
df = pd.DataFrame(grid_sphere, columns=['x','y','z'])
df['distance'] = distances
df = df[df['distance'] == 0]
v1 =  COG(df[['x','y','z']]) -COG(centered_cavity)
print(v1)
line1 = Line(point=COG(centered_cavity), direction=v1)

#https://scikit-spatial.readthedocs.io/en/stable/objects/line.html
"""matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(protein['x'], protein['y'], protein['z'], c='k', s=0.5)
ax.scatter(projected_cavity['x'], projected_cavity['y'], projected_cavity['z'],c='g')
ax.scatter(buriedness[0]['x'], buriedness[0]['y'], buriedness[0]['z'], c='b' )
ax.scatter(buriedness[1]['x'], buriedness[1]['y'], buriedness[1]['z'], c='y' )
ax.scatter(buriedness[2]['x'], buriedness[2]['y'], buriedness[2]['z'], c='r' )
ax.scatter(buriedness[3]['x'], buriedness[3]['y'], buriedness[3]['z'], c='b' )
ax.scatter(buriedness[4]['x'], buriedness[4]['y'], buriedness[4]['z'], c='y' )
plane1.plot_3d(ax, color='r')
plt.show()"""

matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(projection_sphere['x'], projection_sphere['y'], projection_sphere['z'], s=1)
ax.scatter(df['x'], df['y'], df['z'], s=1)
line1.plot_3d(ax)
plt.show()






#1) make line from center protein to center cavity
#2) find furthest point of cavity following that line
#3) take all other points with same distance +- x
#4) project those points on 2D plane
#5) look at that projection and find a good 2D curve to fit and take standard dev or diameter

#calculate narrowness at multiple buriedness intervalls
#first take whole length and devide that into desired number of intervals
#for every interval calculate the narowness
#narowness can be measured