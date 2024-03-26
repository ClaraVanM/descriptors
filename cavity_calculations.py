import os

import matplotlib
import pandas as pd
import numpy as np
from skspatial.objects import Line
from skspatial.objects import Plane
from skspatial.objects import Sphere
from skspatial.objects import Vector
from skspatial.objects import Point
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


def find_cavity_axis(centered_cavity):
    # find the most appropriate line that represents the middle of the cavity.
    # go from center of cavity and make a line with every point, the two neighboring lines with highest angle between them represent the opening
    # can be done in 3D circle,
    sphere1 = Sphere(point=COG(centered_cavity), radius=15)
    projection_sphere = pd.DataFrame(columns=["x", "y", "z"])
    for point in centered_cavity.to_numpy():
        pr = sphere1.project_point(point)
        projection_sphere.loc[len(projection_sphere)] = pr
    grid_sphere = Sphere(point=COG(centered_cavity), radius=15).to_points(n_angles=20).unique()
    # compute pairwise distance
    # output is for every entry of grid_sphere an array of distances to every projection_sphere point
    distances = cdist(grid_sphere, projection_sphere)
    # filter distances with threshold and sum remaining number
    distances = np.sum(distances < 4, axis=1)
    df = pd.DataFrame(grid_sphere, columns=['x', 'y', 'z'])
    df['distance'] = distances
    df = df[df['distance'] == 0]
    vector = COG(df[['x', 'y', 'z']]) - COG(centered_cavity)
    cavity_axis = Line(point=COG(centered_cavity), direction=vector)
    return cavity_axis


def projection(cavity, cavity_axis):
    """
    projects all cavity points on cavity_axis
    :param cavity: x, y and z coordinates of cavity atoms
    :param cavity_axis: axis determined by find_cavity_axis()
    :return: pandas dataframe with projected coordinates
    """
    projection_cavity = pd.DataFrame(columns=['x', 'y', 'z'])
    for point in cavity.to_numpy():
        projection_cavity.loc[len(projection_cavity)] = cavity_axis.project_point(point)
    return projection_cavity


def add_buriedness(cavity, projection,axis):
    #calculate distances and devide based on those lengths
    projection["distance"] = np.dot(cavity - COG(cavity), axis.direction)
    projection["index"] = projection.index
    projection = projection.sort_values(by=['distance'], ascending=True)
    # deepness of cavity (max-min distance)
    deepness = projection['distance'].max() - projection['distance'].min()
    # make 5 intervals
    jumps = deepness / 5
    cavity['buriedness'] = -1
    for i in range(5):
        selection_index = projection[
            (projection["distance"] >= projection['distance'].min() + jumps * (i)) & (
                        projection["distance"] <= projection['distance'].min() + jumps * (i + 1))]['index']
        cavity.loc[selection_index, 'buriedness'] = i
    return cavity, deepness


def narrowness(specific_buriedness, cavity_axis, cavity_center):
    # project points on plane
    # make the plane
    plane1 = Plane(point=cavity_center, normal=cavity_axis.direction)
    # project the points
    plane_pr = pd.DataFrame(columns=["x", "y", "z"])
    for point in specific_buriedness.to_numpy():
        plane_pr.loc[len(plane_pr)] = plane1.project_point(point)
    center = Point(cavity_center)
    plane_pr["distance"] = [center.distance_point(x) for x in plane_pr.to_numpy()]
    #shortest distance gives radius of biggest circle that fits
    return np.min(plane_pr["distance"]), plane_pr


def list_narrowness(df_burriednes, cavity_axis, cavity_center):
    values = [0,1,2,3,4]
    l = list()
    for i in values:
        l.append(narrowness(df_burriednes[df_burriednes["buriedness"] == values[i]][['x','y','z']], cavity_axis, cavity_center)[0])
    boolean = False
    # cut of from back of list to front when value drops under 1.5
    #then cavity is not open anymore
    for i in reversed(l):
        if i > 1.5 and not boolean:
            continue
        else:
            boolean = True
            l[l.index(i)] = 0
    return l





########################################################################################################################
#perform code
cavity = load_pdb("1GNY_cavity.pdb")
protein = load_pdb("1GNY.pdb")
#center data
mean = COG(protein)
protein = protein - mean
centered_cavity = cavity - mean
cavity_mean = COG(centered_cavity)

axis = find_cavity_axis(centered_cavity)
cavity_pr = projection(centered_cavity, axis)
df, depth = add_buriedness(centered_cavity,cavity_pr, axis)
narrow, df_narrow = narrowness(df[df['buriedness'] ==1][['x','y','z']], axis, cavity_mean)
list_narrowness(df, axis, cavity_mean)


########################################################################################################################
#visualize
matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(df['x'], df['y'], df['z'], s=1.5,c = df['buriedness'])
ax.scatter(df_narrow['x'], df_narrow['y'], df_narrow['z'])
axis.plot_3d(ax)
plt.show()