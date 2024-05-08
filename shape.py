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


def COG(data):
    """

    :param data: processed pdb file
    :return: center of the xyz coordinates
    """
    return data[["x","y","z"]].mean()


def find_cavity_axis(cavity, ligand):
    """
    makes a sphere with as center the cavity and then projects all cavity points on the shell of the sphere.
    make a grid sphere with poinst equaly distributed over the surface
    pairwise distances between all projected points and the grid points
    take only distances smaller then 4 and sum them per grid point
    grid points with total distance 0 have no projection points in their approximation
    take the middle of these grid points and make a vector from center of cavity to center of selected grid points
    :param cavity: processed pdb file of cavity
    :return:axis that goes through middle of cavity
    """
    # find the most appropriate line that represents the middle of the cavity.
    sphere1 = Sphere(point=COG(ligand), radius=15)
    projection_sphere = pd.DataFrame(columns=["x", "y", "z"])
    for point in cavity[["x","y","z"]].to_numpy():
        pr = sphere1.project_point(point)
        projection_sphere.loc[len(projection_sphere)] = pr
    grid_sphere = Sphere(point=COG(ligand), radius=15).to_points(n_angles=30).unique()
    # compute pairwise distance
    distances = cdist(grid_sphere, projection_sphere)
    # filter distances with threshold and sum remaining number, so only pr close enough are taken into account (<4 in neighbourhood), and are summed together, --> distance = distance of grid points to all neighbouring pr points
    distances = np.sum(distances < 2, axis=1)
    df = pd.DataFrame(grid_sphere, columns=['x', 'y', 'z'])
    df['distance'] = distances
    df = df[df['distance'] == 0]
    #do clustering and extract the point group with most members as ultimate cavity opening
    df= cluster(df[['x','y','z']])
    vector = COG(df[['x', 'y', 'z']]) - COG(cavity)
    cavity_axis = Line(point=COG(cavity), direction=vector)
    return cavity_axis


def projection(cavity, cavity_axis):
    """
    projects all cavity points on cavity_axis
    :param cavity: pdb file of cavity
    :param cavity_axis: axis determined by find_cavity_axis()
    :return: pandas dataframe with projected coordinates
    """
    projection_cavity = pd.DataFrame(columns=['x', 'y', 'z'])
    for point in cavity[['x','y','z']].to_numpy():
        projection_cavity.loc[len(projection_cavity)] = cavity_axis.project_point(point)
    return projection_cavity


def add_buriedness(cavity, projection, axis):
    """
    split the projection points into 5 intervals, depending on deepness
    :param cavity: processed pdb file of cavity
    :param projection: projection of cavity, output of projection function
    :param axis: axis used for the projection
    :return:pandas dataframe of coordinates of cavity, with extra column of buriedness which has a value from 1 to 5
    5 being the closest to the cavity opening. + deepness of the cavity
    """
    #calculate distances and devide based on those lenths
    projection["distance"] = np.dot(cavity[["x",'y',"z"]] - COG(cavity), axis.direction)
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
    """

    :param specific_buriedness: subset of dataframe based on burriedness
    :param cavity_axis: axis used for projection/middle of cavity
    :param cavity_center: center of cavity
    :return: radius of biggest sphere that fits without having points in the surface
    """
    # project points on plane
    # make the plane
    plane1 = Plane(point=cavity_center, normal=cavity_axis.direction)
    # project the points
    plane_pr = pd.DataFrame(columns=["x", "y", "z"])
    for point in specific_buriedness[["x","y","z"]].to_numpy():
        plane_pr.loc[len(plane_pr)] = plane1.project_point(point)
    center = Point(cavity_center)
    plane_pr["distance"] = [center.distance_point(x) for x in plane_pr.to_numpy()]
    #shortest distance gives radius of biggest circle that fits
    return np.min(plane_pr["distance"])


def list_narrowness(df_burriednes, cavity_axis, cavity_center):
    """

    :param df_burriednes: pandas dataframe with burriedness column
    :param cavity_axis: axis of cavity
    :param cavity_center: center of cavity
    :return: list of narrowness for every degree of burriedness
    """
    values = [0,1,2,3,4]
    l = list()
    for i in values:
        l.append(narrowness(df_burriednes[df_burriednes["buriedness"] == values[i]][['x','y','z']], cavity_axis, cavity_center))
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


def residue_dist_from_axis(df, axis):
    df['dist_from_axis'] = float(0)
    for index, row in df.iterrows():
        df.loc[index, 'dist_from_axis'] = axis.distance_point(Point([row["x"], row["y"], row["z"]]))
    return df

def cluster(opening):
    #do clustering on points that represent opening of cavity in find_cavity_axis.
    opening = opening.reset_index(drop=True)
    dist = cdist(opening, opening,'euclidean')
    point_pairs = np.argwhere(dist<6)
    clusters = []
    for pairs in point_pairs:
        #if point is with itself
        if pairs[0] == pairs[1]:
            b = False
            for list in clusters:
                #check if that point already clustered if so , b = true
                    if pairs[0] in list:
                        b = True
            # if b = true continue, otherwise add
            if not b:
                clusters.append([pairs[0]])
        #now add elements to their cluster, skip doubles by saying <
        elif pairs[0] < pairs[1]:
            for list in clusters:
                if pairs[0] in list and not pairs[1] in list:
                    list.append(pairs[1])
    #put indexes of the points together with their cluster label in a dictionary
    cluster_number = {}
    for group, cluster in enumerate(clusters):
        for i in cluster:
            cluster_number[i] = group
    #add to pandas df
    opening['cluster'] = cluster_number
    opening = opening[opening['cluster'] == opening['cluster'].value_counts().idxmax()]
    return opening[['x','y','z']]

