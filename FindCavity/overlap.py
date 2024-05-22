import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def get_max_min(point_cloud):
    """

    :param point_cloud: coordinates of pdb structure in numpy array
    :return: max and min value
    """
    max_values = np.max(point_cloud, axis=0)
    min_values = np.min(point_cloud, axis=0)
    return max_values, min_values


def intersection(xmin, xmax, ymin, ymax):
    intersect = 0
    if xmin > ymin and xmin < ymax:
        if xmax > ymax:
            intersect = ymax-xmin
        elif xmax < ymax:
            intersect = xmax-xmin
    elif xmin < ymin and xmax > ymin:
        if ymax > xmax:
            intersect =xmax - ymin
        elif xmax > ymax:
            intersect = ymax-ymin
        else:
            intersect = 0
    elif xmin == ymin:
        if xmax == ymax:
            intersect = xmax - xmin
        elif xmax > ymax:
            intersect = ymax - ymin
        else:
            intersect = xmax - xmin
    return intersect


def overlap(intercept, ligand_max, ligand_min):
    return intercept / (ligand_max-ligand_min)


def total_overlap(ligand, cavity):
    ligand, cavity = pca(ligand,cavity)
    cavity_max, cavity_min = get_max_min(cavity)
    ligand_max, ligand_min = get_max_min(ligand)
    x_intercept = intersection(cavity_min["x"], cavity_max["x"], ligand_min["x"], ligand_max["x"])
    y_intercept = intersection(cavity_min["y"], cavity_max["y"], ligand_min["y"], ligand_max["y"])
    z_intercept = intersection(cavity_min["z"], cavity_max["z"], ligand_min["z"], ligand_max["z"])
    x_overlap = overlap(x_intercept, ligand_max["x"], ligand_min["x"])
    y_overlap = overlap(y_intercept, ligand_max["y"], ligand_min["y"])
    z_overlap = overlap(z_intercept, ligand_max["z"], ligand_min["z"])
    if any(x==0 for x in [x_overlap, y_overlap, z_overlap]):
        total_ov = 0
    else:
        total_ov = (x_overlap + y_overlap + z_overlap) / 3
    return total_ov


def pca(ligand, cavity):
    pca = PCA(n_components=3)
    pca.fit(cavity)
    ligand_transformed = pd.DataFrame(pca.transform(ligand), columns = ["x", "y", "z"])
    cavity_transformed = pd.DataFrame(pca.transform(cavity), columns = ["x", "y", "z"])
    return ligand_transformed, cavity_transformed