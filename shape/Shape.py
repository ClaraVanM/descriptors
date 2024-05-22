import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skspatial.objects import Line
from skspatial.objects import Plane
from skspatial.objects import Sphere
from skspatial.objects import Point
from scipy.spatial.distance import cdist

class Shape:

    AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    residues = ['O', 'C', 'N', 'S', 'H']

    def __init__(self, cavity, ligand):
        self.center = self.COG(cavity)
        self.input = cavity
        self.ligand = ligand
        self.axis = self.find_cavity_axis()
        self.cavity, self.depth = self.get_buriedness()


    def get_buriedness(self):
        df = self.residue_dist_from_axis()
        cavity_projection = self.projection(df)
        df, depth = self.add_buriedness(df, cavity_projection)
        return df, depth

    def getDescriptors(self):
        descripors = dict()
        list_narrowness = self.list_narrowness()
        AA_comp = self.AA_per_buriedness()
        exposed = self.exposed_aa(list_narrowness)
        descripors['depth'] = self.depth
        for i in range(len(list_narrowness)):
            descripors['narrowness_' + str(i)] = list_narrowness[i]
        for values in AA_comp.values():
            descripors.update(values)
        descripors.update(exposed)
        return descripors

    @staticmethod
    def COG(df):
        """

        :return: center of the xyz coordinates
        """
        return df[["x", "y", "z"]].mean()

    def find_cavity_axis(self):
        """

        :return:axis that goes through middle of cavity
        """
        # find the most appropriate line that represents the middle of the cavity.
        sphere1 = Sphere(point=Shape.COG(self.ligand), radius=15)
        projection_sphere = pd.DataFrame(columns=["x", "y", "z"])
        for point in self.input[["x", "y", "z"]].to_numpy():
            pr = sphere1.project_point(point)
            projection_sphere.loc[len(projection_sphere)] = pr
        grid_sphere = Sphere(point=Shape.COG(self.ligand), radius=15).to_points(n_angles=30).unique()
        # compute pairwise distance
        distances = cdist(grid_sphere, projection_sphere)
        # filter distances with threshold and sum remaining number, so only pr close enough are taken into account (<4 in neighbourhood), and are summed together, --> distance = distance of grid points to all neighbouring pr points
        distances = np.sum(distances < 2, axis=1)
        df = pd.DataFrame(grid_sphere, columns=['x', 'y', 'z'])
        df['distance'] = distances
        df = df[df['distance'] == 0]
        # do clustering and extract the point group with most members as ultimate cavity opening
        df = Shape.cluster(df[['x', 'y', 'z']])
        vector = Shape.COG(df[['x', 'y', 'z']]) - self.center
        cavity_axis = Line(point=self.center, direction=vector)
        return cavity_axis

    @staticmethod
    def cluster(opening):
        # do clustering on points that represent opening of cavity in find_cavity_axis.
        opening = opening.reset_index(drop=True)
        dist = cdist(opening, opening, 'euclidean')
        point_pairs = np.argwhere(dist < 6)
        clusters = []
        for pairs in point_pairs:
            # if point is with itself
            if pairs[0] == pairs[1]:
                b = False
                for list in clusters:
                    # check if that point already clustered if so , b = true
                    if pairs[0] in list:
                        b = True
                # if b = true continue, otherwise add
                if not b:
                    clusters.append([pairs[0]])
            # now add elements to their cluster, skip doubles by saying <
            elif pairs[0] < pairs[1]:
                for list in clusters:
                    if pairs[0] in list and not pairs[1] in list:
                        list.append(pairs[1])
        # put indexes of the points together with their cluster label in a dictionary
        cluster_number = {}
        for group, cluster in enumerate(clusters):
            for i in cluster:
                cluster_number[i] = group
        # add to pandas df
        opening['cluster'] = cluster_number
        opening = opening[opening['cluster'] == opening['cluster'].value_counts().idxmax()]
        return opening[['x', 'y', 'z']]

    def residue_dist_from_axis(self):
        df = self.input.copy()
        df.loc[:,'dist_from_axis'] = float(0)
        for index, row in df.iterrows():
            df.loc[index, 'dist_from_axis'] = self.axis.distance_point(Point([row["x"], row["y"], row["z"]]))
        # choose closest point of aa to axis as representative for aa
        subset_index = df.groupby('AA_number')['dist_from_axis'].idxmin()
        out = df.loc[subset_index].copy()
        out.reset_index(drop=True, inplace=True)
        return out

    def projection(self, df):
        """
        projects all cavity points on cavity_axis
        """
        projection_cavity = pd.DataFrame(columns=['x', 'y', 'z'])
        for point in df[['x', 'y', 'z']].to_numpy():
            projection_cavity.loc[len(projection_cavity)] = self.axis.project_point(point)
        return projection_cavity

    def add_buriedness(self, df, projection):
        """
        split the projection points into 5 intervals, depending on deepness
        :return:pandas dataframe of coordinates of cavity, with extra column of buriedness which has a value from 1 to 5
        5 being the closest to the cavity opening. + deepness of the cavity
        """
        # calculate distances and devide based on those lenths
        projection["distance"] = np.dot(df[["x", 'y', "z"]] - self.center, self.axis.direction)
        # deepness of cavity (max-min distance)
        deepness = projection['distance'].max() - projection['distance'].min()
        # make 5 intervals
        jumps = deepness / 5
        cavity = df.copy()
        cavity['buriedness'] = -1
        for i in range(6):
            selection_index = projection[
                (projection["distance"] >= projection['distance'].min() + jumps * (i)) & (
                        projection["distance"] < projection['distance'].min() + jumps * (i + 1))].index
            cavity.loc[selection_index, 'buriedness'] = i
        selection = projection[projection['distance']>= projection['distance'].min()+jumps*5].index
        cavity.loc[selection, 'buriedness'] = 4
        return cavity, deepness

    def narrowness(self, selection):
        """

        :return: radius of biggest sphere that fits without having points in the surface
        """
        # project points on plane
        # make the plane
        plane1 = Plane(point=self.center, normal=self.axis.direction)
        # project the points
        plane_pr = pd.DataFrame(columns=["x", "y", "z"])
        for point in selection[["x", "y", "z"]].to_numpy():
            plane_pr.loc[len(plane_pr)] = plane1.project_point(point)
        center = Point(self.center)
        plane_pr["distance"] = [center.distance_point(x) for x in plane_pr.to_numpy()]
        # shortest distance gives radius of biggest circle that fits
        return np.min(plane_pr["distance"])

    def list_narrowness(self):
        """

        :return: list of narrowness for every degree of burriedness
        """
        values = [0, 1, 2, 3, 4]
        l = list()
        for i in values:
            l.append(self.narrowness(self.cavity[self.cavity["buriedness"] == values[i]]))
        boolean = False
        mean = np.nanmean(l)
        for i in range(len(l)):
            if np.isnan(l[i]):
                l[i] = mean
        # cut of from back of list to front when value drops under 1.5
        # then cavity is not open anymore
        for i in reversed(l):
            if i > 1.5 and not boolean:
                continue
            else:
                boolean = True
                l[l.index(i)] = 0
        return l

    def AA_per_buriedness(self):
        """

        :return: dictionary with frequency of amino acids per buriedness level
        """
        total_count = {}
        for i in self.cavity["buriedness"].unique():
            df_temp = self.cavity[self.cavity["buriedness"] == i]
            aa_count = df_temp['AA'].value_counts()
            total_count[i] = {'buried_' + aa + str(i): aa_count.get(aa, 0) for aa in Shape.AA_symb}
        return total_count

    def exposed_aa(self, narrow_list):
        """
        :param narrow_list: list of narrowness per buriedness level
        :return: closest 30% of residues to cavity
        """
        if narrow_list.count(0) != 0:
            index = narrow_list.count(0) - 1
        else:
            index = 0
        depth = list(range(index, 5))
        df = self.cavity[self.cavity['buriedness'].isin(depth)]
        # take closest 5 of aa to axes
        df = df.sort_values(by='dist_from_axis').head(int(len(df) * 0.3))
        d_count = {residue: 0 for residue in Shape.residues}
        for value in df['atom']:
            atom = value[0]
            d_count[atom] += 1
        return d_count
