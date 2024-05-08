import matplotlib.pyplot as plt
import distance_from_ligand
import matplotlib
from skspatial.objects import Sphere
from skspatial.objects import Line
from skspatial.objects import Plane
from skspatial.objects import Point
from scipy.spatial.distance import cdist
import pandas as pd
import numpy as np

import process_file
import shape

###############################################################################################################################################################################
##necessary data
cavity, ligand = process_file.load_pdb("8TC8_neighbor.pdb")
####################################################################################################################################################################################
##sphere projection
"""
sphere1 = Sphere(point=shape.COG(ligand), radius=15)
projection_sphere = pd.DataFrame(columns=["x", "y", "z"])
for point in cavity[["x","y","z"]].to_numpy():
    pr = sphere1.project_point(point)
    projection_sphere.loc[len(projection_sphere)] = pr
grid_sphere = Sphere(point=shape.COG(ligand), radius=15).to_points(n_angles=30).unique()
# compute pairwise distance
distances = cdist(grid_sphere, projection_sphere)
# filter distances with threshold and sum remaining number, so only pr close enough are taken into account (<4 in neighbourhood), and are summed together, --> distance = distance of grid points to all neighbouring pr points
distances = np.sum(distances < 2, axis=1)
df = pd.DataFrame(grid_sphere, columns=['x', 'y', 'z'])
df['distance'] = distances
df = df[df['distance'] == 0]
#do clustering and extract the point group with most members as ultimate cavity opening
df2= shape.cluster((df[['x','y','z']]))
vector = shape.COG(df2[['x', 'y', 'z']] - shape.COG(cavity))
cavity_axis = Line(point=shape.COG(cavity), direction=vector)

matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(projection_sphere["x"], projection_sphere["y"], projection_sphere["z"], s=1)
cavity_axis.plot_3d(ax)
ax.set_axis_off()
ax.patch.set_alpha(0)
plt.show()"""


####################################################################################################################################################################
#####plot buriedness
axis = shape.find_cavity_axis(cavity,ligand)
projection = shape.projection(cavity, axis)
buriedness, deepness = shape.add_buriedness(cavity, projection, axis)

"""matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(buriedness[buriedness["buriedness"] ==0]['x'], buriedness[buriedness["buriedness"] ==0]['y'], buriedness[buriedness["buriedness"] ==0]['z'],c="lightcoral", s=1 )
ax.scatter(buriedness[buriedness["buriedness"] ==1]['x'], buriedness[buriedness["buriedness"] ==1]['y'], buriedness[buriedness["buriedness"] ==1]['z'], c='cornflowerblue',s=1 )
ax.scatter(buriedness[buriedness["buriedness"] ==2]['x'], buriedness[buriedness["buriedness"] ==2]['y'], buriedness[buriedness["buriedness"] ==2]['z'], c='coral' ,s=1)
ax.scatter(buriedness[buriedness["buriedness"] ==3]['x'], buriedness[buriedness["buriedness"] ==3]['y'], buriedness[buriedness["buriedness"] ==3]['z'], c='springgreen',s=1 )
ax.scatter(buriedness[buriedness["buriedness"] ==4]['x'], buriedness[buriedness["buriedness"] ==4]['y'], buriedness[buriedness["buriedness"] ==4]['z'], c='sandybrown',s=1 )
ax.set_axis_off()
ax.patch.set_alpha(0)
plt.show()"""


##################################################################################################################################################################################
##narowness

"""# make the plane
plane1 = Plane(point=shape.COG(cavity), normal=axis.direction)
# project the points
plane_pr = pd.DataFrame(columns=["x", "y", "z"])
specific_buriedness = buriedness[buriedness["buriedness"]==4]
for point in specific_buriedness[["x","y","z"]].to_numpy():
    plane_pr.loc[len(plane_pr)] = plane1.project_point(point)
center = Point(shape.COG(cavity))
plane_pr["distance"] = [center.distance_point(x) for x in plane_pr.to_numpy()]
#shortest distance gives radius of biggest circle that fits
np.min(plane_pr["distance"])
matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(plane_pr["x"], plane_pr["y"], plane_pr["z"])
axis.plot_3d(ax)
ax.set_axis_off()
ax.patch.set_alpha(0)
plt.show()"""


#################################################################################
cavity = distance_from_ligand.dist_from_ligand(cavity, distance_from_ligand.COG(ligand))
cavity = distance_from_ligand.divide_cavity(cavity)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(cavity[cavity["group"] ==0]['x'], cavity[cavity["group"] ==0]['y'], cavity[cavity["group"] ==0]['z'],c="grey", s=1 )
ax.scatter(cavity[cavity["group"] ==1]['x'], cavity[cavity["group"] ==1]['y'], cavity[cavity["group"] ==1]['z'], c='coral',s=1 )
ax.scatter(cavity[cavity["group"] ==2]['x'], cavity[cavity["group"] ==2]['y'], cavity[cavity["group"] ==2]['z'], c='cornflowerblue' ,s=1)
ax.scatter(cavity[cavity["group"] ==3]['x'], cavity[cavity["group"] ==3]['y'], cavity[cavity["group"] ==3]['z'], c='springgreen',s=1 )
ax.scatter(cavity[cavity["group"] ==4]['x'], cavity[cavity["group"] ==4]['y'], cavity[cavity["group"] ==4]['z'], c='sandybrown',s=1 )
axis.plot_3d(ax)
ax.set_axis_off()
ax.patch.set_alpha(0)
plt.show()