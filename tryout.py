import matplotlib
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

centered_cavity = None
tri = Delaunay(centered_cavity.to_numpy())
def edge_length(tri, i):
    p1, p2, p3, p4 = tri.points[tri.simplices[i]]
    return np.linalg.norm(p1-p2), np.linalg.norm(p2-p3), np.linalg.norm(p3-p4), np.linalg.norm(p3-p1), np.linalg.norm(p4-p1), np.linalg.norm(p4-p2)

"""short_simplices = [i for i in range(len(tri.simplices)) if all(length <=5 for length in edge_length(tri,i))]
tri.simplices = tri.simplices[short_simplices]
hull = tri.convex_hull"""

hull = ConvexHull(centered_cavity)
"""for i in hull.simplices:
    print(hull.points)"""
print(hull.points)
#short_simplices = [i for i in range(len(hull.simplices)) if all(length <=5 for length in edge_length(hull,i))]

matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
mesh = Poly3DCollection([centered_cavity.to_numpy()[s] for s in hull.simplices], alpha=0.25, facecolors='cornflowerblue', linewidths=1, edgecolors='lightsteelblue')
ax.add_collection3d(mesh)
ax.set_xlim(-20,20)
ax.set_ylim(-20,20)
ax.set_zlim(-20,20)
plt.show()