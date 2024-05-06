import pandas as pd
import process_file
import shape
import depth_comp
import matplotlib
import matplotlib.pyplot as plt

protein_file = "1O8S.pdb"
fpocket= ("1O8S_out")
cavity_file = "1O8S_neighbor.pdb"
cavity = process_file.load_pdb(cavity_file)
axis = shape.find_cavity_axis(cavity)

matplotlib.use('TkAgg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(cavity["x"], cavity["y"], cavity["z"], s=1)
axis.plot_3d(ax)
plt.show()


