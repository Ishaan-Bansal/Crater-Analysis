import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats

filename = 'Crater_STL_Files/02_09_2022_6Torr_test2.stl'
trimRadius = 250
displayRadius = 250

lowest_point, rotation_matrix, rotation_matrix_point = help.trimming_package(
    filename, trimRadius)
lowest_point = rotation_matrix_point @ lowest_point

# Trim around the lowest_point
new_mesh = Mesh.from_file(
    filename)
help.trimCircleGivenPoint(new_mesh, lowest_point, displayRadius)
mesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
mesh.remove_infinite_values()
# mesh.apply_transform(
#     trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
mesh.apply_transform(rotation_matrix)

mesh.show()

# Histogram of heights
heights = help.zPoints(mesh)
heights -= lowest_point[2]
baseline = stats.mode(heights)
plt.figure(1)
n, bins, patches = plt.hist(
    heights, 300, density=True, facecolor='b', alpha=0.75)
nLine, binsLine = help.histPlot(n, bins)
# plt.plot(binsLine, nLine, 'r--', linewidth=2)
baseline = bins[n.argmax()-1]
plt.axvline(x=baseline, color='g', linestyle='dashed')
# plt.show()

# Getting a slice of the mesh
y, z = help.slicer(mesh)
plt.figure(2)
plt.plot(y, z, 'bs', linewidth=2)
plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.axis('equal')
plt.show()
