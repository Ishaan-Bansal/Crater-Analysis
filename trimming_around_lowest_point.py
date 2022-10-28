import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats

lowest_point = help.lowest_point_file(
    r'Crater_STL_Files/02_09_2022_250Torr_test3.stl')

# Trim around the lowest_point
new_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
help.trimCircleGivenPoint(new_mesh, lowest_point, 200)
mesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
mesh.apply_transform(
    trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
mesh.apply_transform(help.rotation_matrix_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl'))

mesh.show()

# Histogram of heights
heights = help.zPoints(mesh)
heights -= lowest_point[2]
baseline = stats.mode(heights)
plt.figure(1)
n, bins, patches = plt.hist(
    heights, 300, density=True, facecolor='b', alpha=0.75)
nLine, binsLine = help.histPlot(n, bins)
plt.plot(binsLine, nLine, 'r--', linewidth=2)
baseline = bins[n.argmax()-1]
plt.axvline(x=baseline, color='g', linestyle='dashed')
# plt.show()

# Getting a slice of the mesh
y, z = help.slicer(mesh)
plt.figure(2)
plt.plot(y, z, 'bs', linewidth=2)
plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.show()
