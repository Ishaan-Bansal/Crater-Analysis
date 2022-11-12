import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats

# Write the relative path to the file you want to load
filename = 'Crater_STL_Files/02_09_2022_250Torr_test3.stl'
# Radius for trimming_package; Effects the normal vector and rotation matrix
trimRadius = 250
# Radius for trimming the mesh around a point; Effects the histogram and slice
displayRadius = 250

normal_vector, lowest_point, rotation_matrix, rotation_matrix_point = help.trimming_package(
    filename, trimRadius)
print(f"The depth of the crater is {lowest_point}")
# Find the lowest point in the unrotated basis
lowest_point_r = rotation_matrix_point @ lowest_point

# Trim around the lowest_point
new_mesh = Mesh.from_file(
    filename)
help.trimCircleGivenPoint(new_mesh, lowest_point_r, displayRadius)
mesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
mesh.remove_infinite_values()
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
print(f"The mode of the heights is {baseline}")
plt.axvline(x=baseline, color='g', linestyle='dashed')
# plt.show()

# Getting a slice of the mesh in x-z plane
y, z = help.slicer(mesh, plane='x', buffer=2)
# plt.figure(2)
# plt.plot(y, z, 'bs', linewidth=2)
# plt.axhline(y=baseline, color='g', linestyle='dashed')
# plt.axis('equal')
# plt.show()

# Curvature
delta_z = np.zeros(len(z))
for i in range(len(z) - 1):
    change = (z[i+1] - z[i])
    delta_z[i] = change
mask1 = -2 > delta_z
mask2 = 2 < delta_z
mask = []
for i in range(len(mask1)):
    mask.append(mask1[i] or mask2[i])
delta_z = delta_z[mask]
y1 = y[mask]
for i in range(len(delta_z)):
    delta_z[i] = np.sum(delta_z[i:i+3])/3
plt.figure(3)
plt.plot(y1, delta_z, 'bs', linewidth=2)
# plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.axis('equal')


delta2_z = np.zeros(len(delta_z))
for i in range(len(delta_z) - 1):
    change = (delta_z[i+1] - delta_z[i])
    delta2_z[i] = change
plt.figure(4)
plt.plot(y1, delta2_z, 'bs', linewidth=2)
# plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.axis('equal')

max_change = np.where(delta2_z == np.max(delta2_z))
index = np.where(y == y1[max_change])
crater_start = z[index]
print(f"The start of the crater is at height = {crater_start}")


plt.figure(2)
plt.plot(y, z, 'bs', linewidth=2)
plt.axhline(y=crater_start, color='g', linestyle='dashed')
plt.axis('equal')

# print(help.zPoints(mesh))
help.move(mesh, np.array([0, 0, -30]))
# print(lowest_point)
# print(help.zPoints(mesh))
# Getting a slice of the mesh in x-z plane
x, y = help.slicer(mesh, plane='z', buffer=0.5)
plt.figure(5)
plt.plot(x, y, 'bs', linewidth=2)
# plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.axis('equal')
plt.show()
