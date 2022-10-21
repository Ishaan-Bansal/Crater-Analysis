from sympy import centroid
import trimesh
from stl import Mesh
import helper_functions as help
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f, norm
import matplotlib.mlab as mlab
import statistics as stats

your_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
help.trimCircle(your_mesh, 250)
trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))

# Rotate the mesh 90 degrees
trimmed.apply_transform(
    trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))

normalVec = help.unitVectorsToZplane(trimmed.face_normals)
# print(normalVec)

# normalPLane = help.vectorToPlane(normalVec)
# normalPLane.visual.face_colors = trimesh.visual.random_color()

trimmed.apply_transform(help.rotationMatrix(normalVec, np.array([0, 0, 1])))


# Lowest point
lowest_point = help.lowest_point(trimmed, 2)

# sphere = trimesh.creation.uv_sphere()
# help.move(sphere, lowest_point)
# sphere.visual.face_colors = trimesh.visual.random_color()

# Histogram of heights
heights = help.zPoints(trimmed)
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
# slice = trimmed.section(plane_origin=trimmed.centroid,
#                         plane_normal=[1, 0, 0])
# slice.show()
y, z = help.slicer(trimmed)
plt.figure(2)
plt.plot(y, z, 'bs', linewidth=2)
plt.axhline(y=baseline, color='g', linestyle='dashed')
plt.show()

# scene = trimesh.Scene()
# scene.add_geometry(normalPLane)
# scene.add_geometry(trimmed)
# scene.add_geometry(sphere)
# # scene.add_geometry(plane)
# scene.show()
