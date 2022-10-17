from stl import Mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import trimesh
import helper_functions as help
import numpy as np
import statistics as stats

# Load the STL files
your_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')

m = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
# print(m.centroid)
# m.show()

t = help.trimCircle(your_mesh, 250)
trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(t.vectors))
trimmed.show()

# Initialize the array of vertices from the mesh
points = []
data = your_mesh.vectors
for triangle in data:
    points.append(help.triangleToPoint(triangle))
# print(data)

min = points[0][2]
z_points = []

# Add all the z components to the empty array
for point in points:
    z_points.append(point[2])
    if point[2] < min:
        minPoint = point
        min = point[2]

print(f"Lowest point: {minPoint}")
base = stats.median(z_points)
print(f"Baseline: {base}")

depth = - min

print(depth)
# Initialize the array of unit normal vectors from the mesh using a helper function
# vectors = help.mesh_to_vectors(t)
vectors2 = your_mesh.units
vectors3 = your_mesh.normals

# print(len(vectors2))
# print(len(vectors3))

# if np.array_equal(np.sort(vectors.flat), np.sort(vectors2.flat)):
#     print("vectors = vectors2")
# if np.array_equal(np.sort(vectors.flat), np.sort(vectors3.flat)):
#     print("vectors = vectors3")

# for i in range(len(vectors3)):
#     vectors3[i] = vectors3[i]/np.linalg.norm(vectors3[i])

# if np.array_equal(np.sort(vectors3.flat), np.sort(vectors2.flat)):
#     print("vectors3 = vectors2")

# Calculate the normal vector of the Z plane

normalVec = help.unitVectorsToZplane(vectors2)
print(f"normalVec = {normalVec}")
# slice = trimmed.section(plane_origin=trimmed.centroid, plane_normal=normalVec)
# slice.show()
