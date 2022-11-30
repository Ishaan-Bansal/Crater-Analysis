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
# Find the lowest point in the unrotated basis
lowest_point_r = rotation_matrix_point @ lowest_point

# Trim around the lowest_point
new_mesh = Mesh.from_file(
    filename)
help.trimCircleGivenPoint(new_mesh, lowest_point_r, displayRadius)
mesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
mesh.remove_infinite_values()
mesh.apply_transform(rotation_matrix)

help.move(mesh, lowest_point)

ray_directions_array = np.array([])
ray_origins_array = np.array([])
index_tri_array = np.array([])
locations_array = []

for y in range(-150, 150, 10):
    ray_directions = np.ones((30, 3))*[0, 0, -1]
    ray_origins = np.zeros((30, 3))
    for i in range(-150, 150, 10):
        ray_origins[int(i/10)] = np.array([i, y, 100])
    ray_directions_array = np.vstack(
        [ray_directions_array, ray_directions]) if ray_directions_array.size else ray_directions
    # ray_directions_array = np.append(ray_directions_array, ray_directions)
    ray_origins_array = np.vstack(
        [ray_origins_array, ray_origins]) if ray_origins_array.size else ray_origins
    # ray_origins_array = np.append(ray_origins_array, ray_origins)
    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins, ray_directions=ray_directions)
    index_tri_array = np.concatenate((index_tri_array, index_tri), axis=0)
    locations_array.append([locations])

index_tri_array = index_tri_array.astype(int)
locations_array = np.array(locations_array, dtype=object)
# print(locations_array)

ray_visualize = trimesh.load_path(np.hstack((ray_origins_array,
                                             ray_origins_array + ray_directions_array*5.0)).reshape(-1, 2, 3))
mesh.visual.face_colors = [255, 255, 255, 255]
mesh.visual.face_colors[index_tri_array] = [255, 0, 0, 255]
scene = trimesh.Scene([mesh,
                       ray_visualize])
scene.show()

print(locations_array[0])

dz_dy = np.ones((30, 30, 3))*[0, 0, 0]
for i in range(30):
    for j in range(30):
        if j == 0 or j == 29:
            continue
        if locations_array[i].shape[0] == 0:
            continue
        dz_dy[i][j] = np.array(locations_array[i][0]
                               [j+1] - locations_array[i][0][j-1])

dz_dx = np.zeros((30, 30, 3))*[0, 0, 0]
for i in range(30):
    for j in range(30):
        if i == 0 or i == 29:
            continue
        if locations_array[i].shape[0] == 0:
            continue
        dz_dx[i][j] = np.array(locations_array[i+1][0]
                               [j] - locations_array[i-1][0][j])

dz = dz_dy+dz_dx

curvature_visualize = trimesh.load_path(np.hstack((locations_array,
                                                   locations_array + dz*0.00000000001)).reshape(-1, 2, 3))
scene = trimesh.Scene([mesh,
                       curvature_visualize])
scene.show()
