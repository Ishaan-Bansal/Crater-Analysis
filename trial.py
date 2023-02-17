import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh

filename = "Crater_STL_Files/02_09_2022_250Torr_test3.stl"
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

# The bounds of the square in which the rays are limited
bounds = 10
# The spacing between each ray
spacing = 5

ray_directions_array = np.array([])
ray_origins_array = np.array([])
index_tri_array = np.array([], dtype=int)
locations_array = np.array([])
x_locations = []
y_locations = []
z_locations = []

for y in range(-bounds, bounds+spacing, spacing):
    locations_row = []
    for x in range(-bounds, bounds+spacing, spacing):
        ray_directions = np.array([[0, 0, -10]])
        ray_origins = np.array([[x, y, 100]])
        locations, index_ray, index_tri = mesh.ray.intersects_location(
            ray_origins=ray_origins, ray_directions=ray_directions)
        if locations.any():
            locations_row.append(locations[0])
            ray_directions_array = np.append(
                ray_directions_array, ray_directions)
            ray_origins_array = np.append(ray_origins_array, ray_origins)
            index_tri_array = np.append(index_tri_array, index_tri)
        else:
            locations_row.append([0, 0, 0])
            ray_directions_array = np.append(
                ray_directions_array, [0, 0, 0])
            ray_origins_array = np.append(ray_origins_array, [0, 0, 0])
            index_tri_array = np.append(index_tri_array, [0, 0, 0])

    locations_row = np.array(locations_row)
    x_locations.append(locations_row[:, 0])
    y_locations.append(locations_row[:, 1])
    z_locations.append(locations_row[:, 2])


x_locations = np.array(x_locations, dtype=float)
y_locations = np.array(y_locations, dtype=float)
z_locations = np.array(z_locations, dtype=float)

ray_visualize = trimesh.load_path(np.hstack((ray_origins_array,
                                             ray_origins_array + ray_directions_array*5.0)).reshape(-1, 2, 3))
# ray_visualize = trimesh.load_path(
#     np.hstack((np.array([10, 10, 100]), np.array([10, 10, 50]))).reshape(-1, 2, 3))
# mesh.visual.face_colors = [255, 255, 255, 255]
# mesh.visual.face_colors[index_tri_array] = [255, 0, 0, 255]
scene = trimesh.Scene([mesh,
                       ray_visualize])
scene.show()

print(z_locations)

vol = help.pixels_within_circle(z_locations, 300, spacing, 0)
print(vol)
