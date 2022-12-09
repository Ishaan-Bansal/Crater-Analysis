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
# The bounds of the square in which the rays are limited
bounds = 150
# The spacing between each ray
spacing = 10
# If the spacing is great, than the resolution of the data is high and vice versa

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
locations_array = np.array([])
x_locations = []
y_locations = []
z_locations = []

divergence_vis = []

for y in range(-bounds, bounds+spacing, spacing):
    # ray_directions = np.ones((10, 3))*[0, 0, -1]
    # ray_origins = np.zeros((10, 3))
    locations_row = []
    for x in range(-bounds, bounds+spacing, spacing):
        ray_directions = np.array([[0, 0, -1]])
        ray_origins = np.array([[x, y, 100]])
        divergence_vis.append(np.array([0, 0, 0]))
        locations, index_ray, index_tri = mesh.ray.intersects_location(
            ray_origins=ray_origins, ray_directions=ray_directions)
        locations_row.append(locations[0])
    locations_row = np.array(locations_row)
    x_locations.append(locations_row[:, 0])
    y_locations.append(locations_row[:, 1])
    z_locations.append(locations_row[:, 2])


# index_tri_array = index_tri_array.astype(int)
x_locations = np.array(x_locations, dtype=float)
y_locations = np.array(y_locations, dtype=float)
z_locations = np.array(z_locations, dtype=float)
divergence_vis = np.array(divergence_vis)

# ray_visualize = trimesh.load_path(np.hstack((ray_origins_array,
#                                              ray_origins_array + ray_directions_array*5.0)).reshape(-1, 2, 3))
# mesh.visual.face_colors = [255, 255, 255, 255]
# mesh.visual.face_colors[index_tri_array] = [255, 0, 0, 255]
# scene = trimesh.Scene([mesh,
#                        ray_visualize])
# scene.show()

gradient = np.gradient(z_locations)
gradient_norm = np.sqrt(gradient[0]**2 + gradient[1]**2)

print("X Locations:")
print(x_locations)
print("Y Locations:")
print(y_locations)


fig = plt.figure(1)
ax = plt.axes()
CS = ax.contour(x_locations, y_locations, gradient_norm, cmap='turbo')
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('2D Gradient Contour')
plt.show()

fig = plt.figure(2)
ax = plt.axes(projection='3d')
ax.contour3D(x_locations, y_locations, gradient_norm, 50, cmap='turbo')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('gradient')
ax.set_title('3D Gradient Contour')
plt.show()

fig = plt.figure(3)
img = plt.imshow(z_locations, interpolation='none', cmap='turbo')
# img.setTitle("Height Heat Image")
plt.show()

fig = plt.figure(4)
img = plt.imshow(gradient_norm, interpolation='none', cmap='turbo')
# img.setTitle("Gradient Heat Image")
plt.show()
