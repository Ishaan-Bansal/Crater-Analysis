import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats
import scipy.ndimage

# Write the relative path to the file you want to load
filename = 'Crater_STL_Files/02_09_2022_250Torr_test3.stl'
# filename = "Crater_STL_Files/2022_11_01_50mTorr_h10_1s_032gs_noacrylic.stl"
# Radius for trimming_package; Effects the normal vector and rotation matrix
trimRadius = 250
# Radius for trimming the mesh around a point; Effects the histogram and slice
displayRadius = 250
# The bounds of the square in which the rays are limited
bounds = 150
# The spacing between each ray
spacing = 5
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
# mesh.show()

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

double_gradient = np.gradient(gradient_norm)
double_gradient_norm = np.sqrt(double_gradient[0]**2 + double_gradient[1]**2)

print("X Locations:")
print(x_locations)
print("Y Locations:")
print(y_locations)


fig = plt.figure(1)
ax = plt.axes()
CS = ax.contour(x_locations, y_locations, double_gradient_norm, cmap='turbo')
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('2D Gradient Contour')

fig = plt.figure(2)
ax = plt.axes(projection='3d')
ax.contour3D(x_locations, y_locations, double_gradient_norm, 50, cmap='turbo')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('gradient')
ax.set_title('3D Gradient Contour')

fig = plt.figure(3)
img = plt.imshow(z_locations, interpolation='none', cmap='turbo')

fig = plt.figure(4)
img = plt.imshow(double_gradient_norm, interpolation='none', cmap='turbo')
# plt.show()

ridge_indices = []
for theta_deg in range(0, 360, 10):
    # -- Extract the line...
    # Make a line with "num" points...
    theta_rad = np.radians(theta_deg)
    # These are in _pixel_ coordinates!!
    radius_x = int(z_locations.shape[0]/2)
    radius_y = int(z_locations.shape[1]/2)
    x1, y1 = radius_x - radius_x * \
        np.cos(theta_rad), radius_y - radius_y*np.sin(theta_rad)
    x0, y0 = radius_x, radius_y
    length = int(np.round(np.hypot(x1-x0, y1-y0)))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    mask = (x-radius_x)**2 + (y-radius_y)**2 >= 0.2*(radius_x+radius_y)/2
    x, y = x[mask], y[mask]

    # Extract the values along the line
    x = np.floor(x).astype(int)
    y = np.floor(y).astype(int)
    z_line = double_gradient_norm[x, y]
    # x, y = x + z_locations.shape[0], y + z_locations.shape[1]

    # print("X:")
    # print(x)
    index = np.argmax(z_line)
    pair = np.array([x[index], y[index]])
    ridge_indices.append(pair)


ridge_indices = np.array(ridge_indices)
print(ridge_indices)
ridge_x = ridge_indices[:, 0]
ridge_y = ridge_indices[:, 1]
# plt.figure(5)
plt.plot(ridge_x, ridge_y, 'ko')
plt.show()

# radius_x = int(z_locations.shape[0]/2)
# radius_y = int(z_locations.shape[1]/2)
# x0, y0 = radius_x - radius_x * \
#     np.cos(0), radius_y - radius_y*np.sin(0)
# x1, y1 = radius_x, radius_y
# length = int(np.round(np.hypot(x1-x0, y1-y0)))
# x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
# mask = (x-radius_x)**2 + (y-radius_y)**2 >= 0.2*(radius_x+radius_y)/2
# x, y = x[mask], y[mask]

# # Extract the values along the line
# x = x.round().astype(int)
# y = y.round().astype(int)
# z_line = double_gradient_norm[x, y]
# print(z_line.size)
# print(x.size)


# plt.figure()
# plt.plot(np.sqrt((x-radius_x)**2 + (y-radius_y)**2), z_line)
# plt.show()
