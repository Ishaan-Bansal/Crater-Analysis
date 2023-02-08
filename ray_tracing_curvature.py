import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats
import scipy.ndimage
import ray_tracing_Z_Slices as rtz
import ray_tracing_X_Slices as rtx

# Write the relative path to the file you want to load
filename = 'Crater_STL_Files/02_09_2022_250Torr_test3.stl'
# filename = "Crater_STL_Files/2022_11_01_50mTorr_h10_1s_032gs_noacrylic.stl"
# Radius for trimming_package; Effects the normal vector and rotation matrix
trimRadius = 40
# Radius for trimming the mesh around a point; Effects the histogram and slice
displayRadius = 40
# The bounds of the square in which the rays are limited
bounds = 20
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
mesh.show()

# ray_directions_array = np.array([])
# ray_origins_array = np.array([])
# index_tri_array = np.array([])
locations_array = np.array([])
x_locations = []
y_locations = []
z_locations = []

divergence_vis = []

for y in range(-bounds, bounds+spacing, spacing):
    locations_row = []
    for x in range(-bounds, bounds+spacing, spacing):
        ray_directions = np.array([[0, 0, -1]])
        ray_origins = np.array([[x, y, 100]])
        divergence_vis.append(np.array([0, 0, 0]))
        locations, index_ray, index_tri = mesh.ray.intersects_location(
            ray_origins=ray_origins, ray_directions=ray_directions)
        locations_row.append(locations[0])
        # ray_directions_array = np.append(ray_directions_array, ray_directions)
        # ray_origins_array = np.append(ray_origins_array, ray_origins)
        # index_tri_array = np.append(index_tri_array, index_tri)

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

# fig = plt.figure(1)
# ax = plt.axes()
# CS = ax.contour(x_locations, y_locations, double_gradient_norm, cmap='turbo')
# ax.clabel(CS, inline=True, fontsize=10)
# ax.set_title('2D Gradient Contour')

# fig = plt.figure(2)
# ax = plt.axes(projection='3d')
# ax.contour3D(x_locations, y_locations, double_gradient_norm, 50, cmap='turbo')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('gradient')
# ax.set_title('3D Gradient Contour')

# fig = plt.figure(3)
# img = plt.imshow(z_locations, interpolation='none', cmap='turbo')

fig = plt.figure(4)
img = plt.imshow(double_gradient_norm, interpolation='none', cmap='turbo')
# plt.show()

ridge_indices = []
ridge_z = []
radius_x = int(z_locations.shape[0]/2)
radius_y = int(z_locations.shape[1]/2)
dgn = []

for theta_deg in range(0, 360):
    # -- Extract the line...
    # Make a line with "num" points...
    theta_rad = np.radians(theta_deg)
    # These are in _pixel_ coordinates!!
    x1, y1 = radius_x - radius_x * \
        np.cos(theta_rad), radius_y - radius_y*np.sin(theta_rad)
    x0, y0 = radius_x, radius_y
    length = int(np.round(np.hypot(x1-x0, y1-y0)))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    mask = np.sqrt((x-radius_x)**2 + (y-radius_y) **
                   2) >= (0.1*(radius_x+radius_y)/2)
    x, y = x[mask], y[mask]

    # Extract the values along the line
    x = np.round(x).astype(int)
    y = np.round(y).astype(int)
    z_line = double_gradient_norm[x, y]

#     # print("X:")
#     # print(x)
    index = np.argmax(z_line)
    pair = np.array([y[index], x[index]])
    ridge_indices.append(pair)
    height_line = z_locations[x, y]
    ridge_z.append(height_line[index])

#     # plt.plot(x, y, 'ko')


ridge_indices = np.array(ridge_indices)
ridge_z = np.array(ridge_z)
# print(ridge_indices)
# plt.figure(5)
# plt.plot(ridge_indices[:, 0], ridge_indices[:, 1], 'ko')
# plt.show()

# size = ridge_indices.shape[0]
# ridge_indices_copy = np.ones((size, 3))
# ridge_indices_copy[:, 0] = ridge_indices[:, 0]
# ridge_indices_copy[:, 1] = ridge_indices[:, 1]
# plane_eq = np.linalg.lstsq(ridge_indices_copy, ridge_z)

# a, b, d = plane_eq[0]
# plane = help.create_trimesh_plane(a, b, d)

crater_start = np.mean(ridge_z)
# a, b, d = 0, 0, crater_start
# plane = help.create_trimesh_plane(a, b, d)

# scene = trimesh.Scene([mesh,
#                        plane])
# scene.show()

x, y = [], []
for i in ridge_indices:
    x.append(x_locations[i[0]][i[1]])
    y.append(y_locations[i[0]][i[1]])
x, y = np.array(x), np.array(y)

b = x**2 + y**2
M_a = np.ones((x.size, 3))
M_a[:, 0] = x
M_a[:, 1] = y
sol = np.linalg.lstsq(M_a, b, rcond=None)

x_c = sol[0][0]/2
y_c = sol[0][1]//2
radius = np.sqrt(sol[0][2] + x_c**2 + y_c**2)

# figure, axes = plt.subplots()
# Drawing_uncolored_circle = plt.Circle((x_c, y_c), radius, fill=False)
# plot = plt.plot(x, y, 'k.', linewidth=2)

# axes.set_xlim(-300, 300)
# axes.set_ylim(-300, 300)
# axes.add_patch(Drawing_uncolored_circle)
# plt.title('Circle')
# plt.show()
