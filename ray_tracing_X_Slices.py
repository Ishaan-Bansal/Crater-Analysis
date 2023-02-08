import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats
import os
# import glob


def x_slice_load_file(filename):
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

    ray_directions = np.ones((500, 3))*[0, 0, -1]
    ray_origins = np.zeros((500, 3))
    for i in range(-250, 250):
        ray_origins[i] = np.array([i, 0, 100])

    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions)

    ray_visualize = trimesh.load_path(np.hstack((ray_origins,
                                                 ray_origins + ray_directions*5.0)).reshape(-1, 2, 3))
    mesh.visual.face_colors = [255, 255, 255, 255]
    mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
    scene = trimesh.Scene([mesh,
                           ray_visualize])
    scene.show()
    return locations


def x_slice(mesh):
    ray_directions = np.ones((500, 3))*[0, 0, -1]
    ray_origins = np.zeros((500, 3))
    for i in range(-250, 250):
        ray_origins[i] = np.array([i, 0, 100])

    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions)

    ray_visualize = trimesh.load_path(np.hstack((ray_origins,
                                                 ray_origins + ray_directions*5.0)).reshape(-1, 2, 3))
    mesh.visual.face_colors = [255, 255, 255, 255]
    mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
    scene = trimesh.Scene([mesh,
                           ray_visualize])
    scene.show()
    return locations


# Write the relative path to the folder you want to load
# path = "Crater_STL_Files"
# cwd = os.getcwd()
# os.chdir(path)
# for filename in os.listdir():
#     if filename[-4:] != ".stl":
#         continue

#     locations = x_slice(filename)
#     x, z = locations[:, 0], locations[:, 2]
#     plt.plot(x, z, 'k.', linewidth=2)
#     plt.axis('equal')
#     plt.title(filename)
#     # plt.show()
#     filename = filename[:-4]
#     os.chdir(cwd)
#     plt.savefig("Crater_Slices/" + filename + "_x-Slice" + ".png")
#     os.chdir(path)
#     plt.close()


# x_slice("Crater_STL_Files/02_09_2022_6Torr_test2.stl")
