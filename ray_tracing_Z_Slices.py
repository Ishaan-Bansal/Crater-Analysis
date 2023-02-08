import helper_functions as help
import trimesh
import numpy as np
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats
import os


def z_slice(mesh, depth):
    
    ray_origins = np.ones((360, 3))*[0, 0, depth]
    ray_directions = np.zeros((360, 3))
    for i in range(0, 360):
        theta = np.radians(i)
        ray_directions[i] = np.array([np.cos(theta), np.sin(theta), 0])

    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions)

    # ray_visualize = trimesh.load_path(np.hstack((ray_origins,
    #                                             ray_origins + ray_directions*5.0)).reshape(-1, 2, 3))
    # mesh.visual.face_colors = [255, 255, 255, 255]
    # mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
    # scene = trimesh.Scene([mesh,
    #                        ray_visualize])
    # scene.show()
    return locations


def z_slice_load_file(filename, depth):
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

    ray_origins = np.ones((360, 3))*[0, 0, depth]
    ray_directions = np.zeros((360, 3))
    for i in range(0, 360):
        theta = np.radians(i)
        ray_directions[i] = np.array([np.cos(theta), np.sin(theta), 0])

    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions)

    # ray_visualize = trimesh.load_path(np.hstack((ray_origins,
    #                                             ray_origins + ray_directions*5.0)).reshape(-1, 2, 3))
    # mesh.visual.face_colors = [255, 255, 255, 255]
    # mesh.visual.face_colors[index_tri] = [255, 0, 0, 255]
    # scene = trimesh.Scene([mesh,
    #                        ray_visualize])
    # scene.show()
    return locations


# locations = z_slice("Crater_STL_Files/02_09_2022_6Torr_test2.stl", 20)
# x, y = locations[:, 0], locations[:, 1]
# plt.plot(x, y, 'k.', linewidth=2)
# plt.axis('equal')
# plt.show()

# # Write the relative path to the folder you want to load
# path = "Crater_STL_Files"
# cwd = os.getcwd()
# os.chdir(path)
# for filename in os.listdir():
#     if filename[-4:] != ".stl":
#         continue
#     locations = z_slice(filename, 20)
#     x, y = locations[:, 0], locations[:, 1]
#     plt.plot(x, y, 'k.', linewidth=2)
#     plt.axis('equal')
#     # plt.show()
#     filename = filename[:-4]
#     os.chdir(cwd)
#     plt.savefig("Crater_Slices/" + filename + "_Z-Slice" + ".png")
#     os.chdir(path)
#     plt.close()
