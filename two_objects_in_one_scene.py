import trimesh
# from trimesh.voxel import creation
from trimesh.proximity import signed_distance
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
# import pyacvd
# import pyvista
from skimage import measure
import skimage
import scipy
from stl import Mesh
import helper_functions as help

# mesh = trimesh.load_mesh(
#     r"C:\Users\ishaa\Desktop\Research Work\Crater_STL_Files\02_09_2022_250Torr_test3.STL")

# Trim the crater mesh using numpy stl
your_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
t = help.trimCircle(your_mesh, 250)
trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(t.vectors))

sphere = trimesh.load_mesh(
    r'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/sphere.stl')

# sphere = trimesh.icosphere(3, 1, (3,))

# Change sphere color
sphere.visual.face_colors = trimesh.visual.random_color()

# PLane
plane = help.plane()
planeTrimesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(plane.vectors))
planeTrimesh.visual.face_colors = trimesh.visual.random_color()

# Adding two meshs to the same scene
scene = trimesh.Scene()
scene.add_geometry(sphere)
scene.add_geometry(planeTrimesh)
scene.add_geometry(trimmed)
scene.show()

# # Matplotlib version two
# figure1 = plt.figure()
# axes1 = figure1.add_subplot(111, projection='3d')

# axes1.add_collection3d(mplot3d.art3d.Poly3DCollection(trimmed.triangles,
#                                                       facecolor=np.array(
#                                                           [(1., 1., 1., 1)]),
#                                                       edgecolor=np.array([(0., 0., 0., 1)])))

# axes1.plot_trisurf(trimmed.vertices[:, 0], trimmed.vertices[:, 1],
#                    triangles=trimmed.faces, Z=trimmed.vertices[:, 2])
# plt.show()
