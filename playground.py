import trimesh
from stl import Mesh
import helper_functions as help
import numpy as np

your_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
help.trimCircle(your_mesh, 250)
trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))

# Rotate the mesh 90 degrees
trimmed.apply_transform(
    trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))

normalVec = help.unitVectorsToZplane(trimmed.face_normals)
print(normalVec)

normalPLane = help.vectorToPlane(normalVec)
normalPLane.visual.face_colors = trimesh.visual.random_color()

trimmed.apply_transform(help.rotationMatrix(normalVec, np.array([0, 0, 1])))


zerozeroPlane = help.plane()
zerozeroPlane.visual.face_colors = trimesh.visual.random_color()

# Lowest point
lowest_point = help.lowest_point(trimmed, 2)
print(lowest_point)

# Trim around the lowest_point
# new_mesh = Mesh.from_file(
#     'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
# # help.trimCircleGivenPoint(new_mesh, lowest_point, 200)
# newTrim = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
# newTrim.apply_transform(
#     trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
# newTrim.apply_transform(help.rotationMatrix(normalVec, np.array([0, 0, 1])))

sphere = trimesh.creation.uv_sphere()
help.move(sphere, lowest_point)
sphere.visual.face_colors = trimesh.visual.random_color()


slice = trimmed.section(plane_origin=lowest_point,
                        plane_normal=[0, 0, 1])
# plane = help.vectorToPlane([0, 0, 1])
# help.move(plane, lowest_point)

scene = trimesh.Scene()
scene.add_geometry(normalPLane)
scene.add_geometry(trimmed)
scene.add_geometry(sphere)
# scene.add_geometry(plane)
scene.show()

slice.show()

# hello world
