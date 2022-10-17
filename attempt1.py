import trimesh
from stl import Mesh
import helper_functions as help
import numpy as np

your_mesh = Mesh.from_file(
    'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')
t = help.trimCircle(your_mesh, 250)
trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(t.vectors))

# Rotate the mesh 90 degrees
trimmed.apply_transform(
    trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))

normalVec = help.unitVectorsToZplane(trimmed.face_normals)
print(normalVec/np.linalg.norm(normalVec))

normalPLane = help.vectorToPlane(normalVec)
normalPLane.visual.face_colors = trimesh.visual.random_color()


# zerozeroPlane = help.plane()
# zerozeroPlane

# Rotate the mesh so that the normal vector faces only in the z-direction
# trimmed.apply_transform(
#     trimesh.transformations.rotation_matrix(
#         np.pi/2-help.rotationAngle(normalVec, [1, 0, 0]), [1, 0, 0]))
# trimmed.apply_transform(
#     trimesh.transformations.rotation_matrix(
#         np.pi/2-help.rotationAngle(normalVec, [0, 1, 0]), [0, 1, 0]))


# Lowest point
lowest_point = help.lowest_point(trimmed, 2)
print(lowest_point)
# Trim around the lowest_point
# newMesh = help.trimCircleGivenPoint(t, lowest_point, 200)
# newTrim = trimesh.Trimesh(**trimesh.triangles.to_kwargs(newMesh.vectors))

scene = trimesh.Scene()
scene.add_geometry(normalPLane)
scene.add_geometry(trimmed)
# scene.add_geometry(zerozeroPlane)
scene.show()
# help.lowest_point(trimmed, 2)
