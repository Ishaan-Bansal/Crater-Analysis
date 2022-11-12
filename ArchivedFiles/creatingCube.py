import imp
from stl import mesh
import math
import numpy
import trimesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create 3 faces of a cube
data = numpy.zeros(6, dtype=mesh.Mesh.dtype)
# Top of the cube
data['vectors'][0] = numpy.array([[0, 1, 1],
                                  [1, 0, 1],
                                  [0, 0, 1]])
data['vectors'][1] = numpy.array([[1, 0, 1],
                                  [0, 1, 1],
                                  [1, 1, 1]])
# Front face
data['vectors'][2] = numpy.array([[1, 0, 0],
                                  [1, 0, 1],
                                  [1, 1, 0]])
data['vectors'][3] = numpy.array([[1, 1, 1],
                                  [1, 0, 1],
                                  [1, 1, 0]])
# Left face
data['vectors'][4] = numpy.array([[0, 0, 0],
                                  [1, 0, 0],
                                  [1, 0, 1]])
data['vectors'][5] = numpy.array([[0, 0, 0],
                                  [0, 0, 1],
                                  [1, 0, 1]])
# Since the cube faces are from 0 to 1 we can move it to the middle by
# substracting .5
# data['vectors'] -= .5

print(data)

########################################
# Create a new mesh

your_mesh = mesh.Mesh(data, remove_empty_areas=False)

m = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
m.show()


########################################
# pyplot.show() method (does not work)

# # Create a new plot
# figure = pyplot.figure()
# axes = mplot3d.Axes3D(figure)
# # Render the cube
# axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))
# # Auto scale to the mesh size
# scale = your_mesh.points.flatten()
# axes.auto_scale_xyz(scale, scale, scale)
# # Show the plot to the screen
# pyplot.show()
