import trimesh
import numpy as np
import helper_functions as help
from stl import Mesh

################################################################

# Heat map for crater visualization

# Write the relative path to the file you want to load
filename = 'Crater_STL_Files/02_09_2022_6Torr_test2.stl'
# Radius for trimming_package; Effects the normal vector and rotation matrix
trimRadius = 250
# Radius for trimming the mesh around a point; Effects the histogram and slice
displayRadius = 250

normal_vector, lowest_point, rotation_matrix, rotation_matrix_point = help.trimming_package(
    filename, trimRadius)

# Load the mesh
your_mesh = Mesh.from_file(filename)
# Trim the mesh for display
help.trimCircleGivenPoint(your_mesh, lowest_point, displayRadius)
# Convert to Trimesh object
crater = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
# Rotate the the mesh to the x-y plane
crater.apply_transform(rotation_matrix)

zpoints = help.zPoints(crater)  # - lowest_point[2] + 1
log_zpoints = np.log10(np.abs(zpoints))
crater.visual.vertex_colors = trimesh.visual.interpolate(
    log_zpoints, color_map='turbo')

crater.show()
