from vedo import trimesh2vedo, show
import trimesh
import numpy as np
import helper_functions as help

################################################################

# Heat map for crater visualization

crater = trimesh.load_mesh(
    r'C:/Users/ishaa/Desktop/Research Work/Crater_STL_Files/02_09_2022_6Torr_test2.stl')

normalVec = np.array(
    [0.03792327795568258, 0.019775532970647217, 0.4442363538833091])
crater.apply_transform(
    trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
crater.apply_transform(help.rotationMatrix(normalVec, np.array([0, 0, 1])))

lowest_point = np.array([4.11343818,  12.56165876, - 47.95400997])


def lowest_vector(vertex):
    vertex[:, 2] -= 47.95400997
    return vertex


print((lowest_vector(crater.vertices)))

radii = np.linalg.norm(lowest_vector(crater.vertices), axis=1)
crater.visual.vertex_colors = trimesh.visual.interpolate(
    radii, color_map='viridis')

crater.show()
