import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mesh = trimesh.load_mesh(
    r"C:\Users\ishaa\Desktop\Research Work\Crater_STL_Files\02_09_2022_250Torr_test3.STL")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(mesh.vertices[:, 0], mesh.vertices[:, 1],
                triangles=mesh.faces, Z=mesh.vertices[:, 2])
plt.show()
