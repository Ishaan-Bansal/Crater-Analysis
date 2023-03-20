import ray_tracing_curvature as rtc
import ray_tracing_X_Slices as rtx
import ray_tracing_Z_Slices as rtz
import os
import matplotlib.pyplot as plt
import trimesh
from stl import Mesh
import helper_functions as help
import numpy as np

# Write the relative path to the folder you want to load
path = "Test Crater Stls"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    if filename[-4:] != ".stl":
        continue
    print("Working on: " + filename)
    # mesh = help.load_mesh(filename)
    mesh = trimesh.load_mesh(filename)
    mesh.show()
    filename = filename[:-4]
    os.chdir(cwd)
    os.chdir("Crater_Properties")
    os.mkdir(filename)
    os.chdir(filename)

    depth, diameter, volume = rtc.crater_properties(mesh)

    file = open(filename + ".txt", 'w')
    info = f"Depth: {depth} , Diameter: {diameter} , Volume: {volume}"
    file.write(info)
    file.close()

    # if int(depth/10) > 0:
    #     for z in range(int(depth/10), int(depth), int(depth/10)):
    #         locations = rtz.z_slice(mesh, z)
    #         figz = plt.figure()
    #         x, y = locations[:, 0], locations[:, 1]
    #         plt.plot(x, y, 'k.', linewidth=2)
    #         plt.axis('equal')

    #         figz.savefig(filename + "_Z-Slice_" + str(z) + ".png")
    #         plt.close()
    #         np.savetxt(filename + "_Z-Slice_" + str(z) + ".csv", locations, delimeter=',')

    # locations = rtz.z_slice(mesh, depth)
    # figz = plt.figure()
    # x, y = locations[:, 0], locations[:, 1]
    # plt.plot(x, y, 'k.', linewidth=2)
    # plt.axis('equal')

    # figz.savefig(filename + "_Z-Slice_" + str(depth) + ".png")    
    # plt.close()
    # np.savetxt(filename + "_Z-Slice_" + str(depth) + ".csv", locations, delimeter=',')
    
    # locations2 = rtx.x_slice(mesh)
    # figx = plt.figure()
    # x, z = locations2[:, 0], locations2[:, 2]
    # plt.plot(x, z, 'k.', linewidth=2)
    # plt.axis('equal')

    # figx.savefig(filename + "_X-Slice" + ".png")    
    # plt.close()
    # np.savetxt(filename + "_X-Slice" + ".csv", locations2, delimeter=',')


    os.chdir(cwd)
    os.chdir(path)

print("Finished")