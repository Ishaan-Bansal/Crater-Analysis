import ray_tracing_curvature as rtc
import ray_tracing_X_Slices as rtx
import ray_tracing_Z_Slices as rtz
import os
import matplotlib.pyplot as plt
import trimesh
from stl import Mesh
import helper_functions as help
import numpy as np
import pandas as pd

df = pd.DataFrame()

# Write the relative path to the folder you want to load
path = "Crater_STL_Files"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    if filename[-4:] != ".stl":
        continue
    print("Working on: " + filename)

    mesh = help.load_mesh(filename)

    # angle = np.pi/2*0
    # direction = [1, 0, 0]
    # center = [0, 0, 0]

    # rot_matrix = trimesh.transformations.rotation_matrix(angle, direction, center)

    # mesh.apply_transform(rot_matrix)
    
    filename = filename[:-4]

    # radius_row, depth_col = help.from_filename(filename)

    # if depth_col not in df.columns:
    #     df[depth_col] = pd.Series(dtype=object)
    # if not df.index.isin([radius_row]).any():
    #     df.loc[radius_row] = [str()] * len(df.columns) 

    os.chdir(cwd)
    os.chdir("Crater_Properties_2")
    os.mkdir(filename)
    os.chdir(filename)

    depth, diameter, volume = rtc.crater_properties(mesh)

    file = open(filename + ".txt", 'w')
    info = f"Depth: {depth} , Diameter: {diameter} , Volume: {volume}"
    file.write(info)
    file.close()

    # df.loc[radius_row, depth_col] = info

    if int(depth/10) > 0:
        for z in range(int(depth/10), int(depth), int(depth/10)):
            locations = rtz.z_slice(mesh, z)
            figz = plt.figure()
            x, y = locations[:, 0], locations[:, 1]
            plt.plot(x, y, 'k.', linewidth=2)
            plt.axis('equal')

            figz.savefig(filename + "_Z-Slice_" + str(z) + ".png")
            plt.close()
            np.savetxt(filename + "_Z-Slice_" + str(z) + ".csv", locations, delimiter=',')

    locations = rtz.z_slice(mesh, depth)
    figz = plt.figure()
    x, y = locations[:, 0], locations[:, 1]
    plt.plot(x, y, 'k.', linewidth=2)
    plt.axis('equal')

    figz.savefig(filename + "_Z-Slice_" + str(depth) + ".png")    
    plt.close()
    np.savetxt(filename + "_Z-Slice_" + str(depth) + ".csv", locations, delimiter=',')
    
    locations2 = rtx.x_slice(mesh)
    figx = plt.figure()
    x, z = locations2[:, 0], locations2[:, 2]
    plt.plot(x, z, 'k.', linewidth=2)
    plt.axis('equal')

    figx.savefig(filename + "_X-Slice" + ".png")    
    plt.close()
    np.savetxt(filename + "_X-Slice" + ".csv", locations2, delimiter=',')


    os.chdir(cwd)
    os.chdir(path)

os.chdir(cwd)
df.to_csv("Crater_Properties/test_craters.csv")

print("Finished")