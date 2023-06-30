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

columns = ['ID', 'Depth', 'Diameter', 'Volume', "Radius of curvature", "Radius of Sphere", "Radius of Sphere Fit", "Standard deviation"]
df = pd.DataFrame(columns=columns)

# Write the relative path to the folder you want to load
path = "Testing/Foam Sphere STLS Post Proccessed/All"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    if filename[-4:] != ".stl":
        continue
    print("Working on: " + filename)

    mesh = help.load_mesh(filename, trimRadius=100, displayRadius=100)
    
    filename = filename[:-4]
    
    os.chdir(cwd)
    os.chdir("Testing/Foam Sphere Analysis/All - mean")
    os.mkdir(filename)
    os.chdir(filename)

    depth, diameter, volume, ridge_height = rtc.crater_properties(mesh, bounds=70)

    # if int(depth/10) > 0:
    #     for z in range(int(depth/10), int(depth), int(depth/10)):
    #         locations = rtz.z_slice(mesh, z)
    #         figz = plt.figure()
    #         x, y = locations[:, 0], locations[:, 1]
    #         plt.plot(x, y, 'k.', linewidth=2)
    #         plt.axis('equal')

    #         figz.savefig(filename + "_Z-Slice_" + str(z) + ".png")
    #         plt.close()
    #         np.savetxt(filename + "_Z-Slice_" + str(z) + ".csv", locations, delimiter=',')

    # locations = rtz.z_slice(mesh, depth)
    # figz = plt.figure()
    # x, y = locations[:, 0], locations[:, 1]
    # plt.plot(x, y, 'k.', linewidth=2)
    # plt.axis('equal')

    # figz.savefig(filename + "_Z-Slice_" + str(depth) + ".png")    
    # plt.close()
    # np.savetxt(filename + "_Z-Slice_" + str(depth) + ".csv", locations, delimiter=',')
    
    locations2 = rtx.x_slice(mesh)
    # figx = plt.figure()
    # x, z = locations2[:, 0], locations2[:, 2]
    # plt.plot(x, z, 'k.', linewidth=2)
    # plt.axis('equal')

    # figx.savefig(filename + "_X-Slice" + ".png")    
    # plt.close()
    # np.savetxt(filename + "_X-Slice" + ".csv", locations2, delimiter=',')

    # Foam testing stuff
    roc = help.radius_of_curvature(locations2, depth=depth)

    ros = help.radius_of_sphere(depth=depth, diameter=diameter)

    roS, cx, cy, cz, err = help.sphereFit(mesh, radius=diameter/2, depth=depth, R=74.93)
    # end

    file = open(filename + ".txt", 'w')
    info = f"Depth: {depth} , Diameter: {diameter} , Volume: {volume}, Ridge Height: {ridge_height}, Radius of curvature : {roc}, Radius of Sphere : {ros}, Radius of Sphere Fit : {roS}, Standard deviation: {err}"
    file.write(info)
    file.close()

    new_row = {'ID' : filename, 'Depth' : depth, 'Diameter' : diameter, 'Volume' : volume, 'Ridge Height': ridge_height, "Radius of curvature" : roc, "Radius of Sphere" : ros, "Radius of Sphere Fit" : roS, "Standard deviation" : err}
    df.loc[len(df)] = new_row

    os.chdir(cwd)
    os.chdir(path)

os.chdir(cwd)
df.to_csv("Testing/Foam Sphere Analysis/All - mean/foam_sphere_analysis.csv")

print("Finished")