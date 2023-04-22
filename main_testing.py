from helper_functions import lowest_point
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

dimDF = {}

# Write the relative path to the folder you want to load
path = "Test Crater Type 2 STLs"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    if filename[-4:] != ".stl":
        continue
    print("Working on: " + filename)

    mesh = trimesh.load_mesh(filename)
    help.move(mesh, lowest_point(mesh))
    
    filename = filename[:-4]

    radius_row, depth_col, curv_3dim = help.from_filename_t2(filename)

    if not help.find_in_dict(dimDF, curv_3dim):
        dimDF[curv_3dim] = pd.DataFrame()

    if depth_col not in dimDF[curv_3dim].columns:
        dimDF[curv_3dim][depth_col] = pd.Series(dtype=object)
    if not dimDF[curv_3dim].index.isin([radius_row]).any():
        dimDF[curv_3dim].loc[radius_row] = [str()] * len(dimDF[curv_3dim].columns) 

    os.chdir(cwd)
    os.chdir("Test_Crater_Props")
    os.mkdir(filename)
    os.chdir(filename)

    depth, diameter, volume = rtc.crater_properties(mesh)

    file = open(filename + ".txt", 'w')
    info = f"Depth: {depth} , Diameter: {diameter} , Volume: {volume}"
    file.write(info)
    file.close()

    dimDF[curv_3dim].loc[radius_row, depth_col] = info

    os.chdir(cwd)
    os.chdir(path)

os.chdir(cwd)

for key in dimDF.keys():
    dimDF[key].to_csv("Crater_Properties/c" + key + "_test_t2_craters.csv")

print("Finished")