# on_csv(filename, savepath): Interpolate algorithm error due to analyis in mm
# total_error(filename, savepath): Combine algorithm and camera error in mm


import pandas as pd
import numpy as np
from scipy.interpolate import griddata

# Define the depth and diameter values
depth_values = [5, 10, 15, 20, 25, 30, 35]
diameter_values = [40, 80, 120, 160, 200, 240]

# Read the CSV file into a DataFrame
depth_table = pd.read_csv("Error Lookup function/depth_error.csv", header=None)

# Set the column and row names
depth_table.columns = depth_values
depth_table.index = diameter_values

# Read the CSV file into a DataFrame
diameter_table = pd.read_csv("Error Lookup function/diameter_error.csv", header=None)

# Set the column and row names
diameter_table.columns = depth_values
diameter_table.index = diameter_values

# Read the CSV file into a DataFrame
volume_table = pd.read_csv("Error Lookup function/volume_error.csv", header=None)

# Set the column and row names
volume_table.columns = depth_values
volume_table.index = diameter_values


def linear_interpolation(depth, diameter):
    # Create the grid of diameter and depth values
    grid_diameter, grid_depth = np.meshgrid(diameter_table.index, diameter_table.columns)

    # Flatten the grid and corresponding values
    points = np.column_stack((grid_diameter.ravel(), grid_depth.ravel()))
    values_depth = depth_table.values.ravel()
    values_diameter = diameter_table.values.ravel()
    values_volume = volume_table.values.ravel()

    # Perform linear interpolation for depth
    depth_error = griddata(points, values_depth, (diameter, depth), method='linear')

    # Perform linear interpolation for diameter
    diameter_error = griddata(points, values_diameter, (diameter, depth), method='linear')

    # Perform linear interpolation for volume
    volume_error = griddata(points, values_volume, (diameter, depth), method='linear')

    return abs(float(depth_error)), abs(float(diameter_error)), abs(float(volume_error))


def on_csv(filename, savepath):
    df = pd.read_csv(filename)
    df["Depth Error"] = ''
    df["Diamter Error"] = ''
    df["Volume Error"] = ''
    df["Ridge Height Error"] = ''
    for index, row in df.iterrows():
        dep = row["Depth"]
        dia = row["Diameter"]
        if dep < 5:
            dep = 5
        elif dep > 35:
            dep = 35
        if dia < 40:
            dia = 40
        elif dia > 240:
            dia = 240
        depth, diameter, volume = linear_interpolation(dep, dia)
        df.loc[index, "Depth Error"], df.loc[index, "Diameter Error"], df.loc[index, "Volume Error"] = depth, diameter, volume
        df.loc[index, "Ridge Height Error"] = depth
    df.to_csv(filename, index_col=False)
    print("File Generated")

def total_error(filename, savepath):
    df = pd.read_csv(filename) # for analysis+error.csv
    df["Depth Total Error"] = ''
    df["Diamter Total Error"] = ''
    df["Volume Total Error"] = ''
    df["Ridge Height Total Error"] = ''
    for index, row in df.iterrows():
        dep = row["Depth"]
        dia = row["Diameter"]
        vol = row["Volume"]
        rh = row["Ridge Height"]
        dep_err = row["Depth Error"]
        dia_err = row["Diameter Error"]
        vol_err = row["Volume Error"]
        rh_err = row["Ridge Height Error"]

        delta_camera = 0.7 # Calculated value for Microsoft Azure Kinect DK

        delta_dep = delta_camera 
        delta_dia = 2*delta_camera
        delta_volume = (3*delta_dep/dep)*vol
        delta_rh = delta_camera

        df.loc[index,"Depth Total Error"] = np.sqrt(delta_dep**2 + dep_err**2)
        df.loc[index,"Diamter Total Error"] = np.sqrt(dia_err**2 + delta_dia**2)
        df.loc[index,"Volume Total Error"] = np.sqrt(vol_err**2 + delta_volume**2)
        df.loc[index,"Ridge Height Total Error"] = np.sqrt(delta_rh**2 + rh_err**2)
    df.to_csv(savepath, index_col=False)
    print("File Generated")

# # Generate algorithm error for each instance
# file = "Lab Craters\September 2023 Results/analysis.csv"
# savepath = "Lab Craters\September 2023 Results/errors.csv"
# on_csv(file, savepath)
# ############################################

# Combine algorithm error and camera error for each instance
file = "Lab Craters\September 2023 Results/analysis.csv"
savepath = "Lab Craters\September 2023 Results/final.csv"
total_error(file, savepath)
#############################################
