# Interpolate error due to analyis in mm


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


def on_csv(filename):
    df = pd.read_csv(filename)
    err_df = df.copy(deep=True)
    for index, row in df.iterrows():
        dep = row["Depth"]
        dia = row["Diameter"]
        if dep < 5:
            dep = 5
        if dia < 40:
            dia = 40
        depth, diameter, volume = linear_interpolation(dep, dia)
        err_df.loc[index, "Depth"], err_df.loc[index, "Diameter"], err_df.loc[index, "Volume"] = depth, diameter, volume
        err_df.loc[index, "Ridge Height"] = depth
    err_df.to_csv("Lab Craters\Combined Results/errors.csv")
    print("File Generated")

# file = "Lab Craters\Combined Results/analysis.csv"
# on_csv(file)

linear_interpolation(20, 30)