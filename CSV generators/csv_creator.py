import os
import pandas as pd
import helper_functions as help

df = pd.DataFrame()

# Write the relative path to the folder you want to load
path = "Test Craters TXTs"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    print("Working on: " + filename)

    filename = filename[:-4]

    radius_row, depth_col = help.from_filename(filename)

    if depth_col not in df.columns:
        df[depth_col] = pd.Series(dtype=object)
    if not df.index.isin([radius_row]).any():
        df.loc[radius_row] = [str()] * len(df.columns) 


    file = open(filename + ".txt", 'r')

    df.loc[radius_row, depth_col] = file.read()

os.chdir(cwd)
df.to_csv("test_craters.csv")

print("Finished")