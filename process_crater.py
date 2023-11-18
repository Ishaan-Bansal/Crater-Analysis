import os
import helper_functions as help

# Write the relative path to the folder you want to load
path = "Lab Craters/November 2023 STLs"
save_path = "Processed"
cwd = os.getcwd()
os.chdir(path)
for filename in os.listdir():
    if filename[-4:] != ".stl":
        continue
    print("Working on: " + filename)

    mesh = help.load_mesh(filename, trimRadius=250, displayRadius=250)
    mesh.export(save_path + "/" + filename)

print("Finished")