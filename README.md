# CraterAnalysis
This purpose of this project is to provide meaningful analysis of craters. The primary requirement is interpretable .stl file. 

The program will provide the depth, width (in production) and slices of the crater as well as a 3D togographical heat map of the crater.

For the depth, width (in production) and slices, run "trimming_around_lowest_point.py".
"trimming_around_lowest_point.py"
Inputs:
    filename (str): Relative path of the .stl file containing the crater
    trimRadius (int/float): Radius for trimming_package; Effects the normal vector and rotation matrix
    dislpayRadius (int/float): Radius for trimming the mesh around a point; Effects the histogram and slice
Outputs:
    depth (float): Printed in terminal
    width (float): Printed in terminal (in production)
    slices: Matplotlib graphs

For the togographical heat map, run "togographic_map.py".
"togographic_map.py"
Inputs:
    filename (str): Relative path of the .stl file containing the crater 
    trimRadius (int/float): Radius for trimming_package; Effects the normal vector and rotation matrix
    dislpayRadius (int/float): Radius for trimming the mesh around a point; Effects the histogram and slice
Outputs:
    Togographical heat map: Pyglet window