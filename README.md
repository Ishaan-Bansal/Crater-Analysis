# Crater Analysis (In Development)

The purpose of this project is to provide meaningful analysis of craters. The primary requirement is an interpretable .stl file. 

The program will provide the depth, width (in production) and slices of the crater as well as a 3D togographical heat map of the crater. 

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the following package which are required for this project.

```bash
pip install trimesh
pip install numpy-stl
pip install numpy
```

## Usage

For the depth, width (in production) and slices, run "primary_analysis.py".  
"primary_analysis.py"   
- Inputs:  
    - filename (str): Relative path of the .stl file containing the crater
    - trimRadius (int/float): Radius for trimming_package; Effects the normal vector and rotation matrix
    - dislpayRadius (int/float): Radius for trimming the mesh around a point; Effects the histogram and slice
- Outputs:  
    - depth (float): Printed in terminal
    - width (float): Printed in terminal (in production)
    - slices: Matplotlib graphs

For the togographical heat map, run "togographic_map.py".  
"togographic_map.py"  
- Inputs:
    - filename (str): Relative path of the .stl file containing the crater 
    - trimRadius (int/float): Radius for trimming_package; Effects the normal vector and rotation matrix
    - dislpayRadius (int/float): Radius for trimming the mesh around a point; Effects the histogram and slice
- Outputs:
    - Togographical heat map: Pyglet window

## Implementation

These programs use helper functions in order to analyse the crater. Read "helper_functions.py" as well as [numpy-stl documentation](https://numpy-stl.readthedocs.io/en/latest/) and [Trimesh documentation](https://trimsh.org/index.html) for in-depth understanding of the code.