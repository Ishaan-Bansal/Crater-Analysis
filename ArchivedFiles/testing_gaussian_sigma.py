import helper_functions as help
import trimesh
import numpy as np
import sympy
from stl import Mesh
import matplotlib.pyplot as plt
import statistics as stats
import scipy.ndimage
from scipy.signal import medfilt2d
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_laplace
import ray_tracing_Z_Slices as rtz
import ray_tracing_X_Slices as rtx
from skimage import data, color
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.util import img_as_ubyte
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar

if __name__ == '__main__':
    # Write the relative path to the file you want to load
    filename = 'Lab Craters/November 2023 STLs/Processed/2023_03_15_6Torr_h3_860gs_crater11.stl'
    # Radius for trimming_package; Effects the normal vector and rotation matrix
    trimRadius = 250
    # Radius for trimming the mesh around a point; Effects the histogram and slice
    displayRadius = 250
    # The bounds of the square in which the rays are limited
    bounds = 150
    # The spacing between each ray
    spacing = 1

    mesh = trimesh.load(filename)
    help.move(mesh, help.lowest_point(mesh))

    ray_directions_array = np.array([])
    ray_origins_array = np.array([])
    index_tri_array = np.array([], dtype=int)
    locations_array = np.array([])
    x_locations = []
    y_locations = []
    z_locations = []

    divergence_vis = []

    print("Starting ray tracing...")

    for y in range(-bounds, bounds+spacing, spacing):
        locations_row = []
        for x in range(-bounds, bounds+spacing, spacing):
            ray_directions = np.array([[0, 0, -1]])
            ray_origins = np.array([[x, y, 200]])
            divergence_vis.append(np.array([0, 0, 0]))
            locations, index_ray, index_tri = mesh.ray.intersects_location(
                ray_origins=ray_origins, ray_directions=ray_directions)
            locations_row.append(locations[0])
            ray_directions_array = np.append(ray_directions_array, ray_directions)
            ray_origins_array = np.append(ray_origins_array, ray_origins)
            index_tri_array = np.append(index_tri_array, index_tri)

        locations_row = np.array(locations_row)
        x_locations.append(locations_row[:, 0])
        y_locations.append(locations_row[:, 1])
        z_locations.append(locations_row[:, 2])

    print("Ray tracing finished...")

    # index_tri_array = index_tri_array.astype(int)
    x_locations = np.array(x_locations, dtype=float)
    y_locations = np.array(y_locations, dtype=float)
    z_locations = np.array(z_locations, dtype=float)
    divergence_vis = np.array(divergence_vis)

    diameter_list = []
    sigma_list = np.linspace(0, 8, 100)
    for s in sigma_list:
        z_locations = gaussian_filter(z_locations, sigma=s)

        gradient = np.gradient(z_locations)
        gradient_norm = np.sqrt(gradient[0]**2 + gradient[1]**2)

        d2z_dx2 = np.gradient(gradient[0])[0]
        d2z_dy2 = np.gradient(gradient[1])[1]
        double_gradient_norm = d2z_dx2 + d2z_dy2
        ridge_indices = []
        ridge_z = []
        radius_x = int(z_locations.shape[0]/2)
        radius_y = int(z_locations.shape[1]/2)
        dgn = []

        for theta_deg in range(0, 360):
            # -- Extract the line...
            # Make a line with "num" points...
            theta_rad = np.radians(theta_deg)
            # These are in _pixel_ coordinates!!
            x1, y1 = radius_x - radius_x * \
                np.cos(theta_rad), radius_y - radius_y*np.sin(theta_rad)
            x0, y0 = radius_x, radius_y
            length = int(np.round(np.hypot(x1-x0, y1-y0)))
            x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
            mask = np.sqrt((x-radius_x)**2 + (y-radius_y) **
                        2) >= (0.2*(radius_x+radius_y)/2)
            x, y = x[mask], y[mask]

            # Extract the values along the line
            x = np.round(x).astype(int)
            y = np.round(y).astype(int)
            z_line = double_gradient_norm[x, y]

            index = np.argmin(z_line)
            pair = np.array([y[index], x[index]])
            ridge_indices.append(pair)
            height_line = z_locations[x, y]
            ridge_z.append(height_line[index])


        ridge_indices = np.array(ridge_indices)
        ridge_z = np.array(ridge_z)

        x, y = [], []
        for i in ridge_indices:
            x.append(x_locations[i[0]][i[1]])
            y.append(y_locations[i[0]][i[1]])

        x, y, ridge_indices, ridge_z = help.remove_outliers(x,y, ridge_indices, ridge_z)
        
        x, y = np.array(x), np.array(y)

        crater_start = np.mean(ridge_z)
        a, b, d = 0, 0, crater_start
        plane = help.create_trimesh_plane(a, b, d)

        b = x**2 + y**2
        M_a = np.ones((x.size, 3))
        M_a[:, 0] = x
        M_a[:, 1] = y
        sol = np.linalg.lstsq(M_a, b, rcond=None)

        x_c = sol[0][0]/2
        y_c = sol[0][1]//2
        radius = np.sqrt(sol[0][2] + x_c**2 + y_c**2)
        diameter_list.append(radius*2)

    # fig,ax = plt.subplots()
    # ax.plot(sigma_list, diameter_list, "rx")

    n, bins, patches = plt.hist(
    diameter_list, bins=10, density=True, facecolor='b')
    diameter = bins[n.argmax()]

    # ax.axhline(y=diameter, color='g', linestyle='dashed')
    plt.show()

    # def func(x,a,b,c,d,e):
    #     return a*x**4 + b*x**3 + c*x**2 + d*x + e
    # popt, pcov = curve_fit(func, sigma_list, diameter_list)
    # fig,ax = plt.subplots()
    # ax.plot(sigma_list, diameter_list, "rx")
    # ax.plot(sigma_list, func(sigma_list, *popt), 'r-',
    #         label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f' % tuple(popt))
    # plt.legend()

    # a,b,c,d,e = popt
    # saddle1 = (-3*2*b - np.sqrt((3*2*b)**2 - 4*(4*3*a)*(2*c))) / (2*4*3*a) # solve d2f/dx2 = 0
    # saddle2 = (-3*2*b + np.sqrt((3*2*b)**2 - 4*(4*3*a)*(2*c))) / (2*4*3*a) # solve d2f/dx2 = 0

    # saddle = min(saddle1,saddle2)
    # if 4*a*(saddle1)**3 + 3*b*(saddle1)**2 + 2*c*(saddle1) + d*saddle1 == 0: # check df/dx = 0
    #     saddle = saddle1
    # elif 4*a*(saddle2)**3 + 3*b*(saddle2)**2 + 2*c*(saddle2) + d*saddle2 == 0: # check df/dx = 0 
    #     saddle = saddle2
    # ax.plot(saddle1, func(saddle1, *popt), 'gx')
    # ax.plot(saddle2, func(saddle2, *popt), 'gx')
    # ax.plot(saddle, func(saddle, *popt), 'bx')
    plt.show()

def crater_properties(mesh, bounds):
    # Bounds: Half the lenght of the square in which the rays are limited
    # The spacing between each ray
    spacing = 1

    ray_directions_array = np.array([])
    ray_origins_array = np.array([])
    index_tri_array = np.array([], dtype=int)
    locations_array = np.array([])
    x_locations = []
    y_locations = []
    z_locations = []

    divergence_vis = []

    for y in range(-bounds, bounds+spacing, spacing):
        locations_row = []
        for x in range(-bounds, bounds+spacing, spacing):
            ray_directions = np.array([[0, 0, -1]])
            ray_origins = np.array([[x, y, 200]])
            divergence_vis.append(np.array([0, 0, 0]))
            locations, index_ray, index_tri = mesh.ray.intersects_location(
                ray_origins=ray_origins, ray_directions=ray_directions)
            locations_row.append(locations[0])
            ray_directions_array = np.append(ray_directions_array, ray_directions)
            ray_origins_array = np.append(ray_origins_array, ray_origins)
            index_tri_array = np.append(index_tri_array, index_tri)

        locations_row = np.array(locations_row)
        x_locations.append(locations_row[:, 0])
        y_locations.append(locations_row[:, 1])
        z_locations.append(locations_row[:, 2])

    x_locations = np.array(x_locations, dtype=float)
    y_locations = np.array(y_locations, dtype=float)
    z_locations = np.array(z_locations, dtype=float)
    z_locations = gaussian_filter(z_locations, sigma=5)
    divergence_vis = np.array(divergence_vis)

    gradient = np.gradient(z_locations)
    gradient_norm = np.sqrt(gradient[0]**2 + gradient[1]**2)

    double_gradient = np.gradient(gradient_norm)
    double_gradient_norm = np.sqrt(double_gradient[0]**2 + double_gradient[1]**2)

    fig = plt.figure()
    img = plt.imshow(double_gradient_norm, interpolation='none', cmap='turbo')
    cbar = plt.colorbar(img)

    ridge_indices = []
    ridge_z = []
    radius_x = int(z_locations.shape[0]/2)
    radius_y = int(z_locations.shape[1]/2)

    for theta_deg in range(0, 360):
        # -- Extract the line...
        # Make a line with "num" points...
        theta_rad = np.radians(theta_deg)
        # These are in _pixel_ coordinates!!
        x1, y1 = radius_x - radius_x * \
            np.cos(theta_rad), radius_y - radius_y*np.sin(theta_rad)
        x0, y0 = radius_x, radius_y
        length = int(np.round(np.hypot(x1-x0, y1-y0)))
        x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
        mask = np.sqrt((x-radius_x)**2 + (y-radius_y) **
                    2) >= (0.5*(radius_x+radius_y)/2)
        x, y = x[mask], y[mask]

        # Extract the values along the line
        x = np.round(x).astype(int)
        y = np.round(y).astype(int)
        z_line = double_gradient_norm[x, y]

        index = np.argmax(z_line)
        pair = np.array([y[index], x[index]])
        ridge_indices.append(pair)
        height_line = z_locations[x, y]
        ridge_z.append(height_line[index])

    ridge_indices = np.array(ridge_indices)
    ridge_z = np.array(ridge_z)
    plt.plot(ridge_indices[:, 0], ridge_indices[:, 1], 'r1')

    crater_start = np.mean(ridge_z)

    n, bins, patches = plt.hist(
    mesh.vertices[:,2], bins=300, density=True, facecolor='b', alpha=0.75)
    depth = bins[n.argmax()]

    x, y = [], []
    for i in ridge_indices:
        x.append(x_locations[i[0]][i[1]])
        y.append(y_locations[i[0]][i[1]])

    x, y, ridge_indices, ridge_z = help.remove_outliers(x,y, ridge_indices, ridge_z)
    x, y = np.array(x), np.array(y)

    plt.plot(ridge_indices[:, 0], ridge_indices[:, 1], 'kx')
    plt.title('Ridge Detection')
    plt.savefig("Ridge_Detection.svg")
    plt.close()

    b = x**2 + y**2
    M_a = np.ones((x.size, 3))
    M_a[:, 0] = x
    M_a[:, 1] = y
    sol = np.linalg.lstsq(M_a, b, rcond=None)

    x_c = sol[0][0]/2
    y_c = sol[0][1]//2
    radius = np.sqrt(sol[0][2] + x_c**2 + y_c**2)

    figure, axes = plt.subplots()
    Drawing_uncolored_circle = plt.Circle((x_c, y_c), radius, fill=False)
    plot = plt.plot(x, y, 'k.', linewidth=2)

    axes.set_xlim(-300, 300)
    axes.set_ylim(-300, 300)
    axes.add_patch(Drawing_uncolored_circle)
    plt.title('Circle Fit')
    plt.savefig("Circle Fit.svg")
    plt.close()

    mesh_volume = help.crater_volume_tetra(mesh.triangles, radius, crater_start=depth)

    ridge_height = crater_start - depth
    return depth, 2*radius, mesh_volume, ridge_height
