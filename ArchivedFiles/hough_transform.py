import pandas as pd
from skimage import data, color, io
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.util import img_as_ubyte
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import medfilt2d


df = pd.read_csv('data2.csv', sep=',', header=None)
double_gradient_norm = df.values
double_gradient_norm = medfilt2d(double_gradient_norm, kernel_size=15)

df2 = pd.read_csv('data3.csv', sep=',', header=None)
gradient_norm = df2.values
gradient_norm = medfilt2d(gradient_norm, kernel_size=33)

fig = plt.figure(2)
plt.imshow(gradient_norm, interpolation='none', cmap='turbo')
plt.show()

df3 = pd.read_csv('data4.csv', sep=',', header=None)
z_locations = df3.values

# # Detect circles of different radii
# hough_radii = np.arange(20, 150, 5)
# hough_res = hough_circle(double_gradient_norm, hough_radii)

# # Select the most prominent circle
# accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)

# # Draw the most prominent circle on the elevations array
# circle_elevations = np.copy(double_gradient_norm)
# for center_y, center_x, radius in zip(cy, cx, radii):
#     circy, circx = circle_perimeter(center_y, center_x, radius,
#                                     shape=double_gradient_norm.shape)
#     circle_elevations[circy, circx] = np.nan

# plt.imshow(circle_elevations, interpolation='none', cmap='turbo')
# plt.show()


ridge_indices = []
ridge_z = []
radius_x = int(z_locations.shape[0]/2)
radius_y = int(z_locations.shape[1]/2)
dgn = []

fig = plt.figure(1)
plt.imshow(double_gradient_norm, interpolation='none', cmap='turbo')

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

#     # print("X:")
#     # print(x)
    index = np.argmax(z_line)
    pair = np.array([y[index], x[index]])
    ridge_indices.append(pair)
    height_line = z_locations[x, y]
    ridge_z.append(height_line[index])

    plt.plot(x, y, 'ko')


ridge_indices = np.array(ridge_indices)
ridge_z = np.array(ridge_z)
print(ridge_indices)
plt.figure(5)
plt.plot(ridge_indices[:, 0], ridge_indices[:, 1], 'ko')
plt.show()
