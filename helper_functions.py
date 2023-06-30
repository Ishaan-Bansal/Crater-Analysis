from stl import Mesh
import numpy as np
import trimesh
import statistics as stats
import math
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# def mesh_to_vectors(mesh): # Archived
#     # Convert mesh to vertices of the trianglar faces
#     data = mesh.vectors

#     vectorsOut = []  # Array of vectors to output

#     # Access every triangle in the mesh
#     for triangle in data:
#         # Vector in the plane of the triangle
#         planeVector1 = triangle[1] - triangle[0]
#         # Vector in the plane of the triangle
#         planeVector2 = triangle[2] - triangle[1]

#         # Calculate the normal vector

#         # Normal of the plane of the triangle
#         planeNormalVector = np.cross(planeVector1, planeVector2)
#         # Unit Normal of the plane of the triangle
#         planeUnitNormal = planeNormalVector / np.linalg.norm(planeNormalVector)

#         # Add unit normal vector to output array
#         vectorsOut.append(planeNormalVector)

#     # Convert to np.array
#     vectorsOut = np.array(vectorsOut)

#     return vectorsOut

# Check if a given point is within the given spherical radius from the given point
def withinBounds(point1, centroid, radius):
    point = pointsToPoint(point1)  # Convert N-9 point to N-3 point
    vector = point - centroid  # Distance between point and centroid
    # Return boolean whether the distance is within the radius
    return np.linalg.norm(vector) < radius

# Check if a given point is within the given circular radius from the given point
def withinCircularRadius(point, centroid, radius):
    vector = point - centroid  # Distance between point and centroid
    # Return boolean whether the distance is within the radius
    return np.sqrt(vector[0]**2+vector[1]**2) < radius

# Trim a circle around the centroid of the numpy.stl mesh
def trimCircle(mesh, radius):
    meshObjTrimesh = trimesh.Trimesh(
        **trimesh.triangles.to_kwargs(mesh.vectors))  # Convert to stl.Mesh object to Trimesh object
    # Convert numpy array to python standard array
    meshPointArr = mesh.points.tolist()
    for i in range(len(meshPointArr)):  # Loop through the points
        if not withinBounds(meshPointArr[i], meshObjTrimesh.centroid, radius):
            # If outside bounds then reassign the point to orgin
            meshPointArr[i] = np.array(
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    # Update the vertices of the mesh to the trimmed vertices
    mesh.points = np.array(meshPointArr)
    return mesh  # Returns stl.mesh object


def trimCircleGivenPoint(mesh, point, radius):  # Works with stl.mesh object
    meshPointArr = mesh.points.tolist()
    for i in range(len(meshPointArr)):  # Loop through the points
        if not withinBounds(meshPointArr[i], point, radius):
            # If outside bounds then reassign the point to orgin
            meshPointArr[i] = np.array(
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    # Update the vertices of the mesh to the trimmed vertices
    mesh.points = np.array(meshPointArr)
    return mesh


def pointsToPoint(points):  # N-9 point to N-3 point array for numpy.stl mesh
    x = (points[0] + points[3] + points[6])/3
    y = (points[1] + points[4] + points[7])/3
    z = (points[2] + points[5] + points[8])/3
    return np.array([x, y, z])


# Convert a trianglular face into a one dimensional point
def triangleToPoint(trianle):  # For numpy.stl mesh
    x = trianle[0][0] + trianle[1][0] + trianle[2][0]
    y = trianle[0][1] + trianle[1][1] + trianle[2][1]
    z = trianle[0][2] + trianle[1][2] + trianle[2][2]
    return np.array([x, y, z])


# Input: unit vector array atribute of stl.Mesh object; Output: normal vector to the flat plane
def unitVectorsToZplane(vectors):  # Finds the average normal vector
    # Initialize the empty array components of the unit normal vectors
    x_comp = []
    y_comp = []
    z_comp = []

    # Add all the components to the empty array
    for vector in vectors:
        x_comp.append(vector[0])
        y_comp.append(vector[1])
        z_comp.append(vector[2])

    # Convert the python array to Numpy array
    x_comp = np.array(x_comp)
    y_comp = np.array(y_comp)
    z_comp = np.array(z_comp)

    # Find the average of the components
    y = stats.mean(y_comp)
    z = stats.mean(z_comp)
    x = stats.mean(x_comp)

    # Create the normal vector
    vec = np.array([x, y, z])
    return vec/np.linalg.norm(vec)  # Return the normalized normal vector


def lowest_point(mesh):  # Finds the lowest point in the mesh
    low = mesh.vertices[0]
    for m in mesh.vertices:
        if m[2] < low[2]:
            low = m
    return low


def plane():  # Creates the x-y plane
    data = np.zeros(2, dtype=Mesh.dtype)
    # Plane
    data['vectors'][0] = np.array([[0, 1, 0],
                                  [1, 0, 0],
                                  [0, 0, 0]])
    data['vectors'][1] = np.array([[1, 0, 0],
                                  [0, 1, 0],
                                  [1, 1, 0]])
    data['vectors'] *= 100
    ouput_mesh = Mesh(data, remove_empty_areas=False)
    return trimesh.Trimesh(**trimesh.triangles.to_kwargs(ouput_mesh.vectors))


def vectorToPlane(vector):  # Create an arrow mesh pointing in the direction of the vector
    data = np.zeros(3, dtype=Mesh.dtype)
    # Plane
    data['vectors'][0] = np.array([[0, 0, 0],
                                   [1, 0, 0],
                                   [1, 1, 1]])
    data['vectors'][1] = np.array([[0, 0, 0],
                                   [0, 1, 0],
                                   [1, 1, 1]])
    data['vectors'][2] = np.array([[0, 0, 0],
                                   [0, 0, 1],
                                   [1, 1, 1]])
    data['vectors'] *= vector
    data['vectors'] *= 50
    return trimesh.Trimesh(**trimesh.triangles.to_kwargs(Mesh(data, remove_empty_areas=False).vectors))


def rotationAngle(vector, axis):
    projOntoPlane = vector - vector*np.array(axis)
    sinA = np.dot(projOntoPlane, [0, 0, 1]) / np.linalg.norm(projOntoPlane)
    return np.arcsin(sinA)


def planeVector(mesh):  # Creates a plane from given normal vector
    nV = unitVectorsToZplane(mesh.face_normals)
    A = np.array([[1, 0, 0],
                  [np.cos(np.pi/2), -np.sin(np.pi/2), 0],
                  [np.sin(np.pi/2), np.cos(np.pi/2), 0]])
    return A @ nV


def rotationMatrix(vector1, vector2):  # Rotate from one vector basis to another
    cross = np.cross(vector1, vector2)
    u = cross/np.linalg.norm(cross)  # Rotation Axis

    dot = np.dot(vector1, vector2)
    rotation_angle = math.atan2(np.linalg.norm(cross), dot)  # Rotation Angle

    cos = np.cos(rotation_angle)
    cos_1 = 1 - cos
    sin = np.sin(rotation_angle)

    R = np.array([[cos+(u[0]**2)*(cos_1), u[0]*u[1]*cos_1-u[2]*sin, u[0]*u[2]*cos_1+u[1]*sin, 0],
                  [u[0]*u[1]*cos_1+u[2]*sin, cos +
                      (u[1]**2)*(cos_1), u[1]*u[2]*cos_1-u[0]*sin, 0],
                  [u[2]*u[0]*cos_1-u[1]*sin, u[2]*u[1]*cos_1 +
                      u[0]*sin, cos+(u[2]**2)*(cos_1), 0],
                  [0, 0, 0, 1]])

    return R


# Rotates a point from one basis to another
def rotation_matrix_point(vector1,  vector2):
    cross = np.cross(vector1, vector2)
    u = cross/np.linalg.norm(cross)  # Rotation Axis
    # print("u = " + str(u))

    dot = np.dot(vector1, vector2)
    rotation_angle = math.atan2(np.linalg.norm(cross), dot)  # Rotation Angle
    # print("rotation_angle = " + str(rotation_angle))

    cos = np.cos(rotation_angle)
    cos_1 = 1 - cos
    sin = np.sin(rotation_angle)

    R = np.array([[cos+(u[0]**2)*(cos_1), u[0]*u[1]*cos_1-u[2]*sin, u[0]*u[2]*cos_1+u[1]*sin],
                  [u[0]*u[1]*cos_1+u[2]*sin, cos +
                      (u[1]**2)*(cos_1), u[1]*u[2]*cos_1-u[0]*sin],
                  [u[2]*u[0]*cos_1-u[1]*sin, u[2]*u[1]*cos_1 +
                      u[0]*sin, cos+(u[2]**2)*(cos_1)]])
    return R


# Moves a Trimesh object along a direction vector in its coordinate space
def move(trimeshObj, distance_vector):
    trimeshObj.vertices -= distance_vector

# Moves a Numpy STL object along a direction vector in its coordinate space
def move_np_stl(stlObj, distance_vector):
    stlObj.translate(-distance_vector)

def zPoints(mesh):  # Input: trimeshObj
    output = mesh.vertices
    return output[:, 2]  # Output: Array of z coordinates


# Gets the lowest point of a mesh given the filename
def lowest_point_file(filename, radius):
    your_mesh = Mesh.from_file(filename)  # Load the mesh
    # Trim the mesh around the centroid(Trimesh property)
    trimCircle(your_mesh, radius)
    # Covert to Trimesh object
    trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
    # Rotate the mesh to correct the x-y plane
    trimmed.apply_transform(rotation_matrix_file(filename, radius))
    low_point = lowest_point(trimmed)  # Find the lowest point of the mesh
    return low_point


# Outputs the rotation matrix given a filename
def rotation_matrix_file(filename, radius):
    your_mesh = Mesh.from_file(filename)  # Load the mesh
    trimCircle(your_mesh, radius)  # Trim around the centroid(Trimesh property)
    # Covert to Trimesh object
    trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
    # Find the normal vector using previously defined function
    normalVec = unitVectorsToZplane(trimmed.face_normals)
    # Return the rotation matrix using normal vector and z-axis
    return rotationMatrix(normalVec, np.array([0, 0, 1]))


# Outputs the rotation matrix for a point given a filename
def rotation_matrix__point_file(filename, radius):
    your_mesh = Mesh.from_file(filename)  # Load the mesh
    trimCircle(your_mesh, radius)  # Trim around the centroid(Trimesh property)
    # Covert to Trimesh object
    trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
    normalVec = unitVectorsToZplane(
        trimmed.face_normals)  # Find the normal vector
    # Return the rotation matrix for a point using the Z-axis and normalVec
    return rotation_matrix_point(np.array([0, 0, 1]), normalVec)


def rotation_matrix_file_package(filename, radius):
    your_mesh = Mesh.from_file(filename)  # Load the mesh
    trimCircle(your_mesh, radius)  # Trim around the centroid(Trimesh property)
    # Covert to Trimesh object
    trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
    normalVec = unitVectorsToZplane(
        trimmed.face_normals)  # Find the normal vector
    # Return the rotation matrix for a point using the Z-axis and normalVec
    return rotationMatrix(normalVec, np.array([0, 0, 1])), rotation_matrix_point(np.array([0, 0, 1]), normalVec)

# The goal of this function is to limit the computation time by reducing the number of times the mesh is loaded


# Returns the lowest point in the mesh, rotation matrix, and rotation matrix for a point
def trimming_package(filename, radius):
    your_mesh = Mesh.from_file(filename)  # Load the mesh
    trimCircle(your_mesh, radius)  # Trim around the centroid(Trimesh property)
    # Covert to Trimesh object
    trimmed = trimesh.Trimesh(**trimesh.triangles.to_kwargs(your_mesh.vectors))
    normalVec = unitVectorsToZplane(
        trimmed.face_normals)  # Find the normal vector
    rotation_matrix, rotation_matrix_point = rotation_matrix_file_package(
        filename, radius)
    trimmed.apply_transform(rotation_matrix)
    low_point = lowest_point(trimmed)  # Find the lowest point of the mesh
    return normalVec, low_point, rotation_matrix, rotation_matrix_point


def histPlot(n, bins):  # Rolling mean graphing method for histograms
    nArr = []
    binsArr = []
    for i in range(n.size):
        sum = np.sum(n[i:i+3])/3
        nArr.append(sum)
        binsArr.append(bins[i+1])
    return (np.array(nArr), np.array(binsArr))


def slicer(mesh, plane, buffer):  # Input: Trimesh object
    points = mesh.vertices
    # Boolean array of the within the bounds of the y-z plane
    if plane == 'x':
        planeN = 0
    elif plane == 'y':
        planeN = 1
    elif plane == 'z':
        planeN = 2
    mask1 = -buffer < points[:, planeN]
    mask2 = buffer > points[:, planeN]
    mask = []
    for i in range(len(points)):
        mask.append(mask1[i] and mask2[i])
    points = points[mask]  # Points within the bounds of the y-z plane
    if plane == 'x':
        # Returns: y-coordinates and z-coordinates
        return (points[:, 1], points[:, 2])
    elif plane == 'y':
        # Returns: x-coordinates and z-coordinates
        return (points[:, 0], points[:, 2])
    elif plane == 'z':
        # Returns: x-coordinates and y-coordinates
        return (points[:, 0], points[:, 1])


def create_trimesh_plane(a, b, d):
    # vertices of the cube
    vertices = np.array([250, 250, 250*a+250*b+d, -250, 250, -250*a+250*b+d, 250, -250, 250*a-250*b+d, -250, -250, -250*a-250*b+d],
                        order='C',
                        dtype=object).reshape((-1, 3))
    vertices = np.array(vertices, dtype=np.float64)
    # hardcoded face indices
    faces = [0, 1, 2, 2, 3, 1, 1, 3, 2, 2, 1, 0]
    faces = np.array(faces, order='C', dtype=np.int64).reshape((-1, 3))
    # hardcoded face normals
    face_normals = [a, b, -1, a, b, -1, a, b, 1, a, b, 1]
    face_normals = np.asanyarray(face_normals,
                                 order='C',
                                 dtype=np.float64).reshape(-1, 3)
    # Create Trimesh object
    plane = trimesh.Trimesh(vertices=vertices,
                            faces=faces,
                            face_normals=face_normals,
                            process=False)
    return plane


def crater_volume(array, radius, spacing, crater_top):
    height = array.shape[0]
    width = array.shape[1]
    xc = int(height/2)
    yc = int(width/2)
    volume = 0
    area = spacing**2
    for i in range(height):
        for j in range(width):
            r = np.sqrt((i-yc)**2 + (j-xc)**2)
            if r <= radius:
                volume += area*(crater_top-array[i][j])
    return volume


def cut_top(mesh, max):
    arr = []
    for x in mesh.vertices:
        if x[2] > max:
            x = [np.nan, np.nan, np.nan]
        arr.append(x)
    return arr

def crater_volume_tetra(mesh_triangles, radius, crater_start):
    volume = 0
    once = True
    for triangle in mesh_triangles:
        x, y, z = np.mean(triangle[:,0]), np.mean(triangle[:,1]), np.mean(triangle[:,2])
        if not np.sqrt(x**2 + y**2) < radius:
            continue
        if not z < crater_start:
            continue
        matrix = np.ones((4,4))
        matrix[0:3, 0:3] = triangle
        matrix[3, 0:3] = [0,0,crater_start]
        volume += (np.abs(np.linalg.det(matrix))/6)
        # if once:
        #     print(matrix)
        #     once = False
        #     print(volume)
    return volume


def load_mesh(filename, trimRadius, displayRadius):
    # trimRadius: Radius for trimming_package; Effects the normal vector and rotation matrix
    # displayRadius: Radius for trimming the mesh around a point; Effects the histogram and slice

    normal_vector, lowest_point, rotation_matrix, rotation_matrix_point = trimming_package(
        filename, trimRadius)
    # Find the lowest point in the unrotated basis
    lowest_point_r = rotation_matrix_point @ lowest_point

    # Trim around the lowest_point
    new_mesh = Mesh.from_file(
        filename)
    trimCircleGivenPoint(new_mesh, lowest_point_r, displayRadius)
    mesh = trimesh.Trimesh(**trimesh.triangles.to_kwargs(new_mesh.vectors))
    mesh.remove_infinite_values()
    mesh.apply_transform(rotation_matrix)
    move(mesh, lowest_point)
    return mesh

def from_filename(filename):
    r = ""
    d = ""
    radius = False
    depth = False
    for char in filename:
        if char == "r":
            radius = True
            continue
        elif char == "_":
            radius = False
            continue
        elif char == "d":
            depth = True
            continue
        
        if radius:
            r += char
        elif depth:
            d += char
    return r, d

def from_filename_t2(filename):
    r = ""
    d = ""
    c = ""
    radius = False
    depth = False
    curv = False
    for char in filename:
        if char == "r":
            radius = True
            continue
        elif radius and char == "_":
            radius = False
            continue
        elif char == "d":
            depth = True
            continue
        elif depth and char == "_":
            depth = False
            continue
        elif char == "c":
            curv = True
            continue
        elif curv and char == "_":
            break
        
        if radius:
            r += char
        elif depth:
            d += char
        elif curv:
            c += char
    return r, d, c

def find_in_dict(dict, idx):
    for key in dict.keys():
        if key == idx:
            return True
    return False

def remove_outliers(x , y, ridge_indices, ridge_z):
    r = []
    for i in range(len(x)):
        r.append(np.sqrt(x[i]**2 + y[i]**2))
    x_bar = np.std(r)
    remove_x = []
    remove_y = []
    ind = []

    if (x_bar < 0.2*np.mean(r)):
        for i in range(len(x)):
            if r[i] > np.mean(r) + 2*x_bar or r[i] < np.mean(r) - 2*x_bar:
                remove_x.append(x[i])
                remove_y.append(y[i])
                ind.append(i)
    else: 
        for i in range(len(x)):
            if r[i] > np.mean(r) + x_bar or r[i] < np.mean(r) - x_bar:
                remove_x.append(x[i])
                remove_y.append(y[i])
                ind.append(i)
    for i in remove_x:
        x.remove(i)
    for i in remove_y:
        y.remove(i)
    return x, y, np.delete(ridge_indices, ind, 0), np.delete(ridge_z, ind, 0)

def radius_of_curvature(arr, depth):
    x, y = arr[:,0] , arr[:,2]

    # mask1 = -diameter/2 < x
    # mask2 = diameter/2 > x
    mask1 = -depth < y
    mask2 = depth > y
    mask = []
    for i in range(len(x)):
        mask.append(mask1[i] and mask2[i])

    x, y = x[mask], y[mask]

    b = x**2 + y**2
    M_a = np.ones((x.size, 3))
    M_a[:, 0] = x
    M_a[:, 1] = y
    sol = np.linalg.lstsq(M_a, b, rcond=None)

    x_c = sol[0][0]/2
    y_c = sol[0][1]//2
    radius = np.sqrt(sol[0][2] + x_c**2 + y_c**2)
    return radius

def radius_of_sphere(depth, diameter):
    r = sp.symbols('r')
    eq1 = sp.Eq( (diameter/2)**2 + (r - depth)**2, r**2)
    return sp.solve(eq1, r)

def sphereFit(mesh, radius, depth, R):
    arr = []
    for i in mesh.vertices:
        if i[2] < depth and withinCircularRadius(i, np.array([0,0,0]), radius): 
            arr.append(i)
    arr = np.array(arr)
    r, cx, cy, cz = sFit(arr)
    center = np.array([cx, cy,int(R)])
    errorarr = []
    for i in arr:
        error = np.linalg.norm(i-center) - R 
        errorarr.append(error)
    # fig = plt.figure()
    # ax3D = fig.add_subplot(projection='3d')
    # ax3D.scatter(arr[:,0], arr[:,1], arr[:,2])  
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x=np.cos(u)*np.sin(v)*R
    # y=np.sin(u)*np.sin(v)*R
    # z=np.cos(v)*R
    # x = x + cx
    # y = y + cy
    # z = z + (R)
    # ax3D.plot_wireframe(x, y, z, color="r")
    # ax3D.set_aspect('equal')
    # ax3D.set_zlim3d(-20,100)
    # plt.show()
    return r, cx, cy, cz, np.mean(errorarr)

def sFit(arr):
    #   Assemble the A matrix
    spX = arr[:,0]
    spY = arr[:,1]
    spZ = arr[:,2]
    A = np.zeros((len(spX),4))
    A[:,0] = spX*2
    A[:,1] = spY*2
    A[:,2] = spZ*2
    A[:,3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(spX),1))
    f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
    C, residules, rank, singval = np.linalg.lstsq(A,f, rcond=None)

    #   solve for the radius
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = math.sqrt(t)

    return radius, C[0], C[1], C[2]

if __name__ ==  '__main__':
    filename = "Testing/Foam Sphere STLS Post Proccessed/4 June 2023/04_06_2023_run1_scan1.stl"
    mesh = load_mesh(filename,100,100)
    r, x0, y0, z0, err = sphereFit(mesh, radius=81.49416854901347/2, depth=8.065528520377779, R=74.93)
    print(r)
    print(x0,y0,z0)
    print(err)
    # fig = plt.figure()
    # ax3D = fig.add_subplot(projection='3d')
    # ax3D.scatter(mesh.vertices[:,0], mesh.vertices[:,1], mesh.vertices[:,2])  
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x=np.cos(u)*np.sin(v)*r
    # y=np.sin(u)*np.sin(v)*r
    # z=np.cos(v)*r
    # x = x + x0
    # y = y + y0
    # z = z + (74.93-8.065528520377779)
    # ax3D.plot_wireframe(x, y, z, color="r")
    # ax3D.set_aspect('equal')
    # ax3D.set_zlim3d(-20,100)
    # plt.show()
