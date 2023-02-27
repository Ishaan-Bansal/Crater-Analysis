from audioop import cross
from stl import Mesh
import numpy as np
import trimesh
import statistics as stats
import math


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

# Check if a given point is within the given radius from the given point
def withinBounds(point1, centroid, radius):
    point = pointsToPoint(point1)  # Convert N-9 point to N-3 point
    vector = point - centroid  # Distance between point and centroid
    # Return boolean whether the distance is within the radius
    return np.linalg.norm(vector) < radius


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
