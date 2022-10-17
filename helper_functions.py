from audioop import cross
from stl import Mesh
import numpy as np
import trimesh
import statistics as stats
import math


def mesh_to_vectors(mesh):
    # Convert mesh to vertices of the trianglar faces
    data = mesh.vectors

    vectorsOut = []  # Array of vectors to output

    # Access every triangle in the mesh
    for triangle in data:
        # Vector in the plane of the triangle
        planeVector1 = triangle[1] - triangle[0]
        # Vector in the plane of the triangle
        planeVector2 = triangle[2] - triangle[1]

        # Calculate the normal vector

        # Normal of the plane of the triangle
        planeNormalVector = np.cross(planeVector1, planeVector2)
        # Unit Normal of the plane of the triangle
        planeUnitNormal = planeNormalVector / np.linalg.norm(planeNormalVector)

        # Add unit normal vector to output array
        vectorsOut.append(planeNormalVector)

    # Convert to np.array
    vectorsOut = np.array(vectorsOut)

    return vectorsOut

# Check if a given point is within the given radius from the centroid


def withinBounds(point1, centroid, radius):
    point = pointsToPoint(point1)  # Convert N-9 point to N-3 point
    vector = point - centroid  # Distance between point and centroid
    # Return boolean whether the distance is within the radius
    return np.linalg.norm(vector) < radius

# Check if a given point is within the given radius from the centroid for Trimesh object


def withinBoundsT(point, centroid, radius):
    vector = point - centroid  # Distance between point and centroid
    # Return boolean whether the distance is within the radius
    return np.linalg.norm(vector) < radius

# Trim a circle around the centroid of the mesh


def trimCircle(mesh, radius):
    meshObjTrimesh = trimesh.Trimesh(
        **trimesh.triangles.to_kwargs(mesh.vectors))  # Convert to stl.Mesh object to Trimesh object
    # Convert numpy array to python standard array
    meshPointArr = mesh.points.tolist()
    for i in range(len(meshPointArr)):  # Loop through the points
        if not withinBounds(meshPointArr[i], meshObjTrimesh.centroid, radius):
            # If outside bounds then reassign the point to orgin
            meshPointArr[i] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
    # Update the vertices of the mesh to the trimmed vertices
    mesh.points = np.array(meshPointArr)
    return mesh  # Returns stl.mesh object

# Trim a circle around the centroid of the Trimesh mesh


def trimCircleT(mesh, radius):
    # Convert numpy array to python standard array
    meshPointArr = mesh.vertices
    for i in range(len(meshPointArr)):  # Loop through the points
        if not withinBoundsT(meshPointArr[i], mesh.centroid, radius):
            # If outside bounds then reassign the point to orgin
            meshPointArr[i] = trimesh.caching.tracked_array(
                [0, 0, 0], dtype=None)
    # Update the vertices of the mesh to the trimmed vertices
    mesh.vertices = meshPointArr
    # mesh.remove_infinite_values()``
    return mesh

# Convert N-9 point to N-3 point


def trimCircleGivenPoint(mesh, point, radius):  # Works with stl.mesh object
    meshPointArr = mesh.points.tolist()
    for i in range(len(meshPointArr)):  # Loop through the points
        if not withinBounds(meshPointArr[i], point, radius):
            # If outside bounds then reassign the point to orgin
            meshPointArr[i] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
    # Update the vertices of the mesh to the trimmed vertices
    mesh.points = np.array(meshPointArr)
    return mesh


def pointsToPoint(points):
    x = (points[0] + points[3] + points[6])/3
    y = (points[1] + points[4] + points[7])/3
    z = (points[2] + points[5] + points[8])/3
    return np.array([x, y, z])

# Convert a trianglular face into a one dimensional point


def triangleToPoint(trianle):
    x = trianle[0][0] + trianle[1][0] + trianle[2][0]
    y = trianle[0][1] + trianle[1][1] + trianle[2][1]
    z = trianle[0][2] + trianle[1][2] + trianle[2][2]
    return np.array([x, y, z])

# Input: unit vector array atribute of stl.Mesh object; Output: normal vector to the flat plane


def unitVectorsToZplane(vectors):
    # Initialize the empty array for z components of the unit normal vectors
    x_comp = []
    y_comp = []
    z_comp = []

    # Add all the z components to the empty array
    for vector in vectors:
        x_comp.append(vector[0])
        y_comp.append(vector[1])
        z_comp.append(vector[2])

    x_comp = np.array(x_comp)
    y_comp = np.array(y_comp)
    z_comp = np.array(z_comp)

    x = stats.mean(x_comp)
    y = stats.mean(y_comp)
    z = stats.mean(z_comp)

    return [x, y, z]


# def lowest_point(mesh, plane):
#     points = mesh.vertices
#     min = points[0][plane]
#     z_points = []

#     # Add all the z components to the empty array
#     for point in points:
#         z_points.append(point[2])
#         if point[2] < min:
#             minPoint = point
#             min = point
#     return point

def lowest_point(mesh, plane):
    low = mesh.vertices[0]
    for m in mesh.vertices:
        if m[plane] < low[plane]:
            low = m
            # m.visual.vertex_colors = trimesh.visual.random_color()
    return low


def plane():
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


def vectorToPlane(vector):
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


def planeVector(mesh):
    nV = unitVectorsToZplane(mesh.face_normals)
    A = np.array([[1, 0, 0],
                  [np.cos(np.pi/2), -np.sin(np.pi/2), 0],
                  [np.sin(np.pi/2), np.cos(np.pi/2), 0]])
    return A @ nV


def rotationMatrix(vector1, vector2):
    cross = np.cross(vector1, vector2)
    u = cross/np.linalg.norm(cross)  # Rotation Axis
    # print("u = " + str(u))

    dot = np.dot(vector1, vector2)
    rotation_angle = math.atan2(np.linalg.norm(cross), dot)  # Rotation Angle
    # print("rotation_angle = " + str(rotation_angle))

    cos = np.cos(rotation_angle)
    cos_1 = 1 - cos
    sin = np.sin(rotation_angle)

    R = np.array([[cos+(u[0]**2)*(cos_1), u[0]*u[1]*cos_1-u[2]*sin, u[0]*u[2]*cos_1+u[1]*sin, 0],
                  [u[0]*u[1]*cos_1+u[2]*sin, cos +
                      (u[1]**2)*(cos_1), u[1]*u[2]*cos_1-u[0]*sin, 0],
                  [u[2]*u[0]*cos_1-u[1]*sin, u[2]*u[1]*cos_1 +
                      u[0]*sin, cos+(u[2]**2)*(cos_1), 0],
                  [0, 0, 0, 1]])

    # R = np.array([[0, 0, 2, 0],
    #              [1, 2, 0, 0],
    #              [3, 4, 2, 0],
    #              [0, 0, 0, 1]])
    # angle, direction, point = trimesh.transformations.rotation_from_matrix(R)
    # return trimesh.transformations.rotation_matrix(angle, direction, point)
    return R


def move(trimeshObj, distance_vector):
    trimeshObj.vertices += distance_vector
