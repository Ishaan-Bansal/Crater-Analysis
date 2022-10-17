# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:23:00 2021

@author: Nicolas Rasmont
"""

import trimesh
#from trimesh.voxel import creation
from trimesh.proximity import signed_distance
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
#import pyacvd
#import pyvista
from skimage import measure
import skimage
import scipy


def main():

    # Create a new plot
    # figure1 = plt.figure(1)
    # axes1 = mplot3d.Axes3D(figure1)

    figure1 = plt.figure()
    axes1 = figure1.add_subplot(111, projection='3d')

    # Load the STL files and add the vectors to the plot
    capsule_mesh = trimesh.load_mesh(
        r'C:\Users\ishaa\Desktop\Research Work\Crater_STL_Files\sphere.stl')
    payload_mesh = trimesh.load_mesh(
        r'C:\Users\ishaa\Desktop\Research Work\Crater_STL_Files\cylinder_4.stl')

    cog_payload = np.array(payload_mesh.center_mass).reshape(1, 3)
    mass_payload = np.array([[50000]])

    density_backshell = np.array([[10]])
    density_heatshield = np.array([[50]])

    """MESH FRONT IS TOWARD Z"""
    payload_mesh.apply_transform(
        trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
    #capsule_mesh.apply_transform(trimesh.transformations.rotation_matrix(-np.pi/2, [0,-1,0]))

    # cog = capsule_mesh.center_mass.reshape((1,3))+np.array([[1,0,+2]]) #offset a bit to the fron
    # mass = 30000 #kg

    x0 = np.array([-5, +5, -5, 5, -1, +1, 5, 5, 10, 10])
    bnds = ((-10, 10), (-10, 10), (-10, 10), (-10, 10), (-10, 10),
            (-10, 10), (2, 10), (2, 10), (2, 20), (2, 20))

    #x0 = np.array([0,0,0,5,10])
    #bnds = ((-7,7),(-7,7),(-7,7),(2,10),(2,20))

    # minimized_vect = scipy.optimize.minimize(get_score_from_rbf, x0,
    #                                           args=(payload_mesh,cog_payload,
    #                                                 mass_payload,density_backshell,
    #                                                 density_heatshield),
    #                                           method='Powell', bounds=bnds, tol=None,
    #                                           callback=None,options={'disp': True,'return_all': False})

    # output_vect = minimized_vect.x.flatten()
    # output_vect = np.array([ 1.74982742, 9.99995965, 1.03545331, 9.99995965,
    #                         0.34543452, 9.99995965, 5.11639525, 9.99995617,
    #                         13.7406775, 19.99993369])
    #output_vect = np.array([-0.20101834, -0.83980429,  1.11754894,  6.35887119,  7.46984495])
   # output_vect = np.array([-4.87713703,  1.36872484,  2.26297287,  0.89099608, -6.43089061,  5.27864045,
   # 5.,          5.,         10.,         10.        ])
    output_vect = x0
#    output_vect = np.array([ 1.97786848, 2.30247107, 0.68347034, 5.38309203, 13.66905182])
    #output_vect = np.array([3.29382748, 1.76493875, 1.1228608, 5.35129081, 17.10268396])
#    output_vect = np.array([ 0.6121862, 1.14006256, 2.22089314, 5.54923182, 14.71090453])
#    output_vect = np.array([ 1.20125411,  1.00748952, -2.95322289,  1.69408634, -4.17304105,  5.1103792,
#   6.36664916,  4.47994171,  6.87054156,  3.6243651 ])
    print(output_vect, '\n')

    score = get_score_from_rbf(output_vect, payload_mesh, cog_payload,
                               mass_payload, density_backshell, density_heatshield)
    capsule_mesh = get_mesh_from_grbf(output_vect[0:2], output_vect[2:4], output_vect[4:6], output_vect[6:8],
                                      output_vect[8:10], (-15, 15), (-15, 15), (-15, 15), 15)

    ld_ratio, ballistic_coef, phi_trim, theta_trim, A_ref, Cd, Cl, dCm_long_dtheta, total_mass, total_cog = get_aero_properties(
        capsule_mesh, cog_payload, mass_payload, density_backshell, density_heatshield, 20)

    print('total mass = %.3f' % total_mass, '\n', 'ld_ratio = %.3f' % ld_ratio[0], '\n', 'ballistic_coef = %.3f' % ballistic_coef[0], ' kg/m2\n',
          'A_ref = %.3f' % A_ref, ' m2\n', 'Cd = %.3f' % Cd[0], '\n', 'Cl = %.3f' % Cl[0], '\n', 'dCm_dt = %.3f' % dCm_long_dtheta, '\n', 'score = %.3f' % score)

    velocity_vector = np.array([[np.sin(theta_trim)*np.cos(phi_trim),
                                 np.sin(theta_trim)*np.sin(phi_trim),
                                 np.cos(theta_trim)]])

    force_pos, force_mag = get_newton_force_distrib(
        capsule_mesh, velocity_vector, 1)

    force_mag = force_mag/np.amax(np.linalg.norm(force_mag, 2, 1))*5

    axes1.quiver(force_pos[:, 0]-1*force_mag[:, 0], force_pos[:, 1]-1*force_mag[:, 1],
                 force_pos[:, 2]-1*force_mag[:, 2], force_mag[:, 0],
                 force_mag[:, 1], force_mag[:, 2])
    axes1.quiver(total_cog[:, 0], total_cog[:, 1], total_cog[:, 2],
                 velocity_vector[:, 0]*20,
                 velocity_vector[:, 1]*20,
                 velocity_vector[:, 2]*20, color='red')
    axes1.add_collection3d(mplot3d.art3d.Poly3DCollection(capsule_mesh.triangles,
                                                          facecolor=np.array(
                                                              [(1., 1., 1., 1)]),
                                                          edgecolor=np.array([(0., 0., 0., 1)])))
    capsule_mesh.show()

    axes1.add_collection3d(mplot3d.art3d.Poly3DCollection(payload_mesh.triangles,
                                                          facecolor=np.array(
                                                              [(0., 1., 0., 0.2)]),
                                                          edgecolor=np.array([(0., 0., 0., 0.2)])))
    payload_mesh.show()

    # axes1.scatter(total_cog[0,0], total_cog[0,1], total_cog[0,2],c='red')

    # # Auto scale to the mesh size
    # scale = capsule_mesh.vertices.flatten()
    # axes1.auto_scale_xyz(1.5*scale,1.5*scale,1.5*scale)
    # #axes1.xlim((0,10))
    # # Show the plot to the screen
    # axes1.set_xlabel('x')
    # axes1.set_ylabel('y')
    # axes1.set_zlabel('z')
    # plt.show()
    axes1.plot_trisurf(capsule_mesh.vertices[:, 0], capsule_mesh.vertices[:, 1],
                       triangles=capsule_mesh.faces, Z=capsule_mesh.vertices[:, 2])
    plt.show()


def get_mesh_from_grbf(x_vect, y_vect, z_vect, constant_vect, weight_vect, xlim, ylim, zlim, resolution):
    """fixed threshold at 1"""
    X, Y, Z = np.meshgrid(np.linspace(xlim[0], xlim[1], resolution),
                          np.linspace(ylim[0], ylim[1], resolution),
                          np.linspace(zlim[0], zlim[1], resolution),)

    level = np.zeros(X.shape)

    for i in range(len(x_vect)):
        level += weight_vect[i]*np.exp(-1*((X-x_vect[i])**2 +
                                           (Y-y_vect[i])**2 +
                                           (Z-z_vect[i])**2)/constant_vect[i]**2)

    verts_vect, faces_vect, normals, values = skimage.measure.marching_cubes(
        level, level=1)

    '''scaling marching cube output correctly'''
    verts_vect[:, 0] = xlim[0] + verts_vect[:, 0] * \
        (xlim[1] - xlim[0])/resolution
    verts_vect[:, 1] = ylim[0] + verts_vect[:, 1] * \
        (ylim[1] - ylim[0])/resolution
    verts_vect[:, 2] = zlim[0] + verts_vect[:, 2] * \
        (zlim[1] - zlim[0])/resolution

    trimesh_mesh = trimesh.Trimesh(vertices=verts_vect, faces=faces_vect)
    trimesh.repair.fix_normals(trimesh_mesh)

    return trimesh_mesh


def get_score_from_rbf(x, mesh_payload, cog_payload, mass_payload,
                       density_backshell, density_heatshield):
    print(x)
    x = x.flatten()
    x_vect = x[0:2]
    y_vect = x[2:4]
    z_vect = x[4:6]
    constant_vect = x[6:8]
    weight_vect = x[8:10]

    mesh_shell = get_mesh_from_grbf(
        x_vect, y_vect, z_vect, constant_vect, weight_vect, (-15, 15), (-15, 15), (-15, 15), 15)
    ld_ratio, ballistic_coef, phi_trim, theta_trim, A_ref, Cd, Cl, dCm_long_dtheta, total_mass, total_cog = get_aero_properties(
        mesh_shell, cog_payload, mass_payload, density_backshell, density_heatshield, 20)

    ld_norm = 0.1
    ballistic_norm = 500
    Cm_long_norm = 10
    mass_norm = 20E3
    payload_wrap = signed_distance(mesh_shell, mesh_payload.vertices)
    payload_wrap[payload_wrap > 0] = 0
    payload_wrapping_factor = np.sum(payload_wrap)
    if dCm_long_dtheta < 0:
        dCm_long_dtheta = 0
    return np.sqrt((total_mass/mass_norm)**2+(ballistic_coef/ballistic_norm)**2+(payload_wrapping_factor)**2)/np.sqrt((ld_ratio/ld_norm)**2+(dCm_long_dtheta/Cm_long_norm)**2)


def get_aero_properties(mesh, cog_payload, mass_payload, density_backshell,
                        density_heatshield, resolution):

    phi_trim, theta_trim, dCm_long_dtheta = get_trim(mesh, cog_payload, mass_payload,
                                                     density_backshell, density_heatshield,
                                                     (0, np.pi), (-np.pi/2, np.pi/2), resolution)

    velocity_vector = np.array([[np.sin(theta_trim)*np.cos(phi_trim),
                                 np.sin(theta_trim)*np.sin(phi_trim),
                                 np.cos(theta_trim)]])

    aero_force, aero_moment, total_mass, total_cog = get_aero_fm(mesh, velocity_vector, 1, cog_payload,
                                                                 mass_payload, density_backshell, density_heatshield)
    lift, drag = get_lift_drag(aero_force, velocity_vector)

    """mesh points toward Z"""
    A_ref = trimesh.path.polygons.projected(mesh, normal=[0, 0, 1]).area

    Cd = np.linalg.norm(drag, 2, 1)/(A_ref*0.5)
    Cl = np.linalg.norm(lift, 2, 1)/(A_ref*0.5)

    ld_ratio = Cl/Cd

    ballistic_coef = total_mass/(A_ref*Cd)

    return ld_ratio, ballistic_coef, phi_trim, theta_trim, A_ref, Cd, Cl, dCm_long_dtheta, total_mass, total_cog


def get_plot_quiver(mesh, force_vect, force_loc, cog, cop):

    figure1 = plt.figure()
    axes1 = mplot3d.Axes3D(figure1)
    axes1.quiver(force_loc[:, 0]-1*force_vect[:, 0],
                 force_loc[:, 1]-1*force_vect[:, 1],
                 force_loc[:, 2]-1*force_vect[:, 2],
                 force_vect[:, 0], force_vect[:, 1],
                 force_vect[:, 2])

    axes1.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.triangles,
                                                          facecolor=np.array(
                                                              [(0., 0., 0., 0.)]),
                                                          edgecolor='k'))
    mesh.show()
    axes1.scatter(cog[0, 0], cog[0, 1], cog[0, 2], c='red')
    axes1.scatter(cop[0, 0], cop[0, 1], cop[0, 2], s=100, c='green')
    return None


def get_trim(mesh, cog_payload, mass_payload, density_backshell, density_heatshield, phi_minmax, theta_minmax, n):

    A_ref, Cd_mat, Cl_mat, Cm_long_mat, mass_mat, cog_mat_x, cog_mat_y, cog_mat_z, phi_mat, theta_mat = get_aero_coef(
        mesh, cog_payload, mass_payload, density_backshell, density_heatshield, np.linspace(
            phi_minmax[0], phi_minmax[1], n),
        np.linspace(theta_minmax[0], theta_minmax[1], n))

    '''detect sign change in matrix '''
    Cm_long_zero = get_mask_sign_change(Cm_long_mat)
#    dCm_long_dphi = get_mask_gradient(phi_mat,Cm_long_mat)
    dCm_long_dtheta = get_mask_gradient(theta_mat, Cm_long_mat)

    '''CHECK SIGN'''
#    trim_mask = Cm_long_zero & (dCm_long_dphi>0) & (dCm_long_dtheta>0)
    trim_mask = Cm_long_zero
#    if (trim_mask.any()==False):
#        raise NameError('No stable trim configuration discovered. Try different cog position or wider angular range.')

    phi_trim = phi_mat[trim_mask]
    theta_trim = theta_mat[trim_mask]
    dCm_long_dtheta_trim = dCm_long_dtheta[trim_mask]

    max_trim_index = np.argmax(dCm_long_dtheta_trim)

    return phi_trim[max_trim_index], theta_trim[max_trim_index], dCm_long_dtheta_trim[max_trim_index]


def rotation_matrix_from_vectors(vec1, vec2):

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 /
                                                      np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def get_mask_sign_change(array):

    sign_change_mask = (np.concatenate((np.diff(np.sign(array), axis=0),
                                       np.zeros((1, array.shape[1]))), axis=0) +
                        np.concatenate((np.zeros((1, array.shape[1])),
                                        np.diff(np.sign(array), axis=0)), axis=0) +
                        np.concatenate((np.diff(np.sign(array), axis=1),
                                       np.zeros((array.shape[0], 1))), axis=1) +
                        np.concatenate((np.zeros((array.shape[0], 1)),
                                        np.diff(np.sign(array), axis=1)), axis=1)) > 0
    return sign_change_mask


def get_mask_gradient(array_x, array_y):
    dir_0 = array_x[0, 1]-array_x[0, 0]
    if (dir_0 != 0):
        return np.concatenate((np.diff(array_y, axis=1)/np.diff(array_x, axis=1),
                               np.zeros((array_x.shape[0], 1))), axis=1)
    else:
        return np.concatenate((np.diff(array_y, axis=0)/np.diff(array_x, axis=0),
                               np.zeros((1, array_x.shape[1]))), axis=0)


def get_aeroshell_mass_prop(mesh, density_backshell, density_heatshield, force_distrib):
    '''regions with significant dynamic pressure are heatshield, those without are backshells'''

    threshold = np.amax(np.linalg.norm(force_distrib, 2, 1))/20

    face_density = np.zeros((force_distrib.shape[0], 1))
    face_density[np.linalg.norm(force_distrib, 2, 1)
                 >= threshold] = density_heatshield
    face_density[np.linalg.norm(force_distrib, 2, 1)
                 < threshold] = density_backshell

    mass = np.sum(np.array(np.expand_dims(
        mesh.area_faces, axis=1)*face_density), axis=0)

    cog = np.sum(np.array(np.expand_dims(mesh.area_faces, axis=1) *
                          mesh.triangles_center*face_density),
                 axis=0).reshape((1, 3))/mass

    return cog, mass


def get_aero_coef(mesh, cog_payload, mass_payload, density_backshell, density_heatshield, phi_range, theta_range):
    """projected area"""
    A_ref = trimesh.path.polygons.projected(mesh, normal=[0, 0, 1]).area
    l_ref = np.sqrt(A_ref)

    Cd_mat = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    Cl_mat = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    Cm_long_mat = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    phi_mat, theta_mat = np.meshgrid(phi_range, theta_range)
    mass_mat = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    cog_mat_x = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    cog_mat_y = np.zeros((theta_range.shape[0], phi_range.shape[0]))
    cog_mat_z = np.zeros((theta_range.shape[0], phi_range.shape[0]))

    for i in range(len(theta_range)):
        for j in range(len(phi_range)):

            V_vect = np.array([[np.sin(theta_mat[i, j])*np.cos(phi_mat[i, j]),
                                np.sin(theta_mat[i, j])*np.sin(phi_mat[i, j]),
                                np.cos(theta_mat[i, j])]])

            aero_force, aero_moment, total_mass, total_cog = get_aero_fm(
                mesh, V_vect, 1, cog_payload, mass_payload, density_backshell, density_heatshield)
            lift, drag = get_lift_drag(aero_force, V_vect)

            Cd_mat[i, j] = np.linalg.norm(drag, 2, 1)/(A_ref*0.5)
            Cl_mat[i, j] = np.linalg.norm(lift, 2, 1)/(A_ref*0.5)

            mass_mat[i, j] = total_mass
            cog_mat_x[i, j] = total_cog[:, 0]
            cog_mat_y[i, j] = total_cog[:, 1]
            cog_mat_z[i, j] = total_cog[:, 2]

            """longitudinal coef only"""
            long_rot_axis = unit_vector(
                np.cross(np.array([[0, 0, 1]]), V_vect))
            Cm_long_mat[i, j] = np.sum(
                aero_moment*long_rot_axis)/(A_ref*l_ref*0.5)

    return A_ref, Cd_mat, Cl_mat, Cm_long_mat, mass_mat, cog_mat_x, cog_mat_y, cog_mat_z, phi_mat, theta_mat


def get_aero_fm(mesh, velocity_vector, air_density, cog_payload, mass_payload, density_backshell, density_heatshield):
    """get force distribution from Newtonian theory"""
    force_pos, force_mag = get_newton_force_distrib(
        mesh, velocity_vector, air_density)

    """determine cog and mass"""
    cog_aeroshell, mass_aeroshell = get_aeroshell_mass_prop(
        mesh, density_backshell, density_heatshield, force_mag)

    total_mass = mass_aeroshell + mass_payload

    total_cog = (mass_aeroshell*cog_aeroshell +
                 mass_payload*cog_payload)/total_mass

    """lift-to-drag-ratio, no rotation"""
    aero_force = np.sum(force_mag, axis=0).reshape((1, force_mag.shape[1]))

    """aerodynamic moment at COG check direction of arm!!!!"""
    aero_moment = np.sum(np.cross((force_pos-total_cog), force_mag, axis=1),
                         axis=0).reshape((1, force_mag.shape[1]))

    return aero_force, aero_moment, total_mass, total_cog


def get_lift_drag(aero_force, velocity_vector):

    angle = angle_between(aero_force, velocity_vector)
    unit_velocity = unit_vector(velocity_vector)
    drag = np.linalg.norm(aero_force, 2, 1)*np.cos(angle)*unit_velocity
    lift = aero_force-drag

    return lift, drag


def get_newton_force_distrib(mesh, velocity_vector, air_density):

    lift_position = mesh.triangles_center
    area_vector = np.array(np.expand_dims(
        mesh.area_faces, axis=1)*mesh.face_normals)

    velocity_unit = np.tile(unit_vector(velocity_vector),
                            (area_vector.shape[0], 1))
    angle = angle_between(area_vector, velocity_unit)

    mask = angle > np.pi/2

    area_vector = area_vector - mask*area_vector

    force_vector = -0.5*2*air_density * \
        (np.linalg.norm(velocity_vector, 2, 1))**2*area_vector*(np.cos(angle))**2
    #pressure_vector = 2*mesh.normals*(np.cos(angle))**2
    return (lift_position, force_vector)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector/np.linalg.norm(vector, 2, 1).reshape((vector.shape[0], 1))


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    return np.arccos(np.clip(np.sum(unit_vector(v1)*unit_vector(v2), axis=1), -1.0, 1.0)).reshape((v1.shape[0], 1))


if __name__ == "__main__":
    main()
