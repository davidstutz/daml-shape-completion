#!/usr/bin/env python
"""
Some visualization utilities.
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as mplt
import scipy
import point_cloud

def plot_mesh(mesh, filepath = ''):
    """
    Plot a mesh.

    :param mesh: mesh to plot
    :type mesh: mesh.Mesh
    :param filepath: path to file to save the plot; if empty, the plot is shown
    :type filepath: str
    """

    ax = mplt.Axes3D(plt.figure())
    for i in range(mesh.faces.shape[0]):
        vertex = np.zeros((3, 3))

        for j in range(3):
            vertex[j, :] = mesh.vertices[mesh.faces[i, j], :]

        triangle = mplt.art3d.Poly3DCollection([vertex])
        ax.add_collection3d(triangle)

    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

def plot_volume(volume, filepath = '', step = 1):
    """
    Plot a volume as scatter plot.

    :param volume_a: volume
    :type volume_a: numpy.ndarray
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param step: take every step-th point only
    :type step: int
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    yy, xx, zz = np.where(volume > 0.5)
    xx = xx[::step]
    yy = yy[::step]
    zz = zz[::step]

    ax.scatter(xx, yy, zz, c = zz)

    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

def plot_point_cloud(points, filepath = '', step = 1):
    """
    Plot a point cloud using the given points.

    :param points: N x 3 point matrix
    :type points: numpy.ndarray
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param step: take every step-th point only
    :type step: int
    """

    if isinstance(points, point_cloud.PointCloud):
        points = points.points

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    xx = points[::step, 0]
    yy = points[::step, 1]
    zz = points[::step, 2]

    ax.scatter(xx, yy, zz, c=zz, s=1)

    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()

def plot_point_clouds(point_clouds, filepath = ''):
    assert len(point_clouds) > 0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    c = 0
    for points in point_clouds:
        xx = points[:, 0]
        yy = points[:, 1]
        zz = points[:, 2]

        ax.scatter(xx, yy, zz, c = 0)
        c = c + 1

    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()

def plot_point_cloud_error(point_clouds, filepath = ''):
    assert len(point_clouds) == 2

    points_a = point_clouds[0]
    points_b = point_clouds[1]

    distances = np.zeros((points_a.shape[0], points_b.shape[0]))
    for n in range(points_a.shape[0]):
        points = np.repeat(points_a[n, :].reshape((1, 3)), points_b.shape[0], axis = 0)
        distances[n, :] = np.sum(np.square(points - points_b), axis = 1).transpose()

    min_indices = np.argmin(distances, axis = 1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for n in range(points_a.shape[0]):
        ax.plot(np.array([points_a[n, 0], points_b[min_indices[n], 0]]),
                np.array([points_a[n, 1], points_b[min_indices[n], 1]]),
                np.array([points_a[n, 2], points_b[min_indices[n], 2]]))

    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()

def plot_volumes(volume_a, volume_b, filepath = '', step_a = 1, step_b = 1):
    """
    Plot two volumes in one point cloud using different colors.

    :param volume_a: first volume
    :type volume_a: numpy.ndarray
    :param volume_b: second volume
    :type volume_b: numpy.ndarray
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param step_a: take every step-th point of volume a only
    :type step_a: int
    :param step_b: take every step-th point of volume a only
    :type step_b: int
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    yy_o, xx_o, zz_o = np.where(volume_a > 0.5)
    xx_o = xx_o[::step_a]
    yy_o = yy_o[::step_a]
    zz_o = zz_o[::step_a]

    yy_i, xx_i, zz_i = np.where(volume_b > 0.5)
    xx_i = xx_i[::step_b]
    yy_i = yy_i[::step_b]
    zz_i = zz_i[::step_b]

    c_o = np.ones((xx_o.shape[0], 3))
    c_o[:, 0] = 0
    c_i = np.ones((xx_i.shape[0], 3))
    c_i[:, 1] = 0

    c = np.concatenate((c_o, c_i), axis = 0)
    xx = np.concatenate((xx_o, xx_i), axis = 0)
    yy = np.concatenate((yy_o, yy_i), axis = 0)
    zz = np.concatenate((zz_o, zz_i), axis = 0)

    ax.scatter(xx, yy, zz, c = c)

    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

def plot_volume_angles(volume, angles = [30, 60, 90, 120, 150, 180], filepath = '', step = 1):
    """
    Plots the volume from different angle in several subplots.

    :param volume: volume to plot
    :type volume: numpy.ndarray
    :param anles: list of angles in degrees
    :type angles: [float]
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param step: take every step-th point only
    :type step: int
    """
    n_angles = len(angles)
    assert n_angles > 0

    fig = plt.figure()

    yy, xx, zz = np.where(volume > 0.5)
    xx = xx[::step]
    yy = yy[::step]
    zz = zz[::step]

    elevations = [-30, 0, 30]
    k = 1
    for j in range(3):
        for i in range(n_angles):
            ax = fig.add_subplot(3, n_angles, k, projection = '3d')
            k = k + 1

            ax.scatter(xx, yy, zz, c = zz)
            ax.view_init(elevations[j], angles[i])


    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

def plot_volumes_angles(volume_a, volume_b, angles = [30, 60, 90, 120, 150, 180], filepath = '', step_a = 1, step_b = 1):
    """
    Plots two volumes in the same plot but from different angles in several subplots.

    :param volume_a: first volume
    :type volume_a: numpy.ndarray
    :param volume_b: second volume
    :type: volume_b: numpy.ndarray
    :param angles: list of angles in degrees
    :type angles: [float]
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param step_a: show only every step_a-th point in the first volume
    :type step_a: int
    :param step_b: show only everay step_b-th point in the second volume
    :type step_b: int
    """

    n_angles = len(angles)
    assert n_angles > 0

    fig = plt.figure()

    yy_o, xx_o, zz_o = np.where(volume_a > 0.5)
    xx_o = xx_o[::step_a]
    yy_o = yy_o[::step_a]
    zz_o = zz_o[::step_a]

    yy_i, xx_i, zz_i = np.where(volume_b > 0.5)
    xx_i = xx_i[::step_b]
    yy_i = yy_i[::step_b]
    zz_i = zz_i[::step_b]

    c_o = np.ones((xx_o.shape[0], 3))
    c_o[:, 0] = 0
    c_i = np.ones((xx_i.shape[0], 3))
    c_i[:, 1] = 0

    c = np.concatenate((c_o, c_i), axis = 0)
    xx = np.concatenate((xx_o, xx_i), axis = 0)
    yy = np.concatenate((yy_o, yy_i), axis = 0)
    zz = np.concatenate((zz_o, zz_i), axis = 0)

    elevations = [-30, 0, 30]
    k = 1
    for j in range(3):
        for i in range(n_angles):
            ax = fig.add_subplot(3, n_angles, k, projection = '3d')
            k = k + 1

            ax.scatter(xx, yy, zz, c = c)
            ax.view_init(elevations[j], angles[i])

    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()


def plot_slices(volumes, filepath = '', axis = 0, vmin = 0, vmax = 1):
    """
    Plot two volumes using several planes.

    :param volumes: list of volumes
    :type volumes: [numpy.ndarray]
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param axis: which axis to slice in
    :type axis: int
    :param vmin: minimum value
    :type vmin: float
    :param vmax: maximum value
    :type vmax: float
    """

    if not isinstance(volumes, list):
        assert isinstance(volumes, np.ndarray)
        volumes = [volumes]

    n_volumes = len(volumes)
    assert n_volumes > 0

    axis = int(axis)
    assert axis >= 0 and axis <= 3

    n_slices = volumes[0].shape[axis] // 2
    fig, ax = plt.subplots(n_volumes, n_slices)
    fig.set_size_inches(20, 20)

    for j in range(n_slices):
        for i in range(n_volumes):
            if n_volumes <= 1:
                a = ax[j]
            else:
                a = ax[i, j]

            if axis == 0:
                im = a.imshow(volumes[i][j * volumes[i].shape[axis] // n_slices, :, :], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)
            elif axis == 1:
                im = a.imshow(volumes[i][:, j * volumes[i].shape[axis] // n_slices, :], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)
            elif axis == 2:
                im = a.imshow(volumes[i][:, :, j * volumes[i].shape[axis] // n_slices], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)

    #fig.colorbar(im, ax = ax.ravel().tolist())

    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

def plot_specific_slices(volumes, slices, filepath = '', axis = 0, vmin = 0, vmax = 1):
    """
    Plot two volumes using several planes.

    :param volumes: list of volumes
    :type volumes: [numpy.ndarray]
    :param filepath: path to file to save plot to; plot is shown if empty
    :type filepath: str
    :param axis: which axis to slice in
    :type axis: int
    :param vmin: minimum value
    :type vmin: float
    :param vmax: maximum value
    :type vmax: float
    """

    if not isinstance(volumes, list):
        assert isinstance(volumes, np.ndarray)
        volumes = [volumes]

    n_volumes = len(volumes)
    assert n_volumes > 0

    axis = int(axis)
    assert axis >= 0 and axis <= 3

    assert len(slices) > 0

    fig, ax = plt.subplots(n_volumes, len(slices))
    fig.set_size_inches(20, 20)

    for j in range(len(slices)):
        for i in range(n_volumes):
            if n_volumes <= 1:
                a = ax[j]
            else:
                a = ax[i, j]

            if axis == 0:
                im = a.imshow(volumes[i][slices[j], :, :], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)
            elif axis == 1:
                im = a.imshow(volumes[i][:, slices[j], :], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)
            elif axis == 2:
                im = a.imshow(volumes[i][:, :, slices[j]], cmap = 'coolwarm', interpolation = 'nearest',
                              vmin = vmin, vmax = vmax)

    fig.colorbar(im, ax = ax.ravel().tolist())

    if filepath:
        plt.savefig(filepath, bbox_inches = 'tight')
    else:
        plt.show()

try:
    import mayavi.mlab
    import imageio

    def plot_volume_mayavi(volume, filepath = ''):
        """
        Plot a volume using mayavi.

        :param volume: volume
        :type volume: numpy.ndarray
        :param filepath: path to file to save the figure to; plot is shown if empty
        :type filepath: str
        """

        yy, xx, zz = np.where(volume > 0.5)
        mayavi.mlab.points3d(xx, yy, zz, mode = 'cube', color = (0, 1, 0), scale_factor = 1)

        if filepath:
            mayavi.mlab.savefig(filepath)
        else:
            mayavi.mlab.show()

    def plot_volume_angles_mayavi(volume, angles = [15, 45, 75, 105, 135, 165], filename = ''):
        """
        Plots two volumes in the same plot but from different angles in several subplots.

        :param volume: volume
        :type volume: numpy.ndarray
        :param angles: list of angles in degrees
        :type angles: [float]
        :param filename:filename with path but without extension, multiple images are created; plot is shown if empty
        :type filename: str
        """

        yy, xx, zz = np.where(volume > 0.5)
        height = volume.shape[0]
        width = volume.shape[1]
        depth = volume.shape[2]

        fig = mayavi.mlab.figure(bgcolor = (1, 1, 1))

        filenames = []
        for azimuth in angles:
            mayavi.mlab.points3d(xx, yy, zz, extent = [0, height, 0, width, 0, depth], mode='cube', color=(0,83./255.,159./255.), scale_factor=1)
            mayavi.mlab.view(0, azimuth, 'auto', 'auto')
            #mayavi.mlab.roll(azimuth)

            if filename:
                filenames.append(filename + '_' + str(azimuth) + '.png')
                mayavi.mlab.savefig(filenames[-1])
            else:
                mayavi.mlab.show()
            mayavi.mlab.clf()

        return filenames

    def combine_plots_mayavi(filenames, filepath, delete = False, ignore_first = False, flip_x_y = False):
        """
        Combine the plots int he given filenames in a single image.

        :param filenames: list of filenames
        :type filenames: [str]
        """

        if ignore_first:
            if delete:
                os.unlink(filenames[0])
            filenames = filenames[1:]

        assert len(filenames) > 0
        image = scipy.misc.imread(filenames[0])
        height = image.shape[0]
        width = image.shape[1]
        channels = image.shape[2]

        combined = np.zeros((height, len(filenames)*width, channels), dtype = image.dtype)

        i = 0
        for filename in filenames:
            image = scipy.misc.imread(filename)

            if flip_x_y:
                image = image.transpose((1, 0, 2))

            combined[:, i*width: (i + 1)*width, :] = image

            if delete:
                os.unlink(filename)
            i += 1

        scipy.misc.imsave(filepath, combined)

    def combine_plots_mayavi_gif(filenames, filepath, delete = False, ignore_first = False, flip_x_y = False):
        """
        Combine the plots in the given filenames in a single image.

        :param filenames: list of filenames
        :type filenames: [str]
        """

        if ignore_first:
            if delete:
                os.unlink(filenames[0])

            filenames = filenames[1:]

        assert len(filenames) > 0

        images = []
        for filename in filenames:
            image = imageio.imread(filename)

            if flip_x_y:
                image = image.transpose((1, 0, 2))

            images.append(image)

            if delete:
                os.unlink(filename)

        imageio.mimsave(filepath, images)
except:
    pass