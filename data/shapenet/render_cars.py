#!/usr/bin/env python
"""
Data rendering of generated cuboids.
Based on the configuration files in `config/`.
"""

import numpy as np
import math
import random
import copy
import pyrender
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import mayavi.mlab
import re

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import  utils
from mesh import Mesh

def plot_depth_map(depth_map):
    """
    Show the depth map as image.

    :param depth_map: depth map
    :type depth_map: numpy.ndarray
    """

    plt.imshow(depth_map, cmap = 'viridis', interpolation = 'nearest')
    plt.colorbar()
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    for config_file in os.listdir(config_folder):
        config = utils.read_json(config_folder + config_file)

        image_width = config['image_width']
        image_height = config['image_height']

        intrinsics = np.array([config['focal_length_x'], config['focal_length_y'], config['principal_point_x'], config['principal_point_y']], dtype = np.float64)
        size = np.array([image_height, image_width], dtype = np.int32)
        mesh_center = (config['mesh_center_x'], config['mesh_center_y'], config['mesh_center_z'])
        znf = np.array([config['z_near'], config['z_far']])

        height = config['height']
        width = config['width']
        depth = config['depth']
        suffix = config['suffix']

        off_dir = config['off_dir'] + '_' + str(config['multiplier']) + '_' \
                  + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                  + str(height) + 'x' + str(width) + 'x' + str(depth) + suffix

        depth_file = config['depth_file'] + '_' + str(config['multiplier']) + '_' \
                     + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                     + str(height) + 'x' + str(width) + 'x' + str(depth) + suffix + '.h5'

        angles_file = config['render_orientation_file'] + '_' + str(config['multiplier']) + '_' \
                      + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                      + str(height) + 'x' + str(width) + 'x' + str(depth) + suffix + '.h5'

        print('[Data] checking ' + off_dir)
        files = [off_dir + '/' + filename for filename in os.listdir(off_dir)]
        print('[Data] found ' + str(len(files)) + ' off files')

        # http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically
        def tryint(s):
            try:
                return int(s)
            except:
                return s


        def alphanum_key(s):
            """ Turn a string into a list of string and number chunks.
                "z23a" -> ["z", 23, "a"]
            """
            return [tryint(c) for c in re.split('([0-9]+)', s)]

        files.sort(key = alphanum_key)
        n_files = len(files)

        depth_maps = np.zeros((n_files, image_height, image_width))
        rotations = np.zeros((n_files, 1))

        for i in range(n_files):

            mesh = Mesh.from_off(files[i])

            rotations[i] = (np.random.random() * 2 * config['render_max_y_rotation'] - config['render_max_y_rotation']) / 180.0 * math.pi
            rotation = [
                0,
                rotations[i],
                0
            ]

            mesh.rotate(rotation)
            mesh.translate(mesh_center)

            np_vertices = mesh.vertices.astype(np.float64)
            np_faces = mesh.faces.astype(np.float64)
            np_faces += 1

            print('[Data] rendering %s %d/%d' % (files[i], (i + 1), n_files))
            depth_map, mask, img = pyrender.render(np_vertices.T.copy(), np_faces.T.copy(), intrinsics, znf, size)
            depth_maps[i] = depth_map

        utils.write_hdf5(angles_file, np.array(rotations))
        print('[Data] wrote ' + angles_file)

        utils.write_hdf5(depth_file, depth_maps)
        print('[Data] wrote ' + depth_file)

        #depth_maps = utils.read_hdf5(depth_file)
        #for i in range(n_files):
        #    plot_depth_map(depth_maps[i])
