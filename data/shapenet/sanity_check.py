#!/usr/bin/env python
"""
Some sanity checks of the rendering.
"""

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import numpy as np
import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
import visualization

def filename(key, level = 0):
    """
    Get the real file name by looking up the key in the config and suffixing.

    :param key: key to use in the config
    :type key: str
    :return: filepath
    :rtype: str
    """

    filename = config[key] + '_' + str(config['multiplier']) + '_' \
               + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
               + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth'])

    if level > 0:
        filename = filename + '_' + str(level)

    return filename + config['suffix'] + '.h5'

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    config_file = 'validation.json'
    config = utils.read_json(config_folder + config_file)
    dataset = config_file[:-5]

    depths = utils.read_hdf5(filename('depth_file'))
    print('[Data] read ' + filename('depth_file'))
    filled = utils.read_hdf5(filename('output_file'))
    print('[Data] read ' + filename('output_file'))

    #outputs_df = utils.read_hdf5(filename('output_df_file'))
    #print('[Data] read ' + filename('output_df_file'))
    #outputs_rtdf = utils.read_hdf5(filename('output_rtdf_file'))
    #print('[Data] read ' + filename('output_rtdf_file'))
    inputs = utils.read_hdf5(filename('input_file'))
    print('[Data] read ' + filename('input_file'))
    #inputs_df = utils.read_hdf5(filename('input_df_file'))
    #print('[Data] read ' + filename('input_df_file'))
    #inputs_rtdf = utils.read_hdf5(filename('input_rtdf_file'))
    #print('[Data] read ' + filename('input_rtdf_file'))
    space = utils.read_hdf5(filename('space_file'))
    print('[Data] read ' + filename('space_file'))
    #space_df = utils.read_hdf5(filename('space_df_file'))
    #print('[Data] read ' + filename('space_df_file'))
    #space_rtdf = utils.read_hdf5(filename('space_rtdf_file'))
    #print('[Data] read ' + filename('space_rtdf_file'))

    real_sdf = utils.read_hdf5(filename('sdf_file'))
    input_sdf = utils.read_hdf5(filename('input_sdf_file'))

    #print('[Data] outputs DF: %s' % ' x '.join(map(str, outputs_df.shape)))
    #print('[Data] outputs RTDF: %s' % ' x '.join(map(str, outputs_rtdf.shape)))
    print('[Data] inputs: %s' % ' x '.join(map(str, inputs.shape)))
    #print('[Data] inputs DF: %s' % ' x '.join(map(str, inputs_df.shape)))
    #print('[Data] inputs RTDF: %s' % ' x '.join(map(str, inputs_rtdf.shape)))
    print('[Data] space: %s' % ' x '.join(map(str, space.shape)))
    #print('[Data] space DF: %s' % ' x '.join(map(str, space_df.shape)))
    #print('[Data] space RTDF: %s' % ' x '.join(map(str, space_rtdf.shape)))
    print('[Data] real sdf: %s' % ' x '.join(map(str, real_sdf.shape)))
    print('[Data] input sdf: %s' % ' x '.join(map(str, input_sdf.shape)))

    #sanity_dir = dataset + '_sanity_' \
    #        + '_' + str(height) + '_' + str(width) + '_' + str(depth) + '/'

    #if not os.path.exists(sanity_dir):
    #    os.makedirs(sanity_dir)

    for i in range(min(25, filled.shape[0])):
        #n = random.randint(0, outputs[-1].shape[0])
        n = i
        print('[Data] visualizing %s %d/%d' % (dataset, (n + 1), filled.shape[0]))

        #volumes = [labeled[n][0]]
        #for level in range(config['levels'] + 1):
        #    volumes.append(outputs[level][n][0])
        #    volumes.append(filled[level][n][0])

        #visualization.plot_slices(volumes, dataset + '_' + str(n) + '.png', 0, 0, np.max(labeled))

        #plt.clf()
        #plt.imshow(depths[n])
        #plt.savefig(dataset + '_' + str(n) + '_depth.png')

        slices = []
        for j in range(2, filled[n][0].shape[0] - 2, 1):
            slices.append(j)

        filled_inputs = filled[n][0].copy()
        filled_inputs[inputs[n][0] == 1] = 2
        space_inputs = space[n][0].copy()
        space_inputs[inputs[n][0] == 1] = 2
        mask = np.logical_and(inputs[n][0] == 1, space[n][0] == 1)
        space_inputs[mask] = 3

        real_sdf_bin = real_sdf[n][0].copy()
        real_sdf_bin[real_sdf_bin >= 0] = 0
        real_sdf_bin[real_sdf_bin < 0] = 1
        visualization.plot_specific_slices([filled[n][0], inputs[n][0], space[n][0], filled_inputs, space_inputs, real_sdf_bin], slices, dataset + '_' + str(n) + '.png', 0, 0, 3)

        #visualization.plot_specific_slices([real_sdf[n][0], input_sdf[n][0]], slices, dataset + '_' + str(n) + '_sdf.png', 0, min(np.min(real_sdf[n][0]), np.min(input_sdf[n][0])), max(np.max(real_sdf[n][0]), np.max(input_sdf[n][0])))
