#!/usr/bin/env python
"""
Some sanity checks of the rendering.
"""

from matplotlib import pyplot as plt
import math
import random
import numpy as np
import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
import visualization
import ntpath
import bounding_box
import common

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('Usage!')

    config_file = sys.argv[1]
    assert os.path.exists(config_file)

    print('[Data] reading ' + config_file)
    config = utils.read_json(config_file)
    dataset = ntpath.basename(config_file)[:-5]

    inputs = utils.read_hdf5(common.filename(config, 'input_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_file', '_f.h5', dataset))
    #inputs_sdf = utils.read_hdf5(common.filename(config, 'input_sdf_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'input_sdf_file', '_f.h5', dataset))

    #inputs_combined_gt = utils.read_hdf5(common.filename(config, 'input_combined_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'input_combined_gt_file', '_f.h5', dataset))
    #inputs_sdf_combined_gt = utils.read_hdf5(common.filename(config, 'input_sdf_combined_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'input_sdf_combined_gt_file', '_f.h5', dataset))

    #inputs_gt = utils.read_hdf5(common.filename(config, 'input_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'input_gt_file', '_f.h5', dataset))
    #inputs_sdf_gt = utils.read_hdf5(common.filename(config, 'input_sdf_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'input_sdf_gt_file', '_f.h5', dataset))
    #full_space_gt = utils.read_hdf5(common.filename(config, 'full_space_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'full_space_gt_file', '_f.h5', dataset))
    #part_space_gt = utils.read_hdf5(common.filename(config, 'part_space_gt_file', '_f.h5', dataset))
    #print('[Data] read ' + common.filename(config, 'part_space_gt_file', '_f.h5', dataset))

    full_space = utils.read_hdf5(common.filename(config, 'full_space_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'full_space_file', '_f.h5', dataset))
    part_space = utils.read_hdf5(common.filename(config, 'part_space_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'part_space_file', '_f.h5', dataset))
    bounding_boxes = bounding_box.read_bounding_boxes(common.filename(config, 'bounding_box_file', '.txt', dataset))
    print('[Data] read ' + common.filename(config, 'bounding_box_file', '.txt', dataset))
    #depths = utils.read_hdf5(common.filename(config, 'depth_file'))
    #print('[Data] read ' + common.filename(config, 'depth_file'))
    #cylindral = utils.read_hdf5(common.filename(config, 'cylindral_file'))
    #print('[Data] read ' + common.filename(config, 'cylindral_file'))

    print(inputs.shape)
    print(inputs_sdf.shape)
    print(full_space.shape)
    print(part_space.shape)

    for i in range(min(25, inputs.shape[0])):
        n = i#random.randint(0, inputs.shape[0] - 1)

        print('[Data] visualizing %s %d/%d' % (dataset, (n + 1), inputs.shape[0]))

        #n_voxels = np.sum(inputs[n])
        #n_space = np.sum(space[n])

        #center = bounding_boxes[n].translation
        #distance = math.sqrt(center[0] * center[0] + center[2] * center[2])

        #if n_voxels < config['min_input_voxels'] or n_space < config['min_space_voxels'] or distance > config['max_bounding_box_distance']:
        #    print('   filtered (%d/%d %d/%d %f/%f' % (
        #        n_voxels, config['min_input_voxels'],
        #        n_space, config['min_space_voxels'],
        #        distance, config['max_bounding_box_distance']
        #    ))
        #else:
        #    print('   not filtered (%d/%d %d/%d %f/%f' % (
        #        n_voxels, config['min_input_voxels'],
        #        n_space, config['min_space_voxels'],
        #        distance, config['max_bounding_box_distance']
        #    ))

        slices = []
        for j in range(2, inputs[n][0].shape[0] - 2, 1):
            slices.append(j)

        visualization.plot_specific_slices([inputs[n][0], full_space[n][0], part_space[n][0]], slices, dataset + '_' + str(n) + '.png')
        #visualization.plot_specific_slices([inputs_combined_gt[n][0], inputs_gt[n][0], full_space_gt[n][0], part_space_gt[n][0], inputs_gt[n][1], full_space_gt[n][1], part_space_gt[n][1], inputs_gt[n][2], full_space_gt[n][2], part_space_gt[n][2]], slices, dataset + '_' + str(n) + '_gt.png')
        #visualization.plot_specific_slices([inputs_sdf_combined_gt[n][0], inputs_sdf[n][0]], slices, dataset + '_' + str(n) + '_sdf.png', 0, np.min(inputs_sdf[n]),np.max(inputs_sdf[n]))
        #visualization.plot_volumes(space[n][0], inputs[n][0], '', 10, 1)

        #visualization.plot_volume_mayavi(inputs[n][0])
        #visualization.plot_volume_mayavi(full_space[n][0])
        #visualization.plot_volume_mayavi(part_space[n][0])

        #plt.imshow(depths[n, :, :, 0])
        #plt.imshow(cylindral[n])
        #plt.show()
