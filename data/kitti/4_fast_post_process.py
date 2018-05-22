#!/usr/bin/env python
"""
Post-process computed volumes, i.e. compute the distance transforms.
"""

import numpy as np
from scipy import ndimage

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import  utils
import ntpath
import common

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_file = sys.argv[1]
    assert os.path.exists(config_file)

    print('[Data] reading ' + config_file)
    config = utils.read_json(config_file)
    set = ntpath.basename(config_file)[:-5]

    # ['input_file', 'input_sdf_file', 'full_space_file', 'part_space_file']
    size = 0
    inputs_file = common.filename(config, 'input_file', '.h5', set)
    inputs = utils.read_hdf5(inputs_file)
    print('[Data] read ' + inputs_file)

    inputs_sdf_file = common.filename(config, 'input_sdf_file', '.h5', set)
    inputs_sdf = utils.read_hdf5(inputs_sdf_file)
    print('[Data] read ' + inputs_sdf_file)

    if config['gt']:
        inputs_combined_gt_file = common.filename(config, 'input_combined_gt_file', '.h5', set)
        inputs_combined_gt = utils.read_hdf5(inputs_combined_gt_file)
        print('[Data] read ' + inputs_combined_gt_file)

        #inputs_sdf_combined_gt_file = common.filename(config, 'input_sdf_combined_gt_file', '.h5', set)
        #inputs_sdf_combined_gt = utils.read_hdf5(inputs_sdf_combined_gt_file)

        #inputs_gt_file = common.filename(config, 'input_gt_file', '.h5', set)
        #inputs_gt = utils.read_hdf5(inputs_gt_file)
        #print('[Data] read ' + inputs_gt_file)

        #inputs_sdf_gt_file = common.filename(config, 'input_sdf_gt_file', '.h5', set)
        #inputs_sdf_gt = utils.read_hdf5(inputs_sdf_gt_file)

        #full_space_gt_file = common.filename(config, 'full_space_gt_file', '.h5', set)
        #full_space_gt = utils.read_hdf5(full_space_gt_file)
        #print('[Data] read ' + full_space_gt_file)

        #part_space_gt_file = common.filename(config, 'part_space_gt_file', '.h5', set)
        #part_space_gt = utils.read_hdf5(part_space_gt_file)
        #print('[Data] read ' + part_space_gt_file)

    full_space_file = common.filename(config, 'full_space_file', '.h5', set)
    full_space = utils.read_hdf5(full_space_file)
    print('[Data] read ' + full_space_file)

    part_space_file = common.filename(config, 'part_space_file', '.h5', set)
    part_space = utils.read_hdf5(part_space_file)
    print('[Data] read ' + part_space_file)

    assert len(inputs.shape) >= 4
    if len(inputs.shape) < 5:
        inputs = np.expand_dims(inputs, axis=1)

    assert len(inputs_sdf.shape) >= 4
    if len(inputs_sdf.shape) < 5:
        inputs_sdf = np.expand_dims(inputs_sdf, axis=1)

    if config['gt']:
        assert len(inputs_combined_gt.shape) >= 4
        if len(inputs_combined_gt.shape) < 5:
            inputs_combined_gt = np.expand_dims(inputs_combined_gt, axis=1)

        #assert len(inputs_sdf_combined_gt.shape) >= 4
        #if len(inputs_sdf_combined_gt.shape) < 5:
        #    inputs_sdf_combined_gt = np.expand_dims(inputs_sdf_combined_gt, axis=1)

    assert len(full_space.shape) >= 4
    if len(full_space.shape) < 5:
        full_space = np.expand_dims(full_space, axis=1)

    assert len(part_space.shape) >= 4
    if len(part_space.shape) < 5:
        part_space = np.expand_dims(part_space, axis=1)

    observed = np.sum(np.sum(np.sum(inputs + part_space, axis=4), axis=3), axis=2)
    print(observed.shape)
    inputs_indices = np.nonzero(np.array(observed))
    inputs_indices = inputs_indices[0]
    print(inputs_indices)

    print('[Data] ' + str(inputs_indices.shape[0]) + ' non-zero observations found')
    assert inputs_indices.shape[0] > 0

    print(inputs.shape)
    print(inputs_sdf.shape)
    if config['gt']:
        print(inputs_combined_gt.shape)
        #print(inputs_sdf_combined_gt.shape)
        #print(inputs_gt.shape)
        #print(inputs_sdf_gt.shape)
        #print(full_space_gt.shape)
        #print(part_space_gt.shape)
    print(full_space.shape)
    print(part_space.shape)

    inputs = inputs[inputs_indices]
    inputs_sdf = inputs_sdf[inputs_indices]
    if config['gt']:
        inputs_combined_gt = inputs_combined_gt[inputs_indices]
        #inputs_sdf_combined_gt = inputs_sdf_combined_gt[inputs_indices]
        #inputs_gt = inputs_gt[inputs_indices]
        #inputs_sdf_gt = inputs_sdf_gt[inputs_indices]
        #full_space_gt = full_space_gt[inputs_indices]
        #part_space_gt = part_space_gt[inputs_indices]
    full_space = full_space[inputs_indices]
    part_space = part_space[inputs_indices]

    full_space[inputs == 1] = 0
    part_space[inputs == 1] = 0
    #if config['gt']:
    #    full_space_gt[inputs_gt == 1] = 0
    #    part_space_gt[inputs_gt == 1] = 0

    print(inputs.shape)
    print(inputs_sdf.shape)
    if config['gt']:
        print(inputs_combined_gt.shape)
        #print(inputs_sdf_combined_gt.shape)
        #print(inputs_gt.shape)
        #print(inputs_sdf_gt.shape)
        #print(full_space_gt.shape)
        #print(part_space_gt.shape)
    print(full_space.shape)
    print(part_space.shape)

    inputs_file = common.filename(config, 'input_file', '_f.h5', set)
    utils.write_hdf5(inputs_file, inputs)
    print('[Data] wrote ' + inputs_file)

    inputs_sdf_file = common.filename(config, 'input_sdf_file', '_f.h5', set)
    utils.write_hdf5(inputs_sdf_file, inputs_sdf)
    print('[Data] wrote ' + inputs_sdf_file)

    if config['gt']:
        inputs_combined_gt_file = common.filename(config, 'input_combined_gt_file', '_f.h5', set)
        utils.write_hdf5(inputs_combined_gt_file, inputs_combined_gt)
        print('[Data] wrote ' + inputs_combined_gt_file)

        #inputs_sdf_combined_gt_file = common.filename(config, 'input_sdf_combined_gt_file', '_f.h5', set)
        #utils.write_hdf5(inputs_sdf_combined_gt_file, inputs_sdf_combined_gt)
        #print('[Data] wrote ' + inputs_sdf_combined_gt_file)

        #inputs_gt_file = common.filename(config, 'input_gt_file', '_f.h5', set)
        #utils.write_hdf5(inputs_gt_file, inputs_gt)
        #print('[Data] wrote ' + inputs_gt_file)

        #inputs_sdf_gt_file = common.filename(config, 'input_sdf_gt_file', '_f.h5', set)
        #utils.write_hdf5(inputs_sdf_gt_file, inputs_sdf_gt)
        #print('[Data] wrote ' + inputs_sdf_gt_file)

        #full_space_gt_file = common.filename(config, 'full_space_gt_file', '_f.h5', set)
        #utils.write_hdf5(full_space_gt_file, full_space_gt)
        #print('[Data] wrote ' + full_space_gt_file)

        #part_space_gt_file = common.filename(config, 'part_space_gt_file', '_f.h5', set)
        #utils.write_hdf5(part_space_gt_file, part_space_gt)
        #print('[Data] wrote ' + part_space_gt_file)

    full_space_file = common.filename(config, 'full_space_file', '_f.h5', set)
    utils.write_hdf5(full_space_file, full_space)
    print('[Data] wrote ' + full_space_file)

    part_space_file = common.filename(config, 'part_space_file', '_f.h5', set)
    utils.write_hdf5(part_space_file, part_space)
    print('[Data] wrote ' + part_space_file)
