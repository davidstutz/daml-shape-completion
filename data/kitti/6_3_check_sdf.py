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
    inputs_sdf = utils.read_hdf5(common.filename(config, 'input_sdf_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_sdf_file', '_f.h5', dataset))
    inputs_tsdf = utils.read_hdf5(common.filename(config, 'input_tsdf_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_tsdf_file', '_f.h5', dataset))
    inputs_ltsdf = utils.read_hdf5(common.filename(config, 'input_ltsdf_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_ltsdf_file', '_f.h5', dataset))

    for i in range(min(25, inputs.shape[0])):
        n = i#random.randint(0, inputs.shape[0] - 1)

        print('[Data] visualizing %s %d/%d' % (dataset, (n + 1), inputs.shape[0]))

        slices = []
        for j in range(2, inputs[n][0].shape[0] - 2, 1):
            slices.append(j)

        visualization.plot_specific_slices([inputs[n][0]*np.max(inputs_sdf[n][0]), inputs_sdf[n][0], inputs_tsdf[n][0], inputs_ltsdf[n][0]], slices, dataset + '_' + str(n) + '.png', 0, 0, np.max(inputs_sdf[n][0]))
