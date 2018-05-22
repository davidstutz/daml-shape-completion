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

    for config_file in os.listdir(config_folder):
        print('[Data] reading ' + config_folder + config_file)
        config = utils.read_json(config_folder + config_file)
        
        height = config['height']
        width = config['width']
        depth = config['depth']

        space_file = filename('space_file')
        space = utils.read_hdf5(space_file)
        input_file = filename('input_file')
        input = utils.read_hdf5(input_file)

        space[input == 1] = 0
        if len(space.shape) < 5:
            space = np.expand_dims(space, axis=1)

        utils.write_hdf5(space_file, space)
        print('[Data] wrote ' + space_file)

        for key in ['input', 'space', 'output', 'sdf', 'input_sdf']:
            file = filename(key + '_file')
            volumes = utils.read_hdf5(file)

            volumes = np.squeeze(volumes)
            if len(volumes.shape) < 4:
                volumes = np.expand_dims(volumes, axis = 0)

            if len(volumes.shape) < 5:
                utils.write_hdf5(file, np.expand_dims(volumes, axis = 1))
                print('[Data] wrote ' + file)
