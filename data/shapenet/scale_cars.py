#!/usr/bin/env python
"""
Scale cars.
"""

import numpy as np
from scipy import ndimage
import math

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import  utils
from mesh import Mesh

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    config_files = [config_file for config_file in os.listdir(config_folder)]
    config = utils.read_json(config_folder + config_files[-1])

    simplified_directory = config['simplified_directory'] + '/'
    scaled_directory = config['model_directory'] + '/'

    if not os.path.exists(scaled_directory):
        os.makedirs(scaled_directory)

    filter_models = config['filter_models']
    files = {int(file[:-4]): file for file in os.listdir(simplified_directory) if file[-4:] == '.off'}
    keys = sorted(files.keys(), key=int)

    for key in keys:
        if not key in filter_models:
            in_file = simplified_directory + files[key]
            mesh = Mesh.from_off(in_file)

            # get extents in order to do scaling
            # and translation to (0,0,0)
            min, max = mesh.extents()

            total_min = np.min(np.array(min))
            total_max = np.max(np.array(max))

            centers = ((min[0] + max[0]) / 2, (min[1] + max[1]) / 2, (min[2] + max[2]) / 2)
            sizes = (total_max - total_min, total_max - total_min, total_max - total_min)
            translation = (-centers[0], -centers[1], -centers[2])
            scales = (1 / (sizes[0] + config['padding'] * sizes[0]), 1 / (sizes[1] + config['padding'] * sizes[1]),
                      1 / (sizes[2] + config['padding'] * sizes[2]))

            mesh.translate(translation)
            mesh.scale(scales)
            mesh.rotate((0, config['y_rotation'] / 180. * math.pi, 0))

            print('[Data] %s extents before %f - %f, %f - %f, %f - %f' % (key, min[0], max[0], min[1], max[1], min[2], max[2]))
            min, max = mesh.extents()
            print('[Data] %s extents after %f - %f, %f - %f, %f - %f' % (key, min[0], max[0], min[1], max[1], min[2], max[2]))

            out_file = scaled_directory + '/%04d.off' % key
            mesh.to_off(out_file)