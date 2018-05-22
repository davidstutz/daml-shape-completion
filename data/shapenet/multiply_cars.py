#!/usr/bin/env python
"""
Pre-process modelnet.
"""

import math
import numpy as np
from scipy import ndimage

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import  utils
import mesh
import random

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    for config_file in os.listdir(config_folder):
        config = utils.read_json(config_folder + config_file)

        height = config['height']
        width = config['width']
        depth = config['depth']

        i = 0
        base_directory = config['splitted_directory'] + '_' + str(config['multiplier']) + '_'\
                    + str(config['image_height']) + 'x' + str(config['image_width']) + '_'\
                    + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth']) + str(config['suffix'])\
                    + '/' + config_file[:-5]
        off_dir = config['off_dir'] + '_' + str(config['multiplier']) + '_' \
                  + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                  + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth']) + str(config['suffix'])

        if not os.path.exists(off_dir):
            os.makedirs(off_dir)

        for filename in os.listdir(base_directory):

            in_file = base_directory + '/' + filename
            m = mesh.Mesh.from_off(in_file)

            for j in range(config['multiplier']):
                # each version is done twice, once flipped
                angles = [0.]
                if config['flip']:
                    print('[Data] not flipping')
                    angles = [0., 180.]

                for angle in angles:
                    m_prime = m.copy()

                    rotation = [
                        0,
                        angle/180. * math.pi,
                        0
                    ]
                    m_prime.rotate(rotation)

                    rotation = [
                        (np.random.random()*2*config['max_x_rotation'] - config['max_x_rotation'])/180.0*math.pi,
                        (np.random.random()*2*config['max_y_rotation'] - config['max_y_rotation'])/180.0*math.pi,
                        (np.random.random()*2*config['max_z_rotation'] - config['max_z_rotation'])/180.0*math.pi
                    ]
                    m_prime.rotate(rotation)

                    scale = 1 + np.random.uniform(config['min_scale'], config['max_scale'])
                    scales = [scale, scale, scale]
                    m_prime.scale(scales)

                    translation = [
                        np.random.uniform(config['min_x_translation'], config['max_x_translation']),
                        np.random.uniform(config['min_y_translation'], config['max_y_translation']),
                        np.random.uniform(config['min_z_translation'], config['max_z_translation'])
                    ]
                    m_prime.translate(translation)

                    out_file = off_dir + '/' + str(i) + '.off'
                    m_prime.to_off(out_file)

                    print('[Data] generated %s' % out_file)
                    i += 1