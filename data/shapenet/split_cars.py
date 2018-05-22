#!/usr/bin/env python
"""
Pre-process modelnet.
"""

import random
import shutil

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    taken = []
    files = []
    numbers = {}

    for config_file in os.listdir(config_folder):
        config = utils.read_json(config_folder + config_file)
        model_dir = config['model_directory'] + '/'
        split_dir = config['splitted_directory'] + '_' + str(config['multiplier']) + '_' \
                    + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                    + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth']) + str(config['suffix']) + '/'

        set = config_file[:-5]
        set_dir = split_dir + '/' + set
        numbers[set] = []

        if not os.path.exists(split_dir):
            os.makedirs(split_dir)

        # get all files only once to guarantee the same order
        if len(files) == 0:
            files = os.listdir(model_dir)
            files = [file for file in files if file[-4:] == '.off']
            files = sorted(files)

            taken = [False] * len(files)

        if not os.path.exists(set_dir):
            os.makedirs(set_dir)

        print('[Data] choosing %s in %s' % (set, set_dir))
        n = int(config['fraction']*len(files))

        for i in range(n):
            j = 0
            k = 0
            while taken[j] and k < 1e6: # prevent infinite loop
                j = random.randint(0, len(files) - 1)
                k = k + 1

            taken[j] = True
            numbers[set].append(files[j])

            off_from_file = model_dir + '/' + files[j]
            off_to_file = set_dir + '/' + str(i) + '.off'
            shutil.copy(off_from_file, off_to_file)
            print('[Data] copied')
            print('   %s' % off_from_file)
            print('   %s' % off_to_file)

            # try to read the file to see if ModelNet is flawed
            utils.read_off(off_from_file)
            utils.read_off(off_to_file)

        off_file = config['off_file'] + '_' + str(config['multiplier']) + '_' \
                   + str(config['image_height']) + 'x' + str(config['image_width']) + '_' \
                   + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth']) + str(config['suffix']) + '.txt'
        with open(off_file, 'w') as f:
            for j in numbers[set]:
                f.write(str(j) + '\n')
            print('[Data] wrote' + off_file)