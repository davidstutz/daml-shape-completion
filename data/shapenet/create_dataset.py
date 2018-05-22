#!/usr/bin/env python
"""
Complete dataset creation.
"""

import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import  utils

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_folder = sys.argv[1] + '/'
    assert os.path.exists(config_folder)

    # split sets from original training set
    os.system('python split_cars.py ' + config_folder)

    # random perturbations
    os.system('python multiply_cars.py ' + config_folder)

    # go through the generated cuboids as OFF files and
    # voxelize the meshes using simple cube-triangle intersection
    for config_file in os.listdir(config_folder):
        os.system('./libvoxelizemesh/bin/voxelize_meshs ' + config_folder + config_file)

    # go through the cuboids as OFF files and render the
    # meshes to generate depth files
    os.system('python render_cars.py ' + config_folder)

    # take the depth files and compute the corresponding point clouds
    # and free space volumes, voxelize both
    for config_file in os.listdir(config_folder):
        os.system('./libvoxelizecloud/bin/voxelize_clouds ' + config_folder + config_file)

    # compute (truncated) distance transforms
    os.system('python fast_post_process_cars.py ' + config_folder)
    #os.system('python compute_sdf.py ' + config_folder)


