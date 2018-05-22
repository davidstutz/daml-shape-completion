"""
Visualize results of an experiment.
"""

import numpy as np
import os
import terminaltables
import scipy.ndimage.morphology as morph
from scipy import ndimage
from scipy import misc
import glob
import pickle
import math
import matplotlib
from matplotlib import pyplot as plt
import shutil

import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
from point_cloud import PointCloud
import ntpath
import common

# https://stackoverflow.com/questions/17958485/matplotlib-not-using-latex-font-while-text-usetex-true
matplotlib.rcParams['text.usetex'] = True

if __name__ == '__main__':

    config_file = sys.argv[1]
    dataset = ntpath.basename(config_file)[:-5]
    base_directory = dataset + '/'

    if not os.path.exists(base_directory):
        os.makedirs(base_directory)

    config = utils.read_json(config_file)

    #statistics = utils.read_hdf5('/work/data/shapenet_3dop/real_space_statistics_training_prior.h5')
    #statistics = 1 - statistics
    #print('[Data] read statistics')

    inputs = utils.read_hdf5(common.filename(config, 'input_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_file', '_f.h5', dataset))

    inputs_combined_gt = utils.read_hdf5(common.filename(config, 'input_combined_gt_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'input_combined_gt_file', '_f.h5', dataset))

    space = utils.read_hdf5(common.filename(config, 'part_space_file', '_f.h5', dataset))
    print('[Data] read ' + common.filename(config, 'part_space_file', '_f.h5', dataset))

    #statistics = statistics.reshape(1, 1, statistics.shape[0], statistics.shape[1], statistics.shape[2])
    #statistics = np.repeat(statistics, space.shape[0], axis=0)

    #print(space.shape, statistics.shape)
    #invalid_space = space*statistics

    points = []
    point_dir = common.filename(config, 'bounding_box_txt_directory', '', dataset) + '/'

    for file in os.listdir(point_dir):
        point_file = point_dir + file
        point_cloud = PointCloud.from_txt(point_file)
        #print('[Data] read ' + point_file)

        points.append(point_cloud.points.shape[0])

    points_combined_gt = []
    point_dir = common.filename(config, 'velodyne_gt_txt_directory', '', dataset) + '/'

    for file in os.listdir(point_dir):
        point_file = point_dir + file
        point_cloud = PointCloud.from_txt(point_file)
        #print('[Data] read ' + point_file)

        points_combined_gt.append(point_cloud.points.shape[0])

    frames = [1]*inputs.shape[0]
    frame_dir = common.filename(config, 'velodyne_individual_gt_txt_directory', '', dataset) + '/'

    for i in range(inputs.shape[0]):
        for k in range(-config['gt_range'], config['gt_range'] + 1, config['gt_skip']):
            if k == 0:
                continue;

            txt_file = frame_dir + '%d_%d_%d.txt' % (i, k, frames[i])
            print(txt_file)
            if os.path.exists(txt_file):
                frames[i] += 1

    frames = np.array(frames);
    points = np.array(points)
    points_combined_gt = np.array(points_combined_gt)
    observed_points = np.squeeze(np.sum(np.sum(np.sum(inputs, axis=4), axis=3), axis=2))
    observed_points_combined_gt = np.squeeze(np.sum(np.sum(np.sum(inputs_combined_gt, axis=4), axis=3), axis=2))
    observed_space = np.squeeze(np.sum(np.sum(np.sum(space, axis=4), axis=3), axis=2))
    #observed_invalid_space = np.squeeze(np.sum(np.sum(np.sum(invalid_space, axis=4), axis=3), axis=2))

    mask = np.zeros(inputs.shape)
    mask[inputs == 1] = 1
    mask[space == 1] = 1
    observed_total = float(np.sum(mask))

    statistics_file = common.filename(config, 'statistics_file', '.txt', dataset)
    with open(statistics_file, 'w') as f:
        f.write('[Data] N: ' + str(inputs.shape[0]) + '\n')
        f.write('[Data] pixels: ' + str(inputs.shape[2] * inputs.shape[3] * inputs.shape[4]) + '\n')
        f.write('[Data] points: ' + str(np.mean(points)) + '\n')
        f.write('[Data] observed frames: ' + str(np.sum(frames) / inputs.shape[0]) + '\n')
        f.write('[Data] observed points: ' + str(np.sum(points) / inputs.shape[0]) + '\n')
        f.write('[Data] observed points (combined GT): ' + str(np.sum(points_combined_gt) / inputs.shape[0]) + '\n')
        f.write('[Data] observed points voxels: ' + str(np.sum(observed_points) / inputs.shape[0]) + '\n')
        f.write('[Data] observed points voxels (combined GT): ' + str(np.sum(observed_points_combined_gt) / inputs.shape[0]) + '\n')
        f.write('[Data] observed space voxels: ' + str(np.sum(observed_space) / inputs.shape[0]) + '\n')
        f.write('[Data] observed total voxels: ' + str(observed_total / inputs.shape[0]) + '\n')
        #f.write('[Data] observed invalid space voxels: ' + str(observed_invalid_space / inputs.shape[0]) + '\n')

        print('[Data] N: ' + str(inputs.shape[0]))
        print('[Data] pixels: ' + str(inputs.shape[2]*inputs.shape[3]*inputs.shape[4]))
        print('[Data] points: ' + str(np.mean(points)))
        print('[Data] observed frames: ' + str(np.sum(frames)/inputs.shape[0]))
        print('[Data] observed points: ' + str(np.sum(points) / inputs.shape[0]))
        print('[Data] observed points (combined GT): ' + str(np.sum(points_combined_gt) / inputs.shape[0]))
        print('[Data] observed points voxels: ' + str(np.sum(observed_points)/inputs.shape[0]))
        print('[Data] observed points voxels (combined GT): ' + str(np.sum(observed_points_combined_gt) / inputs.shape[0]))
        print('[Data] observed space voxels: ' + str(np.sum(observed_space)/inputs.shape[0]))
        print('[Data] observed total voxels: ' + str(observed_total/inputs.shape[0]))
        #print('[Data] observed invalid space voxels: ' + str(np.sum(observed_invalid_space) / inputs.shape[0]))
