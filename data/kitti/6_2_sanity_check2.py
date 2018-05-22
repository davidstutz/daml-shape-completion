#!/usr/bin/env python
"""
Some sanity checks of the rendering.
"""

from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as mplt
import random
import numpy as np
import os
import sys
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
import visualization
from mesh import Mesh
from point_cloud import PointCloud
import ntpath
import common

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python create_cuboids.py config_folder')
        exit(1)

    config_file = sys.argv[1]
    assert os.path.exists(config_file)
    config = utils.read_json(config_file)
    dataset = ntpath.basename(config_file)[:-5]

    orig_dir = common.filename(config, 'bounding_box_txt_directory', '', dataset)
    gt_dir = common.filename(config, 'velodyne_gt_txt_directory', '', dataset)
    gt_files = utils.read_ordered_directory(gt_dir)

    for i in range(len(gt_files)):
        orig_file = orig_dir + '/' + ntpath.basename(gt_files[i])
        orig_point_cloud = PointCloud.from_txt(orig_file)
        gt_point_cloud = PointCloud.from_txt(gt_files[i])
        print('[Data] read ' + orig_file)
        print('[Data] read ' + gt_files[i])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        xx = gt_point_cloud.points[:, 0]
        yy = gt_point_cloud.points[:, 1]
        zz = gt_point_cloud.points[:, 2]

        ax.scatter(xx, yy, zz, c='b', marker='+', s=25)

        xx = orig_point_cloud.points[:, 0]
        yy = orig_point_cloud.points[:, 1]
        zz = orig_point_cloud.points[:, 2]
        ax.scatter(xx, yy, zz, c='r', marker='x', s=25)
        plt.show()