import math
import random
import sys
import os
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
from bounding_box import BoundingBox
import ntpath
import common

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python split_kitti.py config_file')
        exit(1)

    config_file = sys.argv[1]
    assert os.path.exists(config_file)

    config = utils.read_json(config_file)

    set = ntpath.basename(config_file)[:-5]
    print('[Data] processing ' + set)

    split = []
    with open(config ['split_file'], 'r') as f:
        split = f.readlines()
        split = [int(line.strip()) for line in split if line.strip() != '']

    bounding_box_directory = common.filename(config, 'bounding_box_directory', '', set)
    label_directory = config['label_directory']
    multiplier = config['multiplier']

    if not os.path.exists(bounding_box_directory):
        os.makedirs(bounding_box_directory)

    num_bounding_boxes = 0
    for i in split:
        bounding_box_file = label_directory + '/%06d.txt' % i
        bounding_boxes = BoundingBox.from_kitti(bounding_box_file)

        multiplied_bounding_boxes = []
        for bounding_box in bounding_boxes:
            if bounding_box.type.lower() == config['object_type']:
                for j in range(multiplier):
                    current_bounding_box = bounding_box.copy()
                    current_bounding_box.size[0] += 2 * config['padding'] * current_bounding_box.size[0]
                    current_bounding_box.size[1] += 2 * config['padding'] * current_bounding_box.size[1]
                    current_bounding_box.size[2] += 2 * config['padding'] * current_bounding_box.size[2]

                    multiplied_bounding_boxes.append(current_bounding_box)

        bounding_box_file = bounding_box_directory + '/%06d.txt' % i
        with open(bounding_box_file, 'w') as f:
            f.write(str(len(multiplied_bounding_boxes)) + '\n')
            for bounding_box in multiplied_bounding_boxes:
                f.write(str(bounding_box) + '\n')

        print('[Data] wrote ' + bounding_box_file)

        num_bounding_boxes += len(multiplied_bounding_boxes)

    num_bounding_box_file = common.filename(config, 'num_bounding_box_file', 'txt', set)
    with open(num_bounding_box_file, 'w') as f:
        f.write(str(num_bounding_boxes))
