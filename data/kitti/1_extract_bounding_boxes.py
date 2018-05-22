import TrackletXML

import math
import random
import sys
import os
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
from bounding_box import BoundingBox
import ntpath
import shutil
import common

def write_frames(filepath, bounding_boxes):

    with open(filepath, 'w') as f:
        f.write(str(len(bounding_boxes)) + '\n')
        for bounding_box in bounding_boxes:
            f.write(str(bounding_box) + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python extract_bounding_boxes.py config_file')
        exit(1)

    config_file = sys.argv[1]
    assert os.path.exists(config_file)

    config = utils.read_json(config_file)

    # read mapping
    # go through mapping
        # for each item, get corresponding XML file
        # TODO also copy the corresponding +-10 velodyne point clouds
        # go through XML file and extract bounding boxes for specific frame and +-10 frames
        # save bounding boxes in separate files per frame

    mapping = []
    with open(config['mapping_file'], 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines if line.strip() != '']

        for line in lines:
            parts = line.split(' ')
            assert len(parts) == 3

            mapping.append({'drive': parts[1], 'frame': parts[2]})
            print('[Data] read ' + parts[1] + ' ' + parts[2])

    rand = []
    with open(config['rand_file'], 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines if line.strip() != '']

        assert len(lines) == 1

        parts = lines[0].split(',')
        rand = [int(part) for part in parts]

    indices = []
    with open(config['split_file'], 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines if line.strip() != '']
        indices = [int(line) for line in lines]

    drives = []
    for j in indices:
        array = mapping[rand[j] - 1]
        print('[Data] processing ' + array['drive'] + ' ' + array['frame'] + ' -> ' + str(j))
        xml_file = config['xml_directory'] + '/' + array['drive'] + '/tracklet_labels.xml'

        if os.path.exists(xml_file):
            print('[Data] reading ' + xml_file)
            tracklets = TrackletXML.parseXML(xml_file)
            bbs_by_frame = {}

            for tracklet in tracklets:
                # sort for cars.
                if tracklet.objectType.lower() == config['object_type']:
                    frame_id = tracklet.firstFrame

                    for i in range(tracklet.nFrames):

                        frame = '%04d' % (frame_id + i)
                        if not frame in bbs_by_frame:
                            bbs_by_frame[frame] = []

                        bounding_box = BoundingBox(tracklet.size, tracklet.trans[i, :], tracklet.rots[i, :])
                        bounding_box.size[0] += 2 * config['padding'] * bounding_box.size[0]
                        bounding_box.size[1] += 2 * config['padding'] * bounding_box.size[1]
                        bounding_box.size[2] += 2 * config['padding'] * bounding_box.size[2]

                        bbs_by_frame[frame].append(bounding_box)

            set = ntpath.basename(config_file)[:-5]
            bb_directory = common.filename(config, 'bounding_box_gt_txt_directory', '', set)
            velodyne_directory = common.filename(config, 'velodyne_gt_bin_directory', '', set)

            if not os.path.exists(bb_directory):
                os.makedirs(bb_directory)
            if not os.path.exists(velodyne_directory):
                os.makedirs(velodyne_directory)

            frame_id = int(array['frame'])
            if ('%04d' % frame_id) in bbs_by_frame.keys():
                for k in range(-config['gt_range'], config['gt_range'] + 1, config['gt_skip']):

                    frame_id = int(array['frame']) + k
                    if ('%04d' % frame_id) in bbs_by_frame.keys():

                        from_bin = config['xml_directory'] + '/' + array['drive'] + '/velodyne_points/data/%010d.bin' % frame_id
                        print('[Data] ' + array['drive'] + ' ' + str(frame_id) + ': ' + from_bin)

                        if os.path.exists(from_bin):
                            to_bin = velodyne_directory + '/%06d_%d.bin' % (j, k)
                            shutil.copy2(from_bin, to_bin)

                            bb_file = bb_directory + '/%06d_%d.txt' % (j, k)
                            print('[Data] ' + array['drive'] + ' ' + str(frame_id) + ': ' + bb_file)
                            write_frames(bb_file, bbs_by_frame[('%04d' % frame_id)])

                if not array['drive'] in drives:
                    drives.append(array['drive'])
            else:
                print('[Data] ' + array['drive'] + ' ' + array['frame'] + ': frame not found')
        else:
            print('[Data] ' + array['drive'] + ' ' + array['frame'] + ': drive xml not found')

    print(drives)