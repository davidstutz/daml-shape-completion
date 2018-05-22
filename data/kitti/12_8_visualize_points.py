import argparse
import sys
import os
sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
from blender_utils import *
import common
import json
import re
import ntpath

def read_json(file):
    """
    Read a JSON file.

    :param file: path to file to read
    :type file: str
    :return: parsed JSON as dict
    :rtype: dict
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        return json.load(fp)

def read_ordered_directory(dir, extension = None):
    """
    Gets a list of file names ordered by integers (if integers are found
    in the file names).

    :param dir: path to directory
    :type dir: str
    :param extension: extension to filter for
    :type extension: str
    :return: list of file names
    :rtype: [str]
    """

    # http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically
    def get_int(value):
        """
        Convert the input value to integer if possible.

        :param value: mixed input value
        :type value: mixed
        :return: value as integer, or value
        :rtype: mixed
        """

        try:
            return int(value)
        except:
            return value

    def alphanum_key(string):
        """
        Turn a string into a list of string and number chunks,
        e.g. "z23a" -> ["z", 23, "a"].

        :param string: input string
        :type string: str
        :return: list of elements
        :rtype: [int|str]
        """

        return [get_int(part) for part in re.split('([0-9]+)', string)]

    def sort_filenames(filenames):
        """
        Sort the given list by integers if integers are found in the element strings.

        :param filenames: file names to sort
        :type filenames: [str]
        """

        filenames.sort(key = alphanum_key)

    assert os.path.exists(dir), 'directory %s not found' % dir

    filenames = [dir + '/' + filename for filename in os.listdir(dir)]
    if extension is not None:
        filenames = [filename for filename in filenames if filename[-len(extension):] == extension]

    sort_filenames(filenames)

    return filenames

if __name__ == '__main__':

    try:
        argv = sys.argv[sys.argv.index("--") + 1:]
    except ValueError:
        log('[Error] "--" not found, call as follows:', LogLevel.ERROR)
        log('[Error] $BLENDER --background --python 12_1_visualize_gt.py -- 1>/dev/null config_folder', LogLevel.ERROR)
        exit()

    if len(argv) < 1:
        log('[Error] not enough parameters, call as follows:', LogLevel.ERROR)
        log('[Error] $BLENDER --background --python 12_1_visualize_gt.py -- 1>/dev/null config_folder', LogLevel.ERROR)
        exit()

    config_file =argv[0]
    assert os.path.exists(config_file), 'file %s does not exist' % config_file

    config = read_json(config_file)
    set = ntpath.basename(config_file)[:-5]

    height = config['height']
    width = config['width']
    depth = config['depth']
    scale = 1./max(height, depth, width)

    txt_directory = common.filename(config, 'bounding_box_txt_directory', '', set)
    txt_files = read_ordered_directory(txt_directory)

    vis_directory = common.filename(config, 'vis_dir', '', set)
    if not os.path.isdir(vis_directory):
        os.makedirs(vis_directory)

    N = 30
    log('[Data] %d samples' % len(txt_files))
    for i in range(N):
        n = i * (len(txt_files) // N)
        txt_file = txt_files[n]

        camera_target = initialize()
        txt_material = make_material('BRC_Material_Point_Cloud', (0.65, 0.23, 0.25), 1, True)

        load_txt(txt_file, 0.0075, txt_material, (-width*scale/2, -depth*scale/2, -height*scale/2), scale, 'xzy')

        rotation = (5, 0, -55)
        distance = 0.35
        png_file = vis_directory + '/%d_points.png' % n
        render(camera_target, png_file, rotation, distance)

        log('[Data] wrote %s' % png_file)