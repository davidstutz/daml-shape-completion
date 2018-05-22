"""
Post-process.
"""

import numpy as np
import os
import sys

sys.path.insert(1, os.path.realpath(__file__ + '../lib/'))
import utils
import ntpath
import common

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('[Data] Usage python 1_post_process.py config_folder')
        exit(1)

    config_file = sys.argv[1]
    assert os.path.exists(config_file), 'file %s does not exist' % config_file

    print('[Data] reading ' + config_file)
    config = utils.read_json(config_file)
    set = ntpath.basename(config_file)[:-5]

    truncation = config['truncation']
    input_sdfs = utils.read_hdf5(common.filename(config, 'input_sdf_file', '_f.h5', set))

    input_tsdfs = input_sdfs.copy()
    input_tsdfs[input_tsdfs > truncation] = truncation
    input_tsdfs[input_tsdfs < -truncation] = -truncation

    input_ltsdfs = input_tsdfs.copy()
    input_ltsdfs[input_ltsdfs > 0] = np.log(input_ltsdfs[input_ltsdfs > 0] + 1)
    input_ltsdfs[input_ltsdfs < 0] = - np.log(np.abs(input_ltsdfs[input_ltsdfs < 0]) + 1)

    input_tsdf_file = common.filename(config, 'input_tsdf_file', '_f.h5', set)
    input_ltsdf_file = common.filename(config, 'input_ltsdf_file', '_f.h5', set)

    utils.write_hdf5(input_tsdf_file, input_tsdfs)
    print('[Data] wrote ' + input_tsdf_file)
    utils.write_hdf5(input_ltsdf_file, input_ltsdfs)
    print('[Data] wrote ' + input_ltsdf_file)
