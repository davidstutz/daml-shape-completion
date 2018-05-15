import os
import sys
import numpy as np

sys.path.insert(1, os.path.realpath(__file__ + '/../../lib/py/'))
import utils
import binvox_rw
import mcubes

def get_val_files(directory, suffix):
    """
    Get files in the test directory matching the given suffix (excluding extension) and return them
    in a dict where the key is the iteration number.

    :param directory: directory to look in
    :type directory: str
    :param suffix: suffix of files to filter
    :type suffix: str
    :return: files as dictionary where the key is the iteration
    :rtype: dict
    """

    files = dict()
    for filename in os.listdir(directory):
        underscore_pos = filename.find('_')
        if underscore_pos <= 0:
            continue

        point_pos = filename.find('.')
        if point_pos <= 0:
            continue

        if filename[underscore_pos + 1:point_pos] == suffix:
            iteration = int(filename[0:underscore_pos])
            files[iteration] = directory + '/' + filename

    return files

if __name__ == '__main__':
    config_file = 'clean.json'
    if len(sys.argv) > 1:
        config_file = sys.argv[1];

    assert os.path.exists(config_file)

    print('[Validation] reading ' + config_file)
    config = utils.read_json(config_file)

    base_directory = config['base_directory'] + '/'
    predictions_files = get_val_files(base_directory, 'predictions')
    predictions_file = predictions_files[sorted(predictions_files, key=int)[-1]]

    print('[Validation] reading ' + predictions_file)
    predictions = utils.read_hdf5(predictions_file)
    predictions = np.squeeze(predictions)
    predictions = predictions[:, 0]

    off_directory = base_directory + config['off_directory'] + '/'
    utils.makedir(off_directory)

    for n in range(predictions.shape[0]):
        vertices, triangles = mcubes.marching_cubes(-predictions[n][1].transpose(1, 0, 2), 0)
        sdf_file = off_directory + str(n) + '.off'
        mcubes.export_off(vertices, triangles, sdf_file)
        print('[Validation] wrote ' + sdf_file)

    binvox_directory = base_directory + config['binvox_directory'] + '/'
    utils.makedir(binvox_directory)

    for n in range(predictions.shape[0]):
        model = binvox_rw.Voxels(predictions[n][0] > 0.5, predictions[n][0].shape, (0, 0, 0), 1)
        binvox_file = binvox_directory + str(n) + '.binvox'
        with open(binvox_file, 'w') as fp:
            model.write(fp)
            print('[Validation] wrote ' + binvox_file)

    if 'validation_outputs' in config.keys() and 'validation_lsdf_outputs' in config.keys():
        data_directory = config['data_directory'] + '/'
        print('[Validation] reading ' + data_directory + config['validation_outputs'])
        targets_1 = utils.read_hdf5(data_directory + config['validation_outputs'])
        print('[Validation] reading ' + data_directory + config['validation_lsdf_outputs'])
        targets_2 = utils.read_hdf5(data_directory + config['validation_lsdf_outputs'])
        targets = np.concatenate((np.expand_dims(targets_1, axis=1), np.expand_dims(targets_2, axis=1)), axis=1)
        targets = np.squeeze(targets)

        occ_predictions = predictions[:, 0]
        sdf_predictions = predictions[:, 1]
        occ_targets = targets[:, 0]
        sdf_targets = targets[:, 1]

        occ_predictions_thresh = occ_predictions.copy()
        occ_predictions_thresh[occ_predictions_thresh > 0.5] = 1
        occ_predictions_thresh[occ_predictions_thresh < 0.5] = 0

        sdf_predictions_thresh = sdf_predictions.copy()
        sdf_predictions_thresh[sdf_predictions_thresh > 0] = 0
        sdf_predictions_thresh[sdf_predictions_thresh < 0] = 1

        occ_error = np.sum(np.abs(occ_predictions - occ_targets)) / (
        predictions.shape[0] * predictions.shape[2] * predictions.shape[3] * predictions.shape[4])
        occ_error_thresh = np.sum(np.abs(occ_predictions_thresh - occ_targets)) / (
        predictions.shape[0] * predictions.shape[2] * predictions.shape[3] * predictions.shape[4])
        sdf_error = np.sum(np.abs(sdf_predictions - sdf_targets)) / (
        predictions.shape[0] * predictions.shape[2] * predictions.shape[3] * predictions.shape[4])
        sdf_error_thresh = np.sum(np.abs(sdf_predictions_thresh - occ_targets)) / (
        predictions.shape[0] * predictions.shape[2] * predictions.shape[3] * predictions.shape[4])

        results_file = base_directory + config['results_file']
        with open(results_file, 'w') as f:
            f.write('[Validation] occ error: ' + str(occ_error) + '\n')
            f.write('[Validation] occ error+thresh: ' + str(occ_error_thresh) + '\n')
            f.write('[Validation] sdf error: ' + str(sdf_error) + '\n')
            f.write('[Validation] sdf error+thresh: ' + str(sdf_error_thresh) + '\n')

            print('[Validation] occ error: ' + str(occ_error))
            print('[Validation] occ error+thresh: ' + str(occ_error_thresh))
            print('[Validation] sdf error: ' + str(sdf_error))
            print('[Validation] sdf error+thresh: ' + str(sdf_error_thresh))

            print('[Validation] INSTRUCTIONS:')
            print('[Validation] You can use the tools in lib/blender to visualize the generated OFF and BINVOX files.')
            print('[Validation] The evaluation results can be found in ' + results_file + '.')
            print('[Validation] To evaluate the mesh-to-mesh distance, please consult the documentation.')
