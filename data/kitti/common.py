import os

def filename(config, key, ext = '.h5', set = ''):
    """
    Get the real file name by looking up the key in the config and suffixing.

    :param key: key to use in the config
    :type key: str
    :param ext: extension to use
    :type ext: str
    :param set: set name
    :type set: str
    :return: filepath
    :rtype: str
    """

    name = config[key] + '_'
    if set:
        name += set + '_'

    name += str(config['multiplier']) + '_' + str(config['height']) + 'x' + str(config['width']) + 'x' + str(config['depth'])\

    if ext:
        name += ext

    return name