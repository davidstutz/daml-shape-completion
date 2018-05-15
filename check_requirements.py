try:
    import os
    import sys
    import re
    import json
    import h5py
    import numpy
    import mcubes
    import zipfile

    print('[Requirements] Seems that all requirements are met. However, this does not guarantee that everything will work.')
    print('[Requirements] Please consult the documentation for details on the requirements and further steps.')
except ImportError as e:
    print('[Requirements] Seems like one or more requirements are not met.')
    print(e)