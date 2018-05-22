import numpy as np
import os

class PointCloud:

    def __init__(self, points=None):
        if points is None:
            self.points = np.zeros((0, 3))
        else:
            self.points = np.array(points)

    @staticmethod
    def from_txt(filepath):
        assert os.path.exists(filepath)

        with open(filepath, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines if line.strip() != '']

            num_points = int(lines[0])
            points = np.zeros((num_points, 3))
            assert num_points == len(lines) - 1

            for i in range(0, num_points):
                line = lines[i + 1]

                parts = line.split(' ')
                assert len(parts) == 3, "invalid line: %s" % line

                for j in range(3):
                    points[i, j] = float(parts[j])

        return PointCloud(points)

    def to_txt(self, filepath):

        with open(filepath, 'w') as f:
            f.write(str(self.points.shape[0]) + '\n')

            for n in range(self.points.shape[0]):
                f.write(str(self.points[n, 0]) + ' '  + str(self.points[n, 1]) + ' ' + str(self.points[n, 2]) + '\n')

    def to_ply(self, filepath):

        with open(filepath, 'w') as f:
            f.write('ply\n')
            f.write('format ascii 1.0\n')
            #f.write('format binary_little_endian 1.0\n')
            #f.write('format binary_big_endian 1.0\n')
            f.write('element vertex ' + str(self.points.shape[0]) + '\n')
            f.write('property float x\n')
            f.write('property float y\n')
            f.write('property float z\n')
            f.write('property uchar red\n')
            f.write('property uchar green\n')
            f.write('property uchar blue\n')
            f.write('end_header\n')

            for n in range(self.points.shape[0]):
                f.write(str(self.points[n, 0]) + ' '  + str(self.points[n, 1]) + ' ' + str(self.points[n, 2]))
                f.write(' 0 0 0\n')