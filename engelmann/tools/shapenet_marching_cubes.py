import numpy as np
import mcubes
import sys
import os

def read_matrix(filepath):
    """
    Read matrix, i.e. from a text file containing one vector.

    :param filepath: file to read
    :type filepath: str
    :return: sdf values
    :rtype: numpy.ndarray
    """

    values = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        lines = [line for line in lines if line != '']

        for line in lines:
            parts = line.split(',')
            parts = [part.strip() for part in parts]
            parts = [float(part) for part in parts]

            values += parts

    return np.array(values)

def read_point_cloud(file):
    """
    Reads point cloud.

    :param file: path to file to read
    :type file: str
    :return: points
    :rtype: [numpy.ndarray]
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        lines = fp.readlines()
        lines = [line.strip() for line in lines if line.strip()]

        num_points = int(lines[0])
        lines = lines[1:]

        points = np.zeros((num_points, 3))

        i = 0
        for line in lines:
            parts = line.split(' ')
            assert len(parts) == 3, 'expected 3 points but got %d (%s)' % (len(parts), file)

            points[i, 0] = float(parts[0])
            points[i, 1] = float(parts[1])
            points[i, 2] = float(parts[2])
            i = i + 1

        return points

class Mesh:
    """
    Represents a mesh.
    """

    def __init__(self, vertices = [[]], faces = [[]]):
        """
        Construct a mesh from vertices and faces.

        :param vertices: list of vertices, or numpy array
        :type vertices: [[float]] or numpy.ndarray
        :param faces: list of faces or numpy array, i.e. the indices of the corresponding vertices per triangular face
        :type faces: [[int]] fo rnumpy.ndarray
        """

        self.vertices = np.array(vertices, dtype = float)
        """ (numpy.ndarray) Vertices. """

        self.faces = np.array(faces, dtype = int)
        """ (numpy.ndarray) Faces. """

        assert self.vertices.shape[1] == 3
        assert self.faces.shape[1] == 3

    def extents(self):
        """
        Get the extents.

        :return: (min_x, min_y, min_z), (max_x, max_y, max_z)
        :rtype: (float, float, float), (float, float, float)
        """

        min = [0]*3
        max = [0]*3

        for i in range(3):
            min[i] = np.min(self.vertices[:, i])
            max[i] = np.max(self.vertices[:, i])

        return tuple(min), tuple(max)

    def switch_axes(self, axis_1, axis_2):
        """
        Switch the two axes, this is usually useful for switching y and z axes.

        :param axis_1: index of first axis
        :type axis_1: int
        :param axis_2: index of second axis
        :type axis_2: int
        """

        temp = np.copy(self.vertices[:, axis_1])
        self.vertices[:, axis_1] = self.vertices[:, axis_2]
        self.vertices[:, axis_2] = temp

    def mirror(self, axis):
        """
        Mirror given axis.

        :param axis: axis to mirror
        :type axis: int
        """

        self.vertices[:, axis] *= -1

    def scale(self, scales):
        """
        Scale the mesh in all dimensions.

        :param scales: tuple of length 3 with scale for (x, y, z)
        :type scales: (float, float, float)
        """

        assert len(scales) == 3

        for i in range(3):
            self.vertices[:, i] *= scales[i]

    def translate(self, translation):
        """
        Translate the mesh.

        :param translation: translation as (x, y, z)
        :type translation: (float, float, float)
        """

        assert len(translation) == 3

        for i in range(3):
            self.vertices[:, i] += translation[i]

    def _rotate(self, R):

        self.vertices = np.dot(R, self.vertices.T)
        self.vertices = self.vertices.T

    def rotate(self, rotation):
        """
        Rotate the mesh.

        :param rotation: rotation in (angle_x, angle_y, angle_z); angles in radians
        :type rotation: (float, float, float
        :return:
        """
        assert len(rotation) == 3

        x = rotation[0]
        y = rotation[1]
        z = rotation[2]

        # rotation around the x axis
        R = np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
        self._rotate(R)

        # rotation around the y axis
        R = np.array([[math.cos(y), 0, math.sin(y)], [0, 1, 0], [-math.sin(y), 0, math.cos(y)]])
        self._rotate(R)

        # rotation around the z axis
        R = np.array([[math.cos(z), -math.sin(z), 0], [math.sin(z), math.cos(z), 0], [0, 0, 1]])
        self._rotate(R)

    def inv_rotate(self, rotation):
        """
        Rotate the mesh.

        :param rotation: rotation in (angle_x, angle_y, angle_z); angles in radians
        :type rotation: (float, float, float
        :return:
        """
        assert len(rotation) == 3

        x = rotation[0]
        y = rotation[1]
        z = rotation[2]

        # rotation around the x axis
        R = np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
        R = R.T
        self._rotate(R)

        # rotation around the y axis
        R = np.array([[math.cos(y), 0, math.sin(y)], [0, 1, 0], [-math.sin(y), 0, math.cos(y)]])
        R = R.T
        self._rotate(R)

        # rotation around the z axis
        R = np.array([[math.cos(z), -math.sin(z), 0], [math.sin(z), math.cos(z), 0], [0, 0, 1]])
        R = R.T
        self._rotate(R)

    def transform(self, transformation):
        """
        Transform a mesh using a 4x4 transformation matrix.

        :param transformation: 4x4 transformation matrix
        :type transformation: numpy.ndarray
        """

        vertices = self.vertices.copy()
        vertices = np.concatenate((vertices, np.ones((vertices.shape[0], 1))), axis=1)

        vertices = np.dot(transformation, vertices.T).T
        vertices[:, 0:3] /= np.repeat(vertices[:, 3].reshape((vertices.shape[0], 1)), 3, axis=1)

        self.vertices = vertices[:, 0:3]

    def copy(self):
        """
        Copy the mesh.

        :return: copy of the mesh
        :rtype: Mesh
        """

        mesh = Mesh(self.vertices.copy(), self.faces.copy())

        return mesh

    def to_off(self, filepath):
        """
        Write mesh to OFF.

        :param filepath: path to write file to
        :type filepath: str
        """

        faces = np.ones((self.faces.shape[0], 4), dtype = int)*3
        faces[:, 1:4] = self.faces[:, :].tolist()
        vertices = self.vertices.tolist()

        num_vertices = len(vertices)
        num_faces = len(faces)

        assert num_vertices > 0
        assert num_faces > 0

        with open(filepath, 'w') as fp:
            fp.write('OFF\n')
            fp.write(str(num_vertices) + ' ' + str(num_faces) + ' 0\n')

            for vertex in vertices:
                assert len(vertex) == 3, 'invalid vertex with %d dimensions found (%s)' % (len(vertex), file)
                fp.write(str(vertex[0]) + ' ' + str(vertex[1]) + ' ' + str(vertex[2]) + '\n')

            for face in faces:
                assert face[0] == 3, 'only triangular faces supported (%s)' % file
                assert len(face) == 4, 'faces need to have 3 vertices, but found %d (%s)' % (len(face), file)

                for i in range(len(face)):
                    assert face[i] >= 0 and face[i] < num_vertices, 'invalid vertex index %d (of %d vertices) (%s)' % (
                    face[i], num_vertices, file)

                    fp.write(str(face[i]))
                    if i < len(face) - 1:
                        fp.write(' ')

                fp.write('\n')

            # add empty line to be sure
            fp.write('\n')

    @staticmethod
    def from_volume(volume):
        """
        Create a mesh from a voxel grid/volume.

        :param volume: volume
        :type volume: numpy.ndarray
        :return: mesh
        :rtype: Mesh
        """

        height = volume.shape[0]
        width = volume.shape[1]
        depth = volume.shape[2]

        vertices = []
        faces = []

        xx, yy, zz = np.where(volume > 0.5)

        for i in range(len(xx)):
            v000 = (yy[i], xx[i], zz[i])                # 0
            v001 = (yy[i], xx[i], zz[i] + 1)            # 1
            v010 = (yy[i], xx[i] + 1, zz[i])            # 2
            v011 = (yy[i], xx[i] + 1, zz[i] + 1)        # 3
            v100 = (yy[i] + 1, xx[i], zz[i])            # 4
            v101 = (yy[i] + 1, xx[i], zz[i] + 1)        # 5
            v110 = (yy[i] + 1, xx[i] + 1, zz[i])        # 6
            v111 = (yy[i] + 1, xx[i] + 1, zz[i] + 1)    # 7

            n = len(vertices)
            f1 = [n + 0, n + 2, n + 4]
            f2 = [n + 4, n + 2, n + 6]
            f3 = [n + 1, n + 3, n + 5]
            f4 = [n + 5, n + 3, n + 7]
            f5 = [n + 0, n + 1, n + 2]
            f6 = [n + 1, n + 3, n + 2]
            f7 = [n + 4, n + 5, n + 7]
            f8 = [n + 4, n + 7, n + 6]
            f9 = [n + 4, n + 0, n + 1]
            f10 = [n + 4, n + 5, n + 1]
            f11 = [n + 2, n + 3, n + 6]
            f12 = [n + 3, n + 7, n + 6]

            vertices.append(v000)
            vertices.append(v001)
            vertices.append(v010)
            vertices.append(v011)
            vertices.append(v100)
            vertices.append(v101)
            vertices.append(v110)
            vertices.append(v111)

            faces.append(f1)
            faces.append(f2)
            faces.append(f3)
            faces.append(f4)
            faces.append(f5)
            faces.append(f6)
            faces.append(f7)
            faces.append(f8)
            faces.append(f9)
            faces.append(f10)
            faces.append(f11)
            faces.append(f12)

        return Mesh(vertices, faces)

if __name__ == '__main__':
    # int dx = 60;
    # int dy = 40;
    # int dz = 60;
    # int x = ((i / dy / dz) % dx) - 30;
    # int y = (i / dz) % dy - 30;
    # int z = (i % dz) - 30;

    sdf_dir = sys.argv[1]
    assert os.path.exists(sdf_dir)

    off_dir = sys.argv[2]
    if not os.path.exists(off_dir):
        os.makedirs(off_dir)

    for i in range(1950):
        sdf_file = sdf_dir + '/' + str(i) + '.txt'
        sdf = read_matrix(sdf_file)
        sdf = sdf.reshape((60, 40, 60))

        vertices, triangles = mcubes.marching_cubes2(sdf, 0)
        mesh = Mesh(vertices, triangles)

        pose_file = sdf_dir + '/' + str(i) + '_pose.txt'
        pose = read_matrix(pose_file)
        pose = pose.reshape((4, 4))

        mesh.scale((0.1, 0.1, 0.1))
        mesh.translate((-3, -3, -3))

        mesh.transform(pose)

        # Eigen::Matrix4d trafo = Eigen::Matrix4d::Identity();
        # trafo(0,3) -= params.width/2;
        # trafo(1,3) -= params.height/2;
        # trafo(2,3) -= params.depth/2;
        # pointcloud->transform(trafo);
        # trafo = Eigen::Matrix4d::Identity()*6/params.width;
        # trafo(1,3) = 1;
        # trafo = Eigen::Matrix4d::Zero();
        # trafo(0,2) = 1;
        # trafo(1,1) = -1; // y points downward in KITTI camera coordinates
        # trafo(2,0) = -1;
        # pointcloud->transform(trafo);

        mesh.scale((1, -1, -1))
        mesh.switch_axes(0, 2)
        mesh.translate((0, -1, 0))
        mesh.scale((54/6, 54/6, 54/6))
        mesh.translate((54//2, 24//2, 24//2))

        off_file = off_dir + '/' + str(i) + '.off'
        mesh.to_off(off_file)
        print('Wrote ' + off_file)

