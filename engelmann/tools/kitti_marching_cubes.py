import numpy as np
import mcubes
import sys
import os

def read_sdf(filepath):
    """
    Read SDF, i.e. from a text file containing one 60*40*60 dim vector.

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

class BoundingBox:
    """
    Represents a bounding box by center, size and orientation. Note that the input will be
    directly read from the KITTI XML files; thus, the z-axis corresponds to height instead
    of depth and the y-axis corresponds to depth instead of height.

    This class will take care of swapping the axes correspondingly.
    """

    # https://gist.github.com/phn/1111712/35e8883de01916f64f7f97da9434622000ac0390
    @staticmethod
    def _normalize_angle(num, lower=0.0, upper=360.0, b=False):
        """Normalize number to range [lower, upper) or [lower, upper].
        Parameters
        ----------
        num : float
            The number to be normalized.
        lower : float
            Lower limit of range. Default is 0.0.
        upper : float
            Upper limit of range. Default is 360.0.
        b : bool
            Type of normalization. See notes.
        Returns
        -------
        n : float
            A number in the range [lower, upper) or [lower, upper].
        Raises
        ------
        ValueError
          If lower >= upper.
        Notes
        -----
        If the keyword `b == False`, the default, then the normalization
        is done in the following way. Consider the numbers to be arranged
        in a circle, with the lower and upper marks sitting on top of each
        other. Moving past one limit, takes the number into the beginning
        of the other end. For example, if range is [0 - 360), then 361
        becomes 1. Negative numbers move from higher to lower
        numbers. So, -1 normalized to [0 - 360) becomes 359.
        If the keyword `b == True` then the given number is considered to
        "bounce" between the two limits. So, -91 normalized to [-90, 90],
        becomes -89, instead of 89. In this case the range is [lower,
        upper]. This code is based on the function `fmt_delta` of `TPM`.
        Range must be symmetric about 0 or lower == 0.
        Examples
        --------
        >>> normalize(-270,-180,180)
        90
        >>> import math
        >>> math.degrees(normalize(-2*math.pi,-math.pi,math.pi))
        0.0
        >>> normalize(181,-180,180)
        -179
        >>> normalize(-180,0,360)
        180
        >>> normalize(36,0,24)
        12
        >>> normalize(368.5,-180,180)
        8.5
        >>> normalize(-100, -90, 90, b=True)
        -80.0
        >>> normalize(100, -90, 90, b=True)
        80.0
        >>> normalize(181, -90, 90, b=True)
        -1.0
        >>> normalize(270, -90, 90, b=True)
        -90.0
        """
        from math import floor, ceil
        # abs(num + upper) and abs(num - lower) are needed, instead of
        # abs(num), since the lower and upper limits need not be 0. We need
        # to add half size of the range, so that the final result is lower +
        # <value> or upper - <value>, respectively.
        res = num
        if not b:
            if lower >= upper:
                raise ValueError("Invalid lower and upper limits: (%s, %s)" % (lower, upper))

            res = num
            if num > upper or num == lower:
                num = lower + abs(num + upper) % (abs(lower) + abs(upper))
            if num < lower or num == upper:
                num = upper - abs(num - lower) % (abs(lower) + abs(upper))

            res = lower if res == upper else num
        else:
            total_length = abs(lower) + abs(upper)
            if num < -total_length:
                num += ceil(num / (-2 * total_length)) * 2 * total_length
            if num > total_length:
                num -= floor(num / (2 * total_length)) * 2 * total_length
            if num > upper:
                num = total_length - num
            if num < lower:
                num = -total_length - num

            res = num * 1.0  # Make all numbers float, to be consistent

        return res

    def __init__(self, size, translation, rotation, source = 'kitti_raw', type = '', meta = ''):
        """
        Constructor.

        :param size: list of length 3 corresponding to (width, height, length)
        :type size: [float]
        :param translation: list of length 3 corresponding to translation (x, y, z
        :type translation: [float]
        :param rotation: list of length 3 corresponding to rotation (x, y, z)
        :type rotation: [float]
        :param kitti: whether to convert from kitti axes
        :type kitti: bool
        :param meta: meta information
        :type meta:
        """

        self.size = []
        """ ([float]) Size of bounding box. """

        if source == 'kitti_raw':
            # original h w l
            # change to l w h
            self.size = [size[2], size[1], size[0]]
        elif source == 'kitti':
            # original h w l
            # change to l w h
            self.size = [size[2], size[1], size[0]]
        else:
            self.size = [size[0], size[1], size[2]]

        # bounding boxes always start at zero height, i.e. its center is
        # determined from height
        self.translation = []
        """ ([float]) Translation, i.e. center, of bounding box. """

        R0_rect = np.array([9.999239000000e-01, 9.837760000000e-03, -7.445048000000e-03,
                            -9.869795000000e-03, 9.999421000000e-01, -4.278459000000e-03,
                            7.402527000000e-03, 4.351614000000e-03, 9.999631000000e-01]).reshape(3, 3)
        R0_rect_inv = np.linalg.inv(R0_rect)
        R_velo_to_cam = np.array([7.533745000000e-03, -9.999714000000e-01, -6.166020000000e-04,
                                  1.480249000000e-02, 7.280733000000e-04, -9.998902000000e-01,
                                  9.998621000000e-01, 7.523790000000e-03, 1.480755000000e-02]).reshape(3, 3)
        R_velo_to_cam_inv = np.linalg.inv(R_velo_to_cam)
        Tr_velo_to_cam = np.array([-4.069766000000e-03,
                                   -7.631618000000e-02,
                                   -2.717806000000e-01]).reshape(3)

        if source == 'kitti_raw':
            # bounding boxes are in camera coordinate system
            # x = right, y = down, z = forward, see paper
            # point clods are in velodyne coordinate system
            # x = forward, y = left, z = up

            # incoming translation is in velodyne coordinates

            self.translation = [-translation[1], translation[2] + self.size[1] / 2, translation[0]]
        elif source == 'kitti':
            # bounding boxes are in camera coordinate system
            # x = right, y = down, z = forward, see paper
            # point clods are in velodyne coordinate system
            # x = forward, y = left, z = up

            # incoming translation is in camera coordinates

            self._translation = np.dot(R_velo_to_cam_inv, np.dot(R0_rect_inv, np.array(translation)) - Tr_velo_to_cam)

            # now in velodyne coordinates

            self.translation = [-self._translation[1], self._translation[2] + self.size[1] / 2, self._translation[0]]
        else:
            self.translation = [translation[0], translation[1], translation[2]]

        self.rotation = []
        """ ([float]) Rotation of bounding box. """

        if source == 'kitti_raw':
            # incoming rotation is in velodyne coordinates
            self.rotation = [rotation[1], self._normalize_angle(rotation[2] + math.pi/2, -math.pi, math.pi), rotation[0]]
        elif source == 'kitti':
            # incoming rotation is in camera coordinates
            self.rotation = [rotation[0], -rotation[2], rotation[1]]
        else:
            self.rotation = [rotation[0], rotation[1], rotation[2]]

        self.type = type
        """ (string) Type. """

        self.meta = meta
        """ (string) Meta information. """

        #print('[Data] (%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f)' \
        #      % (self.size[0], self.size[1], self.size[2],
        #         self.translation[0], self.translation[1], self.translation[2],
        #         self.rotation[0], self.rotation[1], self.rotation[2]))

    def copy(self):
        """
        Copy the bounding box (deep copy).

        :return: copied bounding box
        :rtype:BoundingBox
        """

        return BoundingBox(self.size, self.translation, self.rotation, False, self.meta)

    @staticmethod
    def from_kitti(filepath):
        """
        Read bounding boxes from KITTI.

        :param filepath: bounding box file
        :type filepath: str
        :return: list of bounding boxes
        :rtype: [BoundingBox]
        """

        bounding_boxes = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines if line.strip() != '']

            for line in lines:
                parts = line.split(' ')
                assert len(parts) == 15 or len(parts) == 16, "invalid number of parts in line %s" % filepath

                bounding_boxes.append(BoundingBox([float(parts[8]), float(parts[9]), float(parts[10])],
                                 [float(parts[11]), float(parts[12]), float(parts[13])],
                                 [0., 0., float(parts[14])], 'kitti', parts[0], filepath))

        return bounding_boxes

    def __str__(self):
        """
        Convert to string representation.

        :return: bounding box as string
        :rtype: str
        """

        return '%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %s' % (
            round(self.size[0], 2), round(self.size[1], 2), round(self.size[2], 2),
            round(self.translation[0], 2), round(self.translation[1], 2), round(self.translation[2], 2),
            round(self.rotation[0], 2), round(self.rotation[1], 2), round(self.rotation[2], 2),
            self.meta
        )

def read_bounding_boxes(file):
    """
    Reads bounding boxes.

    :param file: path to file to read
    :type file: str
    :return: bounding boxes
    :rtype: [BoundingBox]
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        lines = fp.readlines()
        lines = [line.strip() for line in lines if line.strip()]
        lines = lines[1:]

        bounding_boxes = []
        for line in lines:
            parts = line.split(' ')
            parts = [part.strip() for part in parts if part]

            assert len(parts) == 9 or len(parts) == 10, 'invalid bounding box line: %s' % line

            size = (float(parts[0]), float(parts[1]), float(parts[2]))
            translation = (float(parts[3]), float(parts[4]), float(parts[5]))
            rotation = (float(parts[6]), float(parts[7]), float(parts[8]))

            meta = ''
            if len(parts) > 9:
                meta = parts[9]

            bounding_boxes.append(BoundingBox(size, translation, rotation, '', '', meta))

        return bounding_boxes

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

    bounding_box_file = sys.argv[3]
    assert os.path.exists(bounding_box_file)

    bounding_boxes = read_bounding_boxes(bounding_box_file)

    off_dir = sys.argv[3]
    if not os.path.exists(off_dir):
        os.makedirs(off_dir)

    for i in range(7118):
        sdf_file = sdf_dir + '/' + str(i) + '.txt'
        sdf = read_sdf(sdf_file)
        sdf = sdf.reshape((60, 40, 60))

        vertices, triangles = mcubes.marching_cubes2(sdf, 0)
        mesh = Mesh(vertices, triangles)

        mesh.scale((0.1, 0.1, 0.1))
        mesh.translate((-3, -3, -3))
        #mesh.scale((1./60, 1./60, 1./60))
        #mesh.scale((1, -1, -1))
        #mesh.switch_axes(0, 2)
        #mesh.scale((54, 54, 54))
        #mesh.translate((54//2, 24//2, 24//2))

        #txt_file = txt_dir + '/' + str(i) + '.txt'
        #points = read_point_cloud(txt_file)
        #points[:, 0] -= 54 // 2
        #points[:, 1] -= 24 // 2
        #points[:, 2] -= 24 // 2
        #points *= bounding_boxes[i].size[0]/54
        #points[:, 1] -= bounding_boxes[i].translation[1]

        #transformation = np.zeros((3, 3))
        #transformation[0, 2] = 1
        #transformation[1, 1] = -1
        #transformation[2, 0] = -1
        #points = np.dot(transformation, points.T).T

        #points[:, 1] -= 7.631618000000e-02
        #mesh.vertices = np.concatenate((mesh.vertices, points), axis=0)

        mesh.translate((0, 7.631618000000e-02, 0))
        mesh.scale((1, -1, -1))
        mesh.switch_axes(0, 2)
        mesh.translate((0, bounding_boxes[i].translation[1], 0))
        mesh.scale((54/bounding_boxes[i].size[0], 54/bounding_boxes[i].size[0], 54/bounding_boxes[i].size[0]))
        mesh.translate((54//2, 24//2, 24//2))

        off_file = off_dir + '/' + str(i) + '.off'
        mesh.to_off(off_file)
        print('Wrote ' + off_file)

