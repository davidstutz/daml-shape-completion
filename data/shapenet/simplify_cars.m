addpath('libsimplifymesh');
addpath('libsimplifymesh/qslim');
addpath('libsimplifymesh/remesher');

OFF_IN_DIR = '/work/data/shapenet/off_cars/';
OFF_OUT_DIR = '/work/data/shapenet/off_simplified/';
mkdir(OFF_OUT_DIR);

files = dir([OFF_IN_DIR '*.off']);
for i = 1: size(files)
    off_file = [OFF_IN_DIR files(i).name];
    [vertices, faces] = read_off(off_file);

    model = {};
    model.mesh = {};
    model.mesh.faces = faces;
    model.mesh.vertices = vertices;
    model = semiConvexHull(model, 0);
    
    off_file = [OFF_OUT_DIR sprintf('%04d', i) '.off'];
    write_off(off_file, model.hull.faces, model.hull.vertices);
    off_file
end

function [vertices, faces] = read_off(filename)
% Read a triangular mesh from OFF.
%
% INPUT:
% filename:     filename to read from
%
% OUTPUT:
% faces:        N x 3 matrix of vertex indices (1-based)
% vertices:     K x 3 matrix of vertices (1-based)

    fid = fopen(filename,'r');
    if fid == - 1
        error('Can''t open the file.');
        return;
    end

    str = fgets(fid);   % -1 if eof
    if ~strcmp(str(1:3), 'OFF')
        error('The file is not a valid OFF one.');    
    end

    str = fgets(fid);
    [a, str] = strtok(str);
    nvert = str2num(a);
    [a, str] = strtok(str);
    nface = str2num(a);

    [A, cnt] = fscanf(fid,'%f %f %f', 3*nvert);
    if cnt ~= 3*nvert
        warning('Problem in reading vertices.');
    end
    A = reshape(A, 3, cnt/3);
    vertices = A';

    [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
    if cnt~=4*nface
        warning('Problem in reading faces.');
    end
    A = reshape(A, 4, cnt/4);
    faces = A(2:4,:)'+1;

    fclose(fid);
end

function write_off(filename, faces, vertices)
% Write a triangular mesh to OFF.
%
% INPUT:
% filename:     filename to write to
% faces:        N x 3 matrix of vertex indices (1-based)
% vertices:     K x 3 matrix of vertices (1-based)
%
% David Stutz <david.stutz@rwth-aachen.de>

    file = fopen(filename, 'w');

    if size(faces, 2) ~= 3
        error('faces should be matrix of size N x 3');
    end

    if size(vertices, 2) ~= 3,
        error('vertices should be matrix of size K x 3');
    end

    n_vertices = size(vertices, 1);
    n_faces = size(faces, 1);
    fprintf (file, 'OFF\n%d %d 0\n', n_vertices, n_faces);

    for v = 1: n_vertices
        fprintf(file, '%f %f %f\n', vertices(v, 1), vertices(v, 2), vertices(v, 3));
    end

    for f = 1: n_faces
        fprintf(file, '%d', 3);
        for v = 1: 3
            % MatLab is 1-based, all others are 0-based!
            fprintf(file, ' %d', faces(f, v) - 1);
        end
        fprintf(file, '\n');
    end

    fclose(file);
end