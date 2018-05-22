function mesh = remesher (mesh,mesh_verts,lloyd_iters)

time_ms = mod(round(now*1e9),1e8);
in_file = sprintf('/tmp/ri_%08d.off',time_ms);
out_file = sprintf('/tmp/ro_%08d.off',time_ms);

write_off(in_file,mesh.vertices,mesh.faces);
unix(sprintf('libsimplifymesh/remesher/remesher %s %s %d %d',in_file,out_file,mesh_verts,lloyd_iters));
[mesh.vertices mesh.faces] = read_off(out_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertex,face] = read_off(filename)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);

[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A';

[A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
if cnt~=4*nface
    warning('Problem in reading faces.');
end
A = reshape(A, 4, cnt/4);
face = A(2:4,:)'+1;

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_off(filename, vertex, face)

if nargin<4
    renormalize = 0;
end

if size(vertex,2)~=3
    vertex=vertex';
end
if size(vertex,2)~=3
    error('vertex does not have the correct format.');
end

if size(face,2)~=3
    face=face';
end
if size(face,2)~=3
    error('face does not have the correct format.');
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% header
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(vertex,1), size(face,1));

% write the points & faces
fprintf(fid, '%f %f %f\n', vertex');
fprintf(fid, '3 %d %d %d\n', face'-1);

fclose(fid);