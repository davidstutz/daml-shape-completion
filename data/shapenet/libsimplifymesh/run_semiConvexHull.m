function run_semiConvexHull()
close all;

addpath('remesher');
addpath('qslim');

% file to process
model_file = 'data/0003.mat';
hull_file = 'data/output.mat';

% parameters
w = 5;
gamma = 0.03;
s = [1.2 1.05];

% load model
model = load(model_file);
model.mesh.faces = int32(model.mesh.faces);
%model.model = data.model.model;
%model.mesh.vertices = data.model.vertices;
%model.mesh.faces = data.model.faces;
%model.scale = data.model.scale;
%model.mesh.vertices = model.mesh.vertices-ones(size(model.mesh.vertices,1),1)*[0 0 0.5];
%model.mesh.vertices = model.mesh.vertices.*(ones(size(model.mesh.vertices,1),1)*(model.scale));

% sample surface points
model.mesh = sampleSurfacePoints(model.mesh,50000);

% compute convex hull
model.hull = model.mesh;  
dt = DelaunayTri(model.mesh.vertices);
model.hull.vertices = dt.X;
model.hull.faces = convexHull(dt); 

% transform convex hull into a uniform mesh
num_vertices = 5000;
model.hull = reducepatch(model.hull,100);
model.hull = remesher(model.hull,num_vertices,100);
model.hull = remesher(model.hull,num_vertices,100);

% compute average edge length
l = 0;
for i=1:size(model.hull.faces,1)
  l = l+norm(model.hull.vertices(model.hull.faces(i,1),:)-model.hull.vertices(model.hull.faces(i,2),:));
  l = l+norm(model.hull.vertices(model.hull.faces(i,2),:)-model.hull.vertices(model.hull.faces(i,3),:));
  l = l+norm(model.hull.vertices(model.hull.faces(i,3),:)-model.hull.vertices(model.hull.faces(i,1),:));
end
l = l/(size(model.hull.faces,1)*3);
l = 0.8*l;
l2 = l^2;

% optimize hull
figure; view(145,15); axis tight; hold on;
for i=1:length(s)
  model.hull.vertices = s(i)*model.hull.vertices;
  for j=1:15
    fprintf('Iteration: %d - %d (s: %.3f)\n',i,j,s(i));
    [E,dE] = modelEnergy(model,l2,w);
    model.hull.vertices = model.hull.vertices-gamma*dE;
    plotHull(model.hull);
    pause(0.1); refresh;
  end
  model.hull = remesher(model.hull,num_vertices,100);
end

% reduce hull
model.hull = qslim(model.hull,'-t',1000,'-m',1000,'-O',3);

% remove sampled points & area
model.mesh = rmfield(model.mesh,'points');
model.mesh = rmfield(model.mesh,'area');

% visualize hull
plotHull(model.hull);
pause(0.1); refresh;

% save hull
save(hull_file,'model');
disp('done!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHull (hull)

cla; hold on;
trisurf(hull.faces,hull.vertices(:,1),hull.vertices(:,2),hull.vertices(:,3), ...
        'FaceColor', 'red','FaceAlpha', 0.5,'EdgeAlpha', 0.7);
axis equal;
axis tight;

