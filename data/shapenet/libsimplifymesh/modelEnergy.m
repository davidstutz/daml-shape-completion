function [E,dE] = modelEnergy (model,l2,w)

% shortcut for vertices and points
x = model.hull.vertices;
y = model.mesh.points;

% shortcuts for vertices and faces
V = model.hull.vertices;
F = model.hull.faces;
n = size(V,1); % number of vertices
m = size(F,1); % number of faces
k = [1 2; 2 3; 3 1];

% nearest neighbors
idx = knnsearch(y,x);

% init energy
E = w*sum(norm2(x-y(idx,:)));
dE = zeros(size(V));

% vertex energy
for i=1:n
  dE(i,:) = w*2*(x(i,:)-y(idx(i),:));
end

% edge energy
for i=1:m
  for j=1:3
    v1 = F(i,k(j,1));
    v2 = F(i,k(j,2));
    E = E + (norm2(V(v1,:)-V(v2,:))-l2)^2;
    dE(v1,:) = dE(v1,:) + 4*(norm2(V(v1,:)-V(v2,:))-l2)*(V(v1,:)-V(v2,:));
    dE(v2,:) = dE(v2,:) - 4*(norm2(V(v1,:)-V(v2,:))-l2)*(V(v1,:)-V(v2,:));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dE = numericalGradient (model,l2,w)

delta = 1e-5;

E = energy(model,l2,w);
dE = zeros(size(model.hull.vertices));

for i=1:size(model.hull.vertices,1)
  for j=1:size(model.hull.vertices,2)
    model2 = model;
    model2.hull.vertices(i,j) = model2.hull.vertices(i,j)+delta;
    dE(i,j) = energy(model2,l2,w)-E;
  end
end
dE = dE/delta; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = norm2 (x)

n = x(:,1).*x(:,1)+x(:,2).*x(:,2)+x(:,3).*x(:,3);
