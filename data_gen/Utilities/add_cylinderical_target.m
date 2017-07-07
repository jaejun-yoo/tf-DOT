function [mesh, index] = add_cylinderical_target(mesh, position, radius, height, prop, dir, reg_label)

% Adds cylindrical (3D) target in (x or y or z direction) 
% to 3D standard mesh properties of target are given in 'prop'
% position = [x y z];
% radius = radius of cylindrical target
% height = [xmin xmax]or [ymin ymax] or [zmin zmax];
% prop = [mua mus ri];
% dir = 1 for x axis, 2 for y axis, 3 for z axis
% reg_label = region label (typically > 0)
% mesh = add_cylinderical_target(mesh, [10 0 50], 7.5, [-10 10], [0.01 1 1.3],1, 1);


% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = loadmesh(mesh);
end

mesh_2d = mesh;
mesh_2d.nodes(:,dir) = 0;
x = position(1);
y = position(2);
z = position(3);
temp = [x y z];
temp(:,dir) = 0;

r = radius;
dist = distance(mesh_2d.nodes(:,1:3),ones(length(mesh_2d.bndvtx),1),temp);
clear mesh_2d;

mua = prop(1);
mus = prop(2);
kappa = 1./(3.*(mua+mus));
ri = prop(3);
c = 3e11/ri;

% ind4 = find(dist<=r);
ind5 = find(mesh.nodes(:,dir) >= height(1));
ind4 = find(mesh.nodes(:,dir) <= height(2));
ind6 = find(dist <= r);
index = intersect(intersect(ind4,ind5),ind6);
mesh.mua(index) = mua;
mesh.mus(index) = mus;
mesh.kappa(index) = kappa;
mesh.ri(index) = ri;
mesh.c(index) = c;
mesh.region(index) = reg_label;


% disp(['Number of nodes modified = ' num2str(length(index))]);