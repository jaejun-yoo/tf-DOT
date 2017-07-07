function mesh = add_spherical_target(mesh, position,prop, radius, reg_label)

% Adds spherical (3D) targets to mesh.
% mesh is the mesh variable or filename.
% position = [x y z];
% radius = radius of cylindrical target
% height = [ymin ymax];
% prop = [mua mus ri];
% reg_label = region label (typically > 0)
% mesh = add_spherical_target(mesh, [10 0 50], [0.01 1 1.3], 7.5, 1)

% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = loadmesh(mesh);
end

dist = distance(mesh.nodes(:,1:3),ones(length(mesh.bndvtx),1),[position(1) position(2) position(3)]);

index = find(dist<=radius);

mesh.mua(index) = prop(1);
mesh.mus(index) = prop(2);
mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
mesh.ri(index) = prop(3);
mesh.c(index)=(3e11/prop(3));
mesh.region(index) = reg_label;

disp(['Number of nodes modified = ' num2str(length(index))]);


