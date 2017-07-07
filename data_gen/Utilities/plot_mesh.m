function plot_mesh(mesh,plane,level,numpix)

% mesh is the mesh variable or filename.
% plane must be either 'x', 'y' or 'z'
% level must be within mesh volume
% numpix number of pixel in x,y,z direction


% load mesh
if ischar(mesh)== 1
    mesh = loadmesh(mesh);
end

x = mesh.nodes(:,1);
y = mesh.nodes(:,2);
z = mesh.nodes(:,3);

xmax = max(x);
xmin = min(x);
xstep=(xmax-xmin)/numpix(1);
ymax = max(y);
ymin = min(y);
ystep=(ymax-ymin)/numpix(2);
zmax = max(z);
zmin = min(z);
zstep=(zmax-zmin)/numpix(3);

if plane == 'x'
  [y,z]=meshgrid(ymin:ystep:ymax,zmin:zstep:zmax);
  y = reshape(y,numel(y),1);
  z = reshape(z,numel(z),1);
  mesh_2d.nodes = [y y z];
  mesh_2d.nodes(:,1) = level;
  mesh_2d.elements = delaunayn([y z]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [y z];
  mesh_2d.nodes(:,3) = 0;
elseif plane == 'y'
  [x,z]=meshgrid(xmin:xstep:xmax,zmin:zstep:zmax);
  x = reshape(x,numel(x),1);
  z = reshape(z,numel(z),1);
  mesh_2d.nodes = [x z z];
  mesh_2d.nodes(:,2) = level;
  mesh_2d.elements = delaunayn([x z]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [x z];
  mesh_2d.nodes(:,3) = 0;
elseif plane == 'z'
  [x,y]=meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
  x = reshape(x,numel(x),1);
  y = reshape(y,numel(y),1);
  mesh_2d.nodes = [x y];
  mesh_2d.nodes(:,3) = level;
  mesh_2d.elements = delaunayn([x y]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [x y];
  mesh_2d.nodes(:,3) = 0;
else
  display('define plane is not valid');
  return;
end

mesh_2d.region = interpolatef2r(mesh,mesh_2d,mesh.region);

figure;
h = trisurf(mesh_2d.elements,...
	mesh_2d.nodes(:,1),...
	mesh_2d.nodes(:,2),...
	mesh_2d.nodes(:,3),...
	mesh_2d.region);

view(2);
shading interp;
axis equal; 
axis off;


function [val_int] = interpolatef2r(fwd_mesh,recon_mesh,val)

% This function interpolates fwd_mesh into recon_mesh
% For the Jacobian it is an integration!
NNC = size(recon_mesh.nodes,1);
for i = 1 : NNC
  if fwd_mesh.fine2coarse(i,1) ~= 0
    if isnan(fwd_mesh.fine2coarse(i,1)) == 1
      val_int(i,1) = NaN;
    else
      val_int(i,1) = (fwd_mesh.fine2coarse(i,2:end) * val(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
    end
  elseif fwd_mesh.fine2coarse(i,1) == 0
    dist = distance(fwd_mesh.nodes, fwd_mesh.bndvtx,[recon_mesh.nodes(i,1:2) 0]);
    mindist = find(dist==min(dist));
    mindist = mindist(1);
    val_int(i,1) = val(mindist);
  end
end

