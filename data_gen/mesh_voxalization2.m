function [pixel,label] = mesh_voxalization2(mesh,numpix)

%% Voxalized the unstructured tetrahedral mesh

% # of voxel 
xmin = min(mesh.nodes(:,1));
xmax = max(mesh.nodes(:,1));
ymin = min(mesh.nodes(:,2));
ymax = max(mesh.nodes(:,2));
zmin = min(mesh.nodes(:,3));
zmax = max(mesh.nodes(:,3));

% Create uniform grid
xstep=(xmax-xmin)/(numpix(1,1)-1);
ystep=(ymax-ymin)/(numpix(1,2)-1);
zstep=(zmax-zmin)/(numpix(1,3)-1);

[X,Y,Z]=meshgrid(xmin:xstep:xmax,...
		 ymin:ystep:ymax,...
		 zmin:zstep:zmax);


pixel.mua = griddata(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.mua,X,Y,Z);
pixel.mus = griddata(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.mus,X,Y,Z);
pixel.ri = griddata(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.ri,X,Y,Z);
pixel.region = griddata(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.region,X,Y,Z);

label = pixel.region;
