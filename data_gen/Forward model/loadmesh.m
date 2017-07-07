function mesh = loadmesh(fn)
% fn is the filename of the mesh (with no extension)

mesh.name = fn;

%% Read mesh nodes
if exist([fn '.node']) == 0
    error([fn '.node file is not present']);
elseif exist([fn '.node']) == 2
    mesh.nodes = load(strcat(fn, '.node'));
    mesh.bndvtx = mesh.nodes(:,1); %sets 1 if boundary node, 0 if internal
    mesh.nodes = mesh.nodes(:,2:end);
end

%% Read mesh parameters

if exist([fn '.param']) == 0
     error('.param file is not present');
elseif exist([fn '.param']) == 2
    param = importdata([fn '.param']);
    if isfield(param,'textdata') == 1
		mesh.type = 'stnd';
		param = param.data;
		mesh.mua = param(:,1);
		mesh.kappa = param(:,2);
		mesh.ri = param(:,3);
		mesh.mus = ((1./mesh.kappa)./3)-mesh.mua; 
    end
end

%% Read mesh element
if exist([fn '.elem']) == 0
    error('.elem file is not present');
elseif exist([fn '.elem']) == 2
    mesh.elements = load(strcat(fn, '.elem'));
    [junk,dim]=size(mesh.elements);
    mesh.dimension = dim-1;
    if mesh.dimension == 2
        mesh.nodes(:,3) = 0;
    end
end

%% Region file
if exist([fn '.region']) ~= 0
    mesh.region = load(strcat(fn, '.region'));
elseif exist([fn '.region']) ~= 2
    mesh.region = zeros(length(mesh.nodes),1);
end


%% Load source locations
if exist([fn '.source']) == 0
    disp([fn '.source file is not present']);
elseif exist([fn '.source']) == 2
    mesh.source.distributed = 0;
    source = importdata([fn '.source']);
  
	mesh.source.fixed = 0;
  
  % Text flags at top of source file with column headings 
	if sum(sum(strcmp(source.textdata,'fixed'))) == 1
		mesh.source.fixed = 1;
		source.textdata = source.textdata(2:end,:);
	end

	mesh.source.num = source.data(:,logical(strcmp(source.textdata,'num')));
	mesh.source.coord(:,1) = source.data(:,logical(strcmp(source.textdata,'x')));
	mesh.source.coord(:,2) = source.data(:,logical(strcmp(source.textdata,'y')));

	if sum(strcmp(source.textdata,'z')) == 1
		mesh.source.coord(:,3) = source.data(:,logical(strcmp(source.textdata,'z')));
	end
	if sum(strcmp(source.textdata,'fwhm')) == 1
		mesh.source.fwhm = source.data(:,logical(strcmp(source.textdata,'fwhm')));
	else
		mesh.source.fwhm = zeros(ns,1);
	end

    % Check and poistion sources
    if mesh.source.fixed == 1;
        disp('Fixed Sources');
        if mesh.dimension == 2
            [ind,int_func] = mytsearchn(mesh,mesh.source.coord(:,1:2));
        elseif mesh.dimension == 3
            [ind,int_func] = mytsearchn(mesh,mesh.source.coord);
        end
        if any(isnan(ind)) == 1
            error('Source(s) outside the mesh; either move them manually or remove ''fixed'' from the source file');
        end
    elseif mesh.source.fixed == 0;
		mus_eff = mesh.mus;
		disp('Moving Sources');
        [mesh]=moving_source(mesh,mus_eff,3);
        clear source mus_eff
    end
end

%% Load detector locations
if exist([fn '.meas']) == 0
    disp([fn '.meas file is not present']);
elseif exist([fn '.meas']) == 2
    meas = importdata([fn '.meas']);
	
	% Text flags at top of source file with column headings 
	mesh.meas.fixed = 0;
	if sum(sum(strcmp(meas.textdata,'fixed'))) == 1
		mesh.meas.fixed = 1;
		meas.textdata = meas.textdata(2,:);
	end
	mesh.meas.num = meas.data(:,logical(strcmp(meas.textdata,'num')));
	mesh.meas.coord(:,1) = meas.data(:,logical(strcmp(meas.textdata,'x')));
	mesh.meas.coord(:,2) = meas.data(:,logical(strcmp(meas.textdata,'y')));
	if sum(strcmp(meas.textdata,'z')) == 1
		mesh.meas.coord(:,3) = meas.data(:,logical(strcmp(meas.textdata,'z')));
	end
    
    
    % Check and position detectors
    if mesh.meas.fixed == 0
        disp('Moving Detectors');
        [mesh]=moving_detector(mesh);
    elseif mesh.meas.fixed == 1
        disp('Fixed Detectors');
    end
    
    if mesh.dimension == 2
        [ind,int_func] = mytsearchn(mesh,mesh.meas.coord(:,1:2));
    elseif mesh.dimension == 3
        [ind,int_func] = mytsearchn(mesh,mesh.meas.coord);
    end
    if any(isnan(ind)) == 1
        warning('Detector(s) outside the mesh');
    else
        mesh.meas.int_func = [ind int_func];
    end
    clear meas
end


%% Load link list for source and detector
if exist([fn '.link']) == 0
    disp([fn '.link file is not present']);
elseif exist([fn '.link']) == 2
 	link = importdata([fn '.link']);
	mesh.link = link.data;
end

%% Load identidity list if exists for the internal RI boundary nodes
if exist([fn '.ident']) == 2
    mesh.ident = load(strcat(fn, '.ident'));
end


%% speed of light in medium
mesh.c=(3e11./mesh.ri);

%% Set boundary coefficient using definition of A using the Fresenel's law:
f=0.9;
Ro=((mesh.ri-1).^2)./((mesh.ri+1).^2);
thetac=asin(1./mesh.ri);
cos_theta_c=abs(cos(asin(1./mesh.ri)));
A=((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
mesh.ksi=1./(2*A);


%%%%%%%%%%%%%%%%%%% Sub function space %%%%%%%%%%%%%%%%%%%%%%

function mesh = moving_detector(mesh)
% Moves the detectors onto the surface of the mesh

remove_last = 0;
if size(mesh.meas.coord,2) == 2
    mesh.meas.coord(:,end+1) = 0;
    remove_last = 1;
end

%% get list of boundary faces
if size(mesh.elements,2) == 4
    faces = [mesh.elements(:,[1,2,3]);
              mesh.elements(:,[1,2,4]);
              mesh.elements(:,[1,3,4]);
              mesh.elements(:,[2,3,4])];
    faces = sort(faces,2);
    faces = unique(faces,'rows');
    faces = faces(sum(mesh.bndvtx(faces),2)==3,:);
elseif size(mesh.elements,2) == 3
    if mesh.dimension == 3
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==3,:);
    elseif mesh.dimension == 2
        faces = [mesh.elements(:,[1,2]);
                  mesh.elements(:,[1,3]);
                  mesh.elements(:,[2,3])];
        faces = sort(faces,2);
        faces = unique(faces,'rows');
        faces = faces(sum(mesh.bndvtx(faces),2)==2,:);
    end
end

%% loop through detectors
for i=1:size(mesh.meas.coord,1)
    
    if mesh.dimension == 2
        
        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.meas.coord(i,:) 0]);
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);
        
        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointLineDistance(mesh.nodes(fi(j,1),:), ...
                mesh.nodes(fi(j,2),:),mesh.meas.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % move detector to the closest point on that face
        mesh.meas.coord(i,1:2) = point(smallest(1),1:2);
        
    elseif mesh.dimension == 3

        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.meas.coord(i,:));
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);

        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
                mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.meas.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % move detector to the closest point on that face
        mesh.meas.coord(i,1:3) = point(smallest(1),1:3);
    
    end
        
end

if remove_last
    mesh.meas.coord(:,end) = [];
end


function mesh = moving_source(mesh,mus_eff,w)
% Moves the sources inside the mesh by 1 scattering distance

failed = 0;

remove_last = 0;
if size(mesh.source.coord,2) == 2
    mesh.source.coord(:,end+1) = 0;
    remove_last = 1;
end

%% check if mus_eff is unrealistic for the mesh size
scatt_dist = 1/mean(mus_eff);
xd = max(mesh.nodes(:,1)) - min(mesh.nodes(:,1));
yd = max(mesh.nodes(:,2)) - min(mesh.nodes(:,2));
if scatt_dist*10 > min(xd,yd)
    scatt_dist = 1;
    ho = findobj('type','figure','name','Warning - Small Mesh');
    if ~isempty(ho)
        close(ho);
    end
    errordlg('Mesh is too small for moving, 1mm will be used for as moving distance','Warning - Small Mesh');
end

%% get list of boundary faces
out_normal = 0;
if size(mesh.elements,2) == 4
    faces = [mesh.elements(:,[1,2,3]);
              mesh.elements(:,[1,2,4]);
              mesh.elements(:,[1,3,4]);
              mesh.elements(:,[2,3,4])];
    faces = sort(faces,2);
    faces = unique(faces,'rows');
    faces = faces(sum(mesh.bndvtx(faces),2)==3,:);
elseif size(mesh.elements,2) == 3
    if mesh.dimension == 3
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==3,:);
        out_normal = 1;
    elseif mesh.dimension == 2
        faces = [mesh.elements(:,[1,2]);
                  mesh.elements(:,[1,3]);
                  mesh.elements(:,[2,3])];
        faces = sort(faces,2);
        faces = unique(faces,'rows');
        faces = faces(sum(mesh.bndvtx(faces),2)==2,:);
    end
end

%% loop through sources
for i=1:size(mesh.source.coord,1)
    
    if mesh.dimension == 2
        
        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,:) 0]);
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);
        
        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointLineDistance(mesh.nodes(fi(j,1),:), ...
                mesh.nodes(fi(j,2),:),mesh.source.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % find normal of that face
        a = mesh.nodes(fi(smallest(1),1),:);
        b = mesh.nodes(fi(smallest(1),2),:);
        n = [a(2)-b(2) b(1)-a(1)];
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos1 = point(smallest(1),1:2) + n*scatt_dist;
        pos2 = point(smallest(1),1:2) - n*scatt_dist;
        ind = mytsearchn(mesh,[pos1;pos2]);
        if ~isnan(ind(1))
            mesh.source.coord(i,1) = pos1(1);
            mesh.source.coord(i,2) = pos1(2);
        elseif ~isnan(ind(2))
            mesh.source.coord(i,1) = pos2(1);
            mesh.source.coord(i,2) = pos2(2);
        else
            failed = 1;
        end
        
    elseif mesh.dimension == 3

        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.source.coord(i,:));
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);

        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
                mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.source.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % find normal of that face
        a = mesh.nodes(fi(smallest(1),1),:);
        b = mesh.nodes(fi(smallest(1),2),:);
        c = mesh.nodes(fi(smallest(1),3),:);
        n = cross(b-a,c-a);
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos2 = point(smallest(1),:) + n*scatt_dist;
        pos1 = point(smallest(1),:) - n*scatt_dist;
        if ~out_normal
            ind = mytsearchn(mesh,[pos1;pos2]);
        end
        if out_normal || ~isnan(ind(1))
            mesh.source.coord(i,1) = pos1(1);
            mesh.source.coord(i,2) = pos1(2);
            mesh.source.coord(i,3) = pos1(3);
        elseif ~isnan(ind(2))
            mesh.source.coord(i,1) = pos2(1);
            mesh.source.coord(i,2) = pos2(2);
            mesh.source.coord(i,3) = pos2(3);
        else
            failed = 1;
        end
    
    end
        
end

if failed == 1
    errordlg('Source(s) could not be moved. The mesh structure may be poor.','Warning');
end

if remove_last
    mesh.source.coord(:,end) = [];
end




function [dist,point] = pointLineDistance(A,B,p)
% finds closest point and distance between a line segment and point

t = dot(p - A,B - A) / dot(B - A,B - A);
if t < 0
    t = 0;
elseif t > 1
    t = 1;
end
point = A + (B - A)*t;

dist = norm(p - point);


function [dist,PP0] = pointTriangleDistance(TRI,P)
% calculate distance between a point and a triangle in 3D

% Copyright (c) 2009, Gwendolyn Fischer
% All rights reserved.

% rewrite triangle in normal form
B = TRI(1,:);
E0 = TRI(2,:)-B;
E1 = TRI(3,:)-B;


D = B - P;
a = dot(E0,E0);
b = dot(E0,E1);
c = dot(E1,E1);
d = dot(E0,D);
e = dot(E1,D);
f = dot(D,D);

det = a*c - b*b;
s   = b*e - c*d;
t   = b*d - a*e;

% Terible tree of conditionals to determine in which region of the diagram
% shown above the projection of the point into the triangle-plane lies.
if (s+t) <= det
  if s < 0
    if t < 0
      %region4
      if (d < 0)
        t = 0;
        if (-d >= a)
          s = 1;
        else
          s = -d/a;
        end
      else
        s = 0;
        if (e >= 0)
          t = 0;
        else
          if (-e >= c)
            t = 1;
          else
            t = -e/c;
          end
        end
      end %of region 4
    else
      % region 3
      s = 0;
      if e >= 0
        t = 0;
      else
        if -e >= c
          t = 1;
        else
          t = -e/c;
        end
      end
    end %of region 3 
  else
    if t < 0
      % region 5
      t = 0;
      if d >= 0
        s = 0;
      else
        if -d >= a
          s = 1;
        else
          s = -d/a;
        end
      end
    else
      % region 0
      invDet = 1/det;
      s = s*invDet;
      t = t*invDet;
    end
  end
else
  if s < 0
    % region 2
    tmp0 = b + d;
    tmp1 = c + e;
    if tmp1 > tmp0 % minimum on edge s+t=1
      numer = tmp1 - tmp0;
      denom = a - 2*b + c;
      if numer >= denom
        s = 1;
        t = 0;
      else
        s = numer/denom;
        t = 1-s;
      end
    else          % minimum on edge s=0
      s = 0;
      if tmp1 <= 0
        t = 1;
      else
        if e >= 0
          t = 0;
        else
          t = -e/c;
        end
      end
    end %of region 2
  else
    if t < 0
      %region6 
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
        numer = tmp1 - tmp0;
        denom = a-2*b+c;
        if (numer >= denom)
          t = 1;
          s = 0;
        else
          t = numer/denom;
          s = 1 - t;
        end
      else  
        t = 0;
        if (tmp1 <= 0)
            s = 1;
        else
          if (d >= 0)
              s = 0;
          else
              s = -d/a;
          end
        end
      end
      %end region 6
    else
      % region 1
      numer = c + e - b - d;
      if numer <= 0
        s = 0;
        t = 1;
      else
        denom = a - 2*b + c;
        if numer >= denom
          s = 1;
          t = 0;
        else
          s = numer/denom;
          t = 1-s;
        end
      end %of region 1
    end
  end
end

PP0 = B + s*E0 + t*E1;

dist = norm(P - PP0);