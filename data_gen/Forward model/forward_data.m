function [data]=forward_data(mesh,frequency)

% Calculates data (phase and amplitude) given mesh and frequency of modulation
% data.phi forward model data at each unknown basis (fem nodes)
% data.aphi forward model data at unknown basis when the sources are placed at detector position (adjoint sources)
% data.complex forward model boundary data (sources and detector pairs) in complex form 

if frequency < 0
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = loadmesh(mesh);
end

% modulation frequency
omega = 2*pi*frequency*1e6;

% Create FEM matricex
if mesh.dimension == 2
    [i,j,s] = MassMatrix_2d(mesh.nodes(:,1:2),sort(mesh.elements')',mesh.bndvtx,mesh.mua,mesh.kappa,mesh.ksi,mesh.c,omega);
elseif mesh.dimension ==3
    [i,j,s] = MassMatrix_3d(mesh.nodes,sort(mesh.elements')',mesh.bndvtx,mesh.mua,mesh.kappa,mesh.ksi,mesh.c,omega);
end

junk = length(find(i==0));
MASS = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
clear junk i j s omega

% Calculate the RHS 

source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);
%% JJ changed
tmp = mesh.link(:,3)==0;
mesh.link(tmp,3) = 1;
%%
% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);

%[nsource,junk]=size(source);
qvec = spalloc(nnodes,nsource,nsource*100);
if mesh.dimension == 2
    for i = 1 : nsource
        s_ind = mesh.source.num == source(i);
        if mesh.source.fwhm(s_ind) == 0
            qvec(:,i) = Source_point(mesh,mesh.source.coord(s_ind,1:2));
        else
            qvec(:,i) = Source(mesh.nodes(:,1:2),sort(mesh.elements')',mesh.dimension,mesh.source.coord(s_ind,1:2),mesh.source.fwhm(s_ind));
        end
    end
elseif mesh.dimension == 3
    for i = 1 : nsource
        s_ind = mesh.source.num == source(i);
        if mesh.source.fwhm(s_ind) == 0
            qvec(:,i) = Source_point(mesh,mesh.source.coord(s_ind,1:3));
        else
            qvec(:,i) = Source(mesh.nodes,sort(mesh.elements')',mesh.dimension,mesh.source.coord(s_ind,:),mesh.source.fwhm(s_ind));
        end
    end
end
clear junk i nnodes nsource w;

% Catch zero frequency (CW) here
if frequency == 0
    MASS = real(MASS);
    qvec = real(qvec);
end

% catch error in source vector
junk = sum(qvec);
junk = find(junk==0);
if ~isempty(junk)
    display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
end
clear junk

% Calculate field for all sources
[data.phi,mesh.R]=solver(MASS,mesh,qvec);
clear qvec;

% Now calculate Adjoint source vector
[qvec] = source_adjoint(mesh);

% Catch zero frequency (CW) here
if frequency == 0
    qvec = real(qvec);
end

% Calculate adjoint field for all detectors
[data.aphi]=solver(conj(MASS),mesh,conj(qvec));
clear qvec MASS;

% Calculate boundary data

[data.complex]=boundary_data(mesh,data.phi);


%%%%%%%%%%%%% sub function space %%%%%%%%%%%%%%%%%%%%%%

function [phi,R]=solver(Mass,mesh,qvec)

% Calculates the field. Uses the biconjugate gradient stabilized method

[nnodes,nsource]=size(qvec);
msg=[];
flag = 0;
phi=zeros(nnodes,nsource);

if isfield(mesh,'R') == 0
	% keeping this out until 
	%alpha = max(sum(abs(Mass),2)./diag(abs(Mass)))-2;
	%alpha =.1;
	%R=ichol(Mass, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
	R = ichol(Mass);
   % R = cholinc(Mass,1e-3);
	%[L,U] = ilu(Mass,struct('type','ilutp','droptol',1e-3,'thresh',0));
	%R = diag(sqrt(abs(diag(U))))\U;

    mesh.R = R;
else
    R = mesh.R;
end

for i = 1 : nsource
    [x,flag] = bicgstabl(Mass,qvec(:,i),1e-12,100,R',R);
    msg = [msg flag];
    phi(:,i) = x;
end

if any(msg==1)
    error('Some solutions did not converge');
elseif any(msg==2)
    error('Some solutions are unusable');
elseif any(msg==3)
    error('Some solutions from stagnated iterations are unusable');

end

function [data] = boundary_data(mesh,phi)

% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));

% source = unique(mesh.link(:,1));

if isfield(mesh.meas,'int_func') == 0
    error('Need to call move_detector on the mesh first');
else
    
    % We don't want contributions from internal nodes on boundary values
	
    bnd_val = mesh.bndvtx(mesh.elements(mesh.meas.int_func(:,1),:));
    [nrow,ncol]=size(bnd_val);
    for i = 1 : nrow
        for j = 1 : ncol
            if bnd_val(i,j) == 0
                mesh.meas.int_func(i,j+1) = 0;
                % make sure the integral is unity here!
                mesh.meas.int_func(i,2:end) = ...
                    1./sum(mesh.meas.int_func(i,2:end)) .* ...
                    mesh.meas.int_func(i,2:end);
            end
        end
    end
    
    data = NaN(size(mesh.link(:,1)));
     
	for i = 1:size(mesh.link,1)
		if mesh.link(i,3) == 1
			sn = source == mesh.link(i,1);
			dn = mesh.meas.num == mesh.link(i,2);
			vtx_ind = mesh.elements(mesh.meas.int_func(dn,1),:);
			data(i) = mesh.meas.int_func(dn,2:end)*phi(vtx_ind,sn);
		end
	end
    
end

function qvec = Source_point(mesh,source)

% Allocate memory
[nnodes,junk]=size(mesh.nodes);
qvec = spalloc(nnodes,1,4);

% find elements for point sources and calculate interpolation functions
if mesh.dimension == 2
    [ind,int_func] = mytsearchn(mesh,source(:,1:2));
elseif mesh.dimension == 3
    [ind,int_func] = mytsearchn(mesh,source);
end
int_func = [ind int_func];

% Go through all measurements and integrate to get nodal values
qvec(mesh.elements(int_func(1,1),:),1) = int_func(1,2:end)'.*complex(cos(0.15),sin(0.15));


function qvec = source_adjoint(mesh)

% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
meas = unique(foo(:,2));
nmeas = length(meas);

[nnodes,junk]=size(mesh.nodes);
qvec = spalloc(nnodes,nmeas,nmeas*5);

% Go through all measurements and integrate to get nodal values
for i = 1 : nmeas
  qvec(mesh.elements(mesh.meas.int_func(mesh.meas.num == meas(i),1),:),i) = ...
      mesh.meas.int_func(mesh.meas.num == meas(i),2:end)' .* ...
      complex(cos(0.15),sin(0.15));
end

