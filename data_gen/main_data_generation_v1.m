%%% 20170705
%%% Simulation data generation
%%% for supervised classification
%%% spherical target only
%%% EXAMPLE
    % three lateral target phantom example in y direction
    %meshb = add_cylinderical_target(mesh, [50 0 20], 2.5, [0 65], [0.02 1 1.3],2, 1);
    %meshb = add_cylinderical_target(meshb, [100 0 20], 5, [0 65], [0.02 1 1.3], 2,2);
    %meshb = add_cylinderical_target(meshb, [150 0 20], 10, [0 65], [0.02 1 1.3], 2,3);
    
    %add spherical target in mesh
    %example (repeat it if you need to add multiple target and change the label number)
    %meshb = add_spherical_target(mesh, [100, 65, 20],[0.01,1,1.33], 7.5, 1);
    
clc;
clear;
close all;

addpath('Meshes','Experimental data','Forward model','Utilities')
%*************************************************************
%load fem meshes for forward model
mesh = loadmesh('mesh(200x130x40)');
const_scale = 2.5e-4;
mua = 0.003;
mus = 0.5;
%change optical property of the background (default (0.01,1,1.33))
mesh.mua(:)= mua;
mesh.mus(:)= mus;
mesh.ri(:)=1.5;

%associated variable update 
mesh.kappa = 1./(3*(mesh.mua+mesh.mus));
mesh.c=(3e11./mesh.ri);
f=0.9;
Ro=((mesh.ri-1).^2)./((mesh.ri+1).^2);
cos_theta_c=abs(cos(asin(1./mesh.ri)));
A=((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
mesh.ksi=1./(2*A);

idx_d = (mesh.meas.coord(:,2)==25 | mesh.meas.coord(:,2)==105 | mesh.meas.coord(:,1)==30 | mesh.meas.coord(:,1)==170);
idx_s = (mesh.source.coord(:,1)==180);
idx1 = ismember(mesh.link(:,1), 61:64);
idx2 = ismember(mesh.link(:,2), 19:40);
% idx_d_del = mesh.meas.num(idx_d);
% idx_s_del = mesh.source.num(idx_s);
% idx1 = ismember(mesh.link(:,1), idx_s_del);
% idx2 = ismember(mesh.link(:,2), idx_d_del);
mesh.link((idx1|idx2),:)=[];

mesh.meas.coord(idx_d,:) = [];
mesh.meas.num = mesh.meas.num(1:18);
mesh.meas.int_func(idx_d,:) = [];
mesh.source.coord(idx_s,:) = [];
mesh.source.num = mesh.source.num(1:60);
mesh.source.fwhm(idx_s,:) = [];



%*************************************************************
nbatch = 10;
inc_x = [60:5:140]; inc_y = [45:5:85]; inc_z = [1:1:40];
label_size = [40,26,20]; % x y z, res_x: 5 mm, res_y:5 mm, res_z: 2 mm
mua_inc = linspace(2*mua,5*mua,100);  mus_inc = linspace(2*mus,5*mus,100);
nx = length(inc_x); ny = length(inc_y); nz = length(inc_z);
reg_label = 1;
default_homo = 0;
for ibatch = nbatch:-1:1
    disp(ibatch)
    meshb = mesh;
    num_incl = randperm(3,1);
    position = zeros(3);
    radius = zeros(3,1);
    prop = zeros(3);
    flag = 1;
    for iter = num_incl:-1:1
        inc_cx = inc_x(randi(nx)); inc_cy = inc_y(randi(ny)); inc_cz = inc_z(randi(nz));
        position(iter,:) = [inc_cx, inc_cy, inc_cz];
        radius(iter) = randperm(10,1)+5;
        cmua_inc = mua_inc(randi(100)); cmus_inc = mus_inc(randi(100));
        prop(iter,:) = [cmua_inc, cmus_inc, 1.33];
        [meshb] = add_spherical_target(meshb, position(iter,:),prop(iter,:), radius(iter,:), reg_label);
    end
    
    %*************************************************************
    % forward model data generation
    %data.phi is complex data inside medium when sources are at sources
    %position. %data.aphi is adjoint data inside medium when sources are at detector
    %position. %data.complex is the boundary data for source and detectors
    %pair
    %use meshb when to genrate hetrogenous data
    %use mesh when to generate homogenous data
    
    if default_homo == 0 && sum(meshb.region(:)) == 0
        data=forward_data(meshb,70); % set 0 for frequency if you need only amplitude data
        data.complex=reshape(data.complex,18,60); %for boundary data
        homo_data = data.complex/const_scale;
        default_homo = 1
    end
    if sum(meshb.region(:)) == 0                
        images.data(:,:,ibatch) = homo_data;
        images.position(:,:,ibatch) = zeros(3);
        images.radius(:,ibatch) = zeros(3,1);
        images.prop(:,:,ibatch) = zeros(3);
        images.labels(:,:,:,ibatch) =  zeros(label_size(2),label_size(1),label_size(3));
        images.num_incl(ibatch) = 0;
        images.flag(:,ibatch) = 0;
        images.meshb = meshb;
    else
        [pixel,label] = mesh_voxalization2(meshb,label_size);
        data=forward_data(meshb,70); % set 0 for frequency if you need only amplitude data
        data.complex=reshape(data.complex,18,60); %for boundary data
        images.data(:,:,ibatch) = data.complex/const_scale;
        images.position(:,:,ibatch) = position;
        images.radius(:,ibatch) = radius;
        images.prop(:,:,ibatch) = prop;
        images.labels(:,:,:,ibatch) = label;
        images.num_incl(ibatch) = num_incl;
        images.flag(:,ibatch) = flag;
        images.meshb{ibatch} = meshb;
    end
end

% figure, for z = 1:20, subplot(5,4,z),imagesc(label(:,:,z)), axis image, title(num2str(z)), end

if default_homo == 0
    for ibatch = 1
        disp(ibatch)
        meshb = mesh;
        num_incl = 0;
        position = zeros(3);
        radius = zeros(3,1);
        prop = zeros(3);        
        [pixel,label] = mesh_voxalization2(meshb,label_size);
        data=forward_data(meshb,70); % set 0 for frequency if you need only amplitude data
        data.complex=reshape(data.complex,18,60); %for boundary data
        homo_data = data.complex/const_scale;
        default_homo = 1
%         images.data(:,:,ibatch) = homo_data;
%         images.position(:,:,ibatch) = zeros(3);
%         images.radius(:,ibatch) = zeros(3,1);
%         images.prop(:,:,ibatch) = zeros(3);
%         images.labels(:,:,:,ibatch) =  zeros(label_size(2),label_size(1),label_size(3));
%         images.num_incl(ibatch) = 0;
%         images.flag(:,ibatch) = 0;
%         images.meshb{ibatch} = meshb;
    end
end

%% INTERPOLATION IN BOTH DET AND SRC
% generate meshes for 2d interpolation in detector and source planes
det_pitch   = 20;    src_pitch_x = 10;    src_pitch_y = 20;
odet_left   = 50;    odet_right  = 150;    odet_top = 85;    odet_bottom = 45;
osrc_left   = 30;    osrc_right  = 170;    osrc_top =  95;    osrc_bottom = 35;
[odetX,odetY]   = meshgrid(odet_left:det_pitch:odet_right,odet_top:-1*det_pitch:odet_bottom);
[osrcX,osrcY]   = meshgrid(osrc_left:src_pitch_x:osrc_right,osrc_top:-1*src_pitch_y:osrc_bottom);
pad_dx         = 21; pad_dy         = 9; pad_sx         = 29; pad_sy          = 4; z_d         = 0; z_s         = 39;
[rdX,rdY]      = meshgrid(linspace(odet_left,odet_right,pad_dx), linspace(odet_top,odet_bottom,pad_dy));
[rsX,rsY]      = meshgrid(linspace(osrc_left,osrc_right,pad_sx), linspace(osrc_top,osrc_bottom,pad_sy));

% *** optical parameters *** %
% Create the background parameters
frequency = 70*1e6; 
mua_bkg = 0.003;  % mm
mus_bkg = 0.5; % mm
ref_bkg = 1.5;
c       = (3e11./ref_bkg); % speed of light in medium
omega   = 2*pi*frequency; %*1e6; % modulation frequency
D0      = 1/3/mus_bkg;
k_sq    = (mua_bkg-1i*omega/c)/D0; % k^2 ; diffuse wave number
G0_rd_rs = cal_G(sqrt(k_sq),D0,odetX',odetY',z_d*ones(size(odetX(:))),osrcX',osrcY',(z_s-1)*ones(size(osrcX(:))));

% for ib = 1:nbatch, if images.flag(:,ib) == 0, idx = ib; break; end,end
U0 = homo_data;

for ib = nbatch:-1:1
    Ut = images.data(:,:,ib);
    U_rytov = -1*log(Ut./U0);
    Y = G0_rd_rs.*U_rytov;
    data = abs(Y);  
    figure(1), for z = 1:20, subplot(5,4,z),imagesc(images.labels(:,:,z,ib)), axis image, title(num2str(z)), end
    drawnow
%     figure(2),
    for s = 60:-1:1
        tmp = interp2(odetX,odetY,transpose(reshape(data(:,s),6,3)),rdX,rdY,'spline');
        tmp2 = flip(tmp,1);
        data_det_interp(:,s,ib) = tmp2(:);
%         subplot(121)        
%         imagesc(flip(reshape(data(:,s),6,3)',1)), axis image;
%         axis image;
%         subplot(122)
%         imagesc(tmp2), axis image;
%         pause
    end
%     figure(3)
    for d = 189:-1:1
        tmp = interp2(osrcX,osrcY,transpose(reshape(data_det_interp(d,:,ib),15,4)),rsX,rsY,'spline');
        tmp2 = flip(tmp,1);
        images.data_interp(d,:,ib) = tmp2(:);
%         subplot(121)        
%         imagesc(flip(reshape(data_det_interp(d,:,ib),15,4)',1)), axis image;
%         axis image;
%         subplot(122)
%         imagesc(flip(tmp,1)), axis image;
%         pause
    end
    ib
end
tmp_data = repmat(reshape(images.data_interp,[],4,29,nbatch),1,1,1,189);

new_data = zeros(9*2,189*nbatch);
new_label = permute(reshape(images.labels(9:17,10:30,:,:),[],20,nbatch),[2,1,3]);
new_label = reshape(new_label,[],nbatch*189);
slayer = {1:2,2:3,3:4};
for ib = 1:nbatch    
    for di = 1:189
        sl = mod(di-1,3)+1;
        ti = mod(ceil(di/3)-1,21)+1;        
        new_data(:,ib) = reshape(tmp_data(di,slayer{sl},ti:ti+8,ib),[],1);        
    end
end
new_nbatch = size(new_label,2);

images.new_data = new_data;
images.new_label = new_label;
images.set = ones(1,new_nbatch);
images.set((floor(2*nbatch/3)+1)*189+1:end) = 3;
save(['imdb_fem_inside_det_num_inc_',num2str(new_nbatch),'_s116_d189_interp.mat'],'images')















