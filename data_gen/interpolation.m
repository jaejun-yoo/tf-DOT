clc;
clear;
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
%*************************************************************
% generate meshes for 2d interpolation in detector plane
det_pitch   = 20;    src_pitch_x = 10;    src_pitch_y = 20;
odet_left   = 30;    odet_right  = 170;    odet_top = 105;    odet_bottom = 25;
osrc_left   = 30;    osrc_right  = 180;    osrc_top =  95;    osrc_bottom = 35;
[odetX,odetY]   = meshgrid(odet_left:det_pitch:odet_right,odet_top:-1*det_pitch:odet_bottom);
[osrcX,osrcY]   = meshgrid(osrc_left:src_pitch_x:osrc_right,osrc_top:-1*src_pitch_y:osrc_bottom);

% phantom_left = 0; phantom_right = 200; phantom_top = 130; phantom_bottom = 0;
phantom_left = 30; phantom_right = 170; phantom_top = 105; phantom_bottom = 25;
pad_x        = 64; pad_y         = 32;
[rX,rY]      = meshgrid(linspace(phantom_left,phantom_right,pad_x), linspace(phantom_top,phantom_bottom,pad_y));
z_d = 0; z_s = 40;
G0_rd_rs = cal_G(sqrt(k_sq),D0,odetX',odetY',z_d*ones(size(odetX(:))),osrcX',osrcY',(z_s-1)*ones(size(osrcX(:))));

% load data
% load('imdb_fem_test_3000.mat')

load('imdb_fem_test_inside_det_rand_num_inc_1000_complex.mat')
nbatch = 200;
for ib = 1:nbatch, if images.flag(:,ib) == 0, idx = ib; break;,end,end
U0 = images.data(:,:,idx);
load('imdb_fem_test_inside_det_1_num_inc_height_25_40_200_complex.mat')
% interpolation
% nbatch = 200;
% for ib = 1:nbatch, if images.flag(:,ib) == 0, idx = ib; break;,end,end
% figure,
for ib = nbatch:-1:1
    Ut = images.data(:,:,ib);
%     U0 = images.data(:,:,idx);
    U_rytov = -1*log(Ut./U0);
    Y = G0_rd_rs.*U_rytov;
    data = abs(Y);
%     
%         
%     for z = 1:20,  subplot(5,4,z),
%         imagesc(images.labels(:,:,z,ib)), axis image,
%         title(['ib: ',num2str(ib),', z:',num2str(z)]);
%     end
%     figure,
%     data = abs(images.data(:,:,ib)-images.data(:,:,8));    
    for s = 64:-1:1
        tmp = interp2(odetX,odetY,transpose(reshape(data(:,s),8,5)),rX,rY,'spline');
        images.data_interp(:,:,s,ib) = flip(tmp,1);
%         subplot(121)        
% %         imagesc(images.labels(:,:,15,ib)); axis image;
%         imagesc(flip(reshape(data(:,s),8,5)',1)), axis image;
%         axis image;
%         subplot(122)
%         imagesc(flip(tmp,1)), axis image;
%         pause
    end
    ib
end

% save('imdb_fem_test_inside_det_rand_num_inc_1000_complex_Rytov.mat','images')
save('imdb_fem_test_inside_det_1_num_inc_height_25_40_200_complex_Rytov.mat','images')
% load('imdb_fem_test_inside_det_3_num_inc_height_25_40_200_complex_Rytov.mat')
for ib = 1:200
    tmp = images.labels(:,:,:,ib);
    if tmp(:) ==0
        disp(ib)
    end
end
%*************************************************************
% check scattering image
delx = 64/200;
dely = 32/130;
delz = 20/40;
inc_x = round([30,170]*delx); inc_y = round([35,95]*dely); inc_z = round([1,40]*delz);

for ib = 1:100
    currlabel = images.labels(inc_y(1):inc_y(2),inc_x(1):inc_x(2),inc_z(1):inc_z(2),ib);
    num = sum(currlabel(:)>0)
    if num == 0
        keyboard
    end
    figure(1),
    for z = 1:20,  subplot(5,4,z),
%         imagesc(images.labels(:,:,z,ib)), axis image,
        imagesc(currlabel(:,:,z)), axis image,
        title(['ib: ',num2str(ib),', z:',num2str(z)]);
        drawnow;
    end
    pause
%     homo = images.data_interp(:,:,:,idx);
    scatter = images.data_interp(:,:,:,ib);
    figure,
    for s = 1:64,
        %     subplot(211),
        %     imagesc(images.labels(:,:,mod(s,20)+1,ib)),
        %     axis image,
        %     subplot(212),
        imagesc(scatter(:,:,s)),
        axis image,
        caxis([min(scatter(:)),max(scatter(:))])
        suptitle(num2str(s)),
        %     colorbar
        pause(0.1);
    end    
end
