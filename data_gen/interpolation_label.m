
phantom_left = 0; phantom_right = 200; phantom_top = 130; phantom_bottom = 0;
pad_x        = 64; pad_y         = 32;
[rX,rY]      = meshgrid(linspace(phantom_left,phantom_right,pad_x), linspace(phantom_top,phantom_bottom,pad_y));
[rX2,rY2]      = meshgrid(linspace(phantom_left,phantom_right,4*pad_x), linspace(phantom_top,phantom_bottom,4*pad_y));

load('imdb_fem_test_inside_det_rand_and_3_num_inc_1200_complex_Rytov_cleansed.mat')
nbatch = 200;
for ib = 1:nbatch, if images.flag(:,ib) == 0, idx = ib; break;,end,end

for ib = 1:size(images.labels,4)%:-1:1
%     figure(1),
%     for z = 1:20,  subplot(5,4,z),
%         imagesc(images.labels(:,:,z,ib)), axis image,
%         title(['ib: ',num2str(ib),', z:',num2str(z)]);
%     end
    
    for z =1:20,  subplot(5,4,z),  
%         figure(2),
        tmp = interp2(rX,rY,single(images.labels(:,:,z,ib)),rX2,rY2,'cubic');        
%         images.labels_interp(:,:,z,ib) = tmp;
%         subplot(121)
%         %         imagesc(images.labels(:,:,15,ib)); axis image;
%         imagesc(images.labels(:,:,z,ib)), axis image;
%         axis image;
%         subplot(122)
%         imagesc(single(tmp(1:4:end,1:4:end,1:4:end))), axis image;
%         suptitle(num2str(z))
%         pause        
        images.labels(:,:,z,ib) = tmp(1:4:end,1:4:end,1:4:end);
    end
    ib
end
images.labels = (images.labels-min(images.labels(:)))/(max(images.labels(:))-min(images.labels(:)));
save('imdb_fem_test_inside_det_rand_and_3_num_inc_1200_complex_Rytov_cleansed_label_smooth2.mat','images')

figure, for ib = 585:-1:1,tmp = images.labels(:,:,end,ib); if sum(tmp(:))>0, ib, imagesc(tmp),title(ib), axis image, pause, end,end