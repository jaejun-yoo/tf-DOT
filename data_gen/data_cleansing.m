clear
load('imdb_fem_test_inside_det_1_num_inc_height_25_40_200_complex_Rytov.mat')
delx = 64/200;
dely = 32/130;
delz = 20/40;
inc_x = round([40,160]*delx); inc_y = round([40,90]*dely); inc_z = round([1,40]*delz);

% figure,
count = 1;
for ib=size(images.labels,4) :-1:1
    currlabel = images.labels(inc_y(1):inc_y(2),inc_x(1):inc_x(2),inc_z(1):inc_z(2),ib);
    num = sum(currlabel(:)>0);
    if ((num <= 100) && (sum(images.flag(:,ib)) ~= 0))    
        idx_ib(count) = ib;
        count = count +1
        
    end
%     for z = 1:20,         
%         subplot(211),imagesc( images.labels(:,:,z,ib)), axis image,
%         subplot(212),imagesc(currlabel(:,:,z)), axis image, 
%         suptitle(['ib: ',num2str(ib),', z:',num2str(z),' ,num: ',num2str(num)]),pause(0.2);
%     end
end


images.data(:,:,idx_ib) = [];
images.position(:,:,idx_ib) = [];
images.radius(:,idx_ib) = [];
images.height(:,idx_ib) = [];
images.prop(:,:,idx_ib) = [];
images.dir(idx_ib) = [];
images.labels(:,:,:,idx_ib) = [];
images.num_incl(idx_ib) = [];
images.flag(:,idx_ib) = [];
images.set(idx_ib) = [];
images.data_interp(:,:,:,idx_ib) = [];

% for ib = 1:200
%     if images.flag(:,ib) ~= 0      
%         currlabel = images.labels(inc_y(1):inc_y(2),inc_x(1):inc_x(2),inc_z(1):inc_z(2),ib);    
%         for z = 1:20,
%             subplot(211),imagesc( images.labels(:,:,z,ib)), axis image,
%             subplot(212),imagesc(currlabel(:,:,z)), axis image,
%             suptitle(['ib: ',num2str(ib),', z:',num2str(z),' ,num: ',num2str(num)]),pause(0.2);
%         end
%     end
% end
save('imdb_fem_test_inside_det_1_num_inc_height_25_40_200_complex_Rytov_cleansed.mat','images')

%%
images2 = images;
load('imdb_fem_test_inside_det_rand_and_3_num_inc_1200_complex_Rytov_cleansed.mat')
images.data = cat(3,images.data,images2.data);
images.position = cat(3,images.position,images2.position);
images.radius = cat(2, images.radius,images2.radius);
images.height = cat(2, images.height,images2.height);
images.prop = cat(3,images.prop,images2.prop);
images.dir = cat(2,images.dir,images2.dir);
images.labels = cat(4,images.labels,images2.labels);
images.num_incl = cat(2, images.num_incl,images2.num_incl);
images.flag = cat(2, images.flag,images2.flag);
images.set = cat(2, images.set,images2.set);
images.data_interp = cat(4,images.data_interp,images2.data_interp);
% delx = 64/200;
% dely = 32/130;
% delz = 20/40;
% inc_x = round([26,173]*delx); inc_y = round([34,94]*dely); inc_z = round([1,40]*delz);
% for ib = size(images.labels,4):-1:1   
%     images.newlabels(:,:,:,ib) = images.labels(inc_y(1):inc_y(2),inc_x(1):inc_x(2),inc_z(1):inc_z(2),ib);
%     ib
% end
% load('imdb_fem_test_inside_det_rand_and_3_num_inc_1200_complex_Rytov_cleansed.mat')


save('imdb_fem_test_inside_det_rand_and_3_num_inc_1400_complex_Rytov_cleansed.mat','images')