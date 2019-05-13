% define functions
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
deconv_tik = @(y,psf,miu) crop2d(real(Ft(F(pad2d(y)).*...
    conj(F(pad2d(psf)))./...
    ((abs(F(pad2d(psf)))).^2+miu))));
reversed = @(x) flip(flip(x,1),2);
% my_xcorr2 = @(x,y) Ft(F(x).*conj(F(y)));
my_xcorr2_pad = @(x,y) Ft(F(pad2d(x)).*conj(F(pad2d(y))));
VEC = @(x) x(:);
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));

% general settings
main_folder = '0305exp/';
sub_folders = {'gfp/','tissue/','tissue_s/','target/'};
num_datasets = 4;
flag_bg = [1,1,1,0];
psf_folder = 'psf/';
psf_sub_folders = {'3ms/','5ms/','10ms/','15ms/','20ms/','30ms/','50ms/'};
num_measure_copies = 32;
num_bg_copies = 16;
num_exposures = 7;
exp_time = [3,5,10,15,20,30,50];
rows = 1944;
cols = 2592;
num_zscan = 21;
psf = zeros(rows,cols,num_zscan);

% load all PSFs and generate HDR PSFs
for i = 1:21
    psf_stack = cell(1,num_exposures);
    for j = 1:num_exposures
        tmp_filename = [main_folder,psf_folder,cell2mat(psf_sub_folders(j)),num2str(i,'%.2d'),'.tif'];
        psf_stack{j} = tmp_filename;
    end
    img_hdr = makehdr(psf_stack,'RelativeExposure',exp_time);
    psf(:,:,i) = img_hdr./sum(img_hdr(:));
end
clear psf_stack tmp_filename img_hdr;
psfs = zeros(rows/2, cols/2, num_zscan);
for i = 1:21
    tmp = psf(:,:,i);
    tmp = average_shrink(tmp,2);
    tmp = tmp./sum(tmp(:));
    psfs(:,:,i) = tmp;
end
clear tmp;

% load all datasets and subtract bg if available
y = zeros(rows,cols,num_datasets);
ys = zeros(rows/2,cols/2,num_datasets);
for i = 1:num_datasets
    if flag_bg(i) == 1
        tmp_dataset = zeros(rows,cols,num_measure_copies);
        tmp_bg = zeros(rows,cols,num_bg_copies);
        for j = 1:num_measure_copies
            tmp_img_file_name = [main_folder,cell2mat(sub_folders(i)),'y (',num2str(j),').tif'];
            tmp = im2double(imread(tmp_img_file_name));
            tmp_dataset(:,:,j) = tmp;
        end
        for j = 1:num_bg_copies
            tmp_bg_file_name = [main_folder,cell2mat(sub_folders(i)),'bg (',num2str(j),').tif'];
            tmp = im2double(imread(tmp_bg_file_name));
            tmp_bg(:,:,j) = tmp;
        end
        tmp_img = mean(tmp_dataset,3);
        tmp_bg = mean(tmp_bg,3);
        tmp_img = tmp_img - tmp_bg;
        tmp_img(tmp_img<=0) = 0;
        y(:,:,i) = tmp_img./max(tmp_img(:));
        tmp_img_s = average_shrink(tmp_img,2);
        ys(:,:,i) = tmp_img_s./max(tmp_img_s(:));
    else
        tmp_dataset = zeros(rows,cols,num_measure_copies);
        for j = 1:num_measure_copies
            tmp_img_file_name = [main_folder,cell2mat(sub_folders(i)),'y (',num2str(j),').tif'];
            tmp = im2double(imread(tmp_img_file_name));
            tmp_dataset(:,:,j) = tmp;
        end
        tmp_img = mean(tmp_dataset,3);
        tmp_img(tmp_img<=0) = 0;
        y(:,:,i) = tmp_img./max(tmp_img(:));
        tmp_img_s = average_shrink(tmp_img,2);
        ys(:,:,i) = tmp_img_s./max(tmp_img_s(:));
    end
    disp([num2str(i),'/',num2str(num_datasets),' datasets done.']);
end
clear tmp_dataset tmp_bg tmp_img_file_name tmp_bg_file_name tmp tmp_img tmp_bg tmp_img_s;

% try tikhonov to find out the best matched PSF
mu = 0.0001;
tik_results = zeros(rows/2, cols/2, num_datasets, num_zscan);
for i = 1:num_datasets
    for j = 1:num_zscan
        tik_results(:,:,i,j) = deconv_tik(ys(:,:,i), psfs(:,:,j), mu);
    end
end
idx = 1;
for i = 1:num_zscan
    figure;imagesc(clip(tik_results(:,:,idx,i),0,1000)), axis image off;truesize;colormap(jet(256));colorbar;
    title(['psf idx: ',num2str(i)]);
end


para = [];
para.mu1 = 1;
para.mu2 = 1;
para.mu3 = 1;
para.mu4 = 1;
para.color = 'gray';
para.maxiter = 30;
para.clip_min= 0;
para.clip_max= 0.1;
para.tv_thres = 0.001;
para.l1_thres = 0.01;
save_folder = 'tuning/';

mu2 = [0.1, 1, 10];
mu3 = [0.1, 1, 10];
mu4 = [0.1, 1, 10];
tv_thres = [0.001, 0.01, 0.1];
l1_thres = [0.001, 0.01, 0.1];
clip_max = [0.01, 0.1];

generate_name = @(a,b,c,d,e,f) [num2str(a),'_',num2str(b),'_',num2str(c),...
    '_',num2str(d),'_',num2str(e),'_',num2str(f),'.mat'];

for a = 3 %1:3
    for b=1:3
        for c=1:3
            for d=1:3
                for e=1:3
                    for f=2 %1:2
                        para.mu2 = mu2(a);
                        para.mu3 = mu3(b);
                        para.mu4 = mu4(c);
                        para.tv_thres = tv_thres(d);
                        para.l1_thres = l1_thres(e);
                        para.clip_max = clip_max(f);
                        [est,~] = LSI_model_ADMM_Solver(tmp_y_bs,tmp_psf2,para);
                        close all;
                        filename = generate_name(a,b,c,d,e,f);
                        save([save_folder,filename],'est');
                        disp(filename);
                    end
                end
            end
        end
    end
end

% normalize PSF according to target measurements
row_range = [1, 300, 700,972];
col_range = [1, 500, 900, 1296];
relative_intensity1 = [66.6, 119.5, 84.5; 85.5, 146.7, 108.6; 63.5, 113.9, 83.8];
relative_intensity1 = relative_intensity1./max(relative_intensity1(:));
relative_intensity2 = [110,123,114;119,135,120;95,107,95];
relative_intensity2 = relative_intensity2./max(relative_intensity2(:));


tmp_psf3 = tmp_psf2;
for i = 1:3
for j = 1:3
patch = tmp_psf3(row_range(i):row_range(i+1), col_range(j):col_range(j+1));
patch = patch./sum(patch(:));
tmp_psf3(row_range(i):row_range(i+1), col_range(j):col_range(j+1)) = patch.*relative_intensity2(i,j);
end
end


%% load two plane measurements gfp and subtract bg
name1 = 'gfp_-1.2 (';
name2 = 'gfp_-0.2 (';
y2 = zeros(1944,2592,2);
tmp = zeros(1944,2592,16);
for i = 1:16
    tmp(:,:,i) = im2double(imread([name1,num2str(i),').tif']));
end
img1 = mean(tmp,3);
for i = 1:16
    tmp(:,:,i) = im2double(imread([name2,num2str(i),').tif']));
end
img2= mean(tmp,3);
for i = 1:16
    tmp(:,:,i) = im2double(imread(['bg (',num2str(i),').tif']));
end
bg = mean(tmp,3);
y2(:,:,1) = imsubtract(img1,bg);
y2(:,:,2) = imsubtract(img2,bg);
y2bs = y2;
y2bs(:,:,1) = bg_removal(y2bs(:,:,1),128);
y2bs(:,:,2) = bg_removal(y2bs(:,:,2),128);





% 
% recons = zeros(2*rows, 2*cols, num_psfs);
% for i = 1:num_psfs
%     [tmp_out,~] = LSI_model_ADMM_Solver(ys_bs(:,:,idx),psfn(:,:,i),para);
%     recons(:,:,i) = tmp_out;
% end
% 
% 
% % process resolution target images
% rt_path = 'resolution target 0304/';
% num_imgs = 64;
% yrt_stack = zeros(rows, cols, num_imgs);
% for i = 1:num_imgs
%     name = [rt_path,'frt (',num2str(i),').tif'];
%     yrt_stack(:,:,i) = im2double(imread(name));
% end
% yrt = mean(yrt_stack, 3);
% yrts = average_shrink(yrt, 2);
% for i = 11:15
%     [est, history] = LSI_model_ADMM_Solver(yrts,psfs(:,:,i),para);
% end






    