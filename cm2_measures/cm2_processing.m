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
my_xcorr2 = @(x,y) Ft(F(x).*conj(F(y)));
VEC = @(x) x(:);

% load images and preprocess background
filenames = {'gfp_array_cross3.tif','gfp_array_cross2.tif','gfp_array_cross.tif','gfp_cross2.tif'...
    ,'gfp_cross.tif','array_gel_exem.tif'};
num_imgs = 6;
rows = 1944;
cols = 2592;
ys = zeros(rows,cols,num_imgs);
bgs = zeros(rows,cols,num_imgs);
ys_bs = zeros(rows,cols,num_imgs);

kernel_size = 250;
for i = 1:num_imgs
    tmp = im2double(imread(cell2mat(filenames(i))));
    ys(:,:,i) = tmp;
    bgs(:,:,i) = imopen(tmp, strel('disk', kernel_size));
    ys_bs(:,:,i) = imsubtract(tmp, bgs(:,:,i));
end
for i = 1:num_imgs
    figure;
    subplot(1,3,1);
    imagesc(ys(:,:,i)),axis image off;colormap parula;
    title('original measure');
    subplot(1,3,2);
    imagesc(bgs(:,:,i)),axis image off;colormap parula;
    title('estimated bg');
    subplot(1,3,3);
    imagesc(ys_bs(:,:,i)),axis image off;colormap parula;
    title('bg subtracted measure');
end
for i = 1:num_imgs
    figure;
    subplot(1,2,1);
    imagesc(my_xcorr2(ys(:,:,i),ys(:,:,i)));axis image off;
    title('xcorr2 of raw measure');
    subplot(1,2,2);
    imagesc(my_xcorr2(ys_bs(:,:,i),ys_bs(:,:,i)));axis image off;
    title('xcorr2 of bg subtracted measure');
end
    
% idx = 4;
% load PSFs from low light and high light
num_psfs = 31;
path_low = 'psf_calibration0228/low/';
path_high = 'psf_calibration0228/high/';
psf_low = zeros(rows,cols,num_psfs);
psf_high = zeros(rows,cols,num_psfs);
for i = 1:num_psfs
    name = [num2str(i,'%.2d'),'.tif'];
    psf_low(:,:,i) = flip(im2double(imread([path_low, name])),1);
    psf_high(:,:,i) = flip(im2double(imread([path_high, name])),1);
end
for i = 1:num_psfs
    figure;
    subplot(1,2,1);
    imagesc(psf_low(:,:,i)),axis image off;colormap gray;title([num2str(i),' low']);
    subplot(1,2,2);
    imagesc(psf_high(:,:,i)),axis image off;colormap gray;title([num2str(i),' high']);
end

% merge high/low light psf together
psf = psf_high;
psf(600:end,1:1800,:) = psf_low(600:end,1:1800,:);
for i = 1:num_psfs
    figure;
    imagesc(psf(:,:,i)),axis image off;colormap gray;truesize;
end

% normalize energy in each sub-psf
radius = 256;
psfn = zeros(rows+2*radius, cols+2*radius, num_psfs);
psf_padded = padarray(psf,[radius,radius,0]);
centers = [80,454; 73, 1322; 948, 472; 932, 1341; 905, 2229; 1814, 486; 1801, 1360; 1785, 2250];
energy_threshold = 10;
for i = 1:num_psfs
%     figure;
%     imagesc(psf_padded(:,:,i)),axis image off,colormap gray;truesize;
%     hold on;
    for j = 1:8
        patch = psf_padded(centers(j,1):centers(j,1)+2*radius, centers(j,2):centers(j,2)+2*radius, i);
        if sum(VEC(patch)) >= energy_threshold
            disp(['i=',num2str(i),'/31, j=',num2str(j),'/31, sum = ',num2str(sum(VEC(patch)))]);
            [row_co,col_co] = locate_mass_center(patch);
%             plot(centers(j,2)+col_co, centers(j,1)+row_co, 'y+');
            pos = [centers(j,2)+col_co-radius,centers(j,1)+row_co-radius,2*radius+1,2*radius+1];
%             rectangle('Position',pos,'edgecolor','r');
            roi = psf_padded(centers(j,1)+row_co-radius:centers(j,1)+row_co+radius, ...
                centers(j,2)+col_co-radius:centers(j,2)+col_co+radius, i);
            roi = roi./sum(roi(:));
            psfn(centers(j,1)+row_co-radius:centers(j,1)+row_co+radius, ...
                centers(j,2)+col_co-radius:centers(j,2)+col_co+radius, i) = roi;
        end
    end
end
for i = 1:num_psfs
    psfn(:,:,i) = psfn(:,:,i)./sum(VEC(psfn(:,:,i)));
end
psfn = psfn(radius+1:radius+rows, radius+1:radius+cols,:);   
for i = 1:num_psfs
figure;
imagesc(psfn(:,:,i));axis image off;colormap gray;truesize;title(i);
end

% start deconvolution
idx = 4;
% mu = 0.001;
% for i = 1:num_psfs
%     tmp = deconv_tik(ys_bs(:,:,idx),psfn(:,:,i),mu);
%     figure;
%     imagesc(clip(tmp,0,1000)),axis image off;colormap gray;truesize;title(i);
% end
para = [];
para.mu1 = 10;
para.mu2 = 0.2;
para.mu3 = 1;
para.mu4 = 0.2;
para.color = 'gray';
para.maxiter = 30;
para.clip_min=0;
para.clip_max=10;
para.tv_thres = 0.02;
para.l1_thres = 0.1;

recons = zeros(2*rows, 2*cols, num_psfs);
for i = 1:num_psfs
    [tmp_out,~] = LSI_model_ADMM_Solver(ys_bs(:,:,idx),psfn(:,:,i),para);
    recons(:,:,i) = tmp_out;
end


% process resolution target images
rt_path = 'resolution target 0304/';
num_imgs = 64;
yrt_stack = zeros(rows, cols, num_imgs);
for i = 1:num_imgs
    name = [rt_path,'frt (',num2str(i),').tif'];
    yrt_stack(:,:,i) = im2double(imread(name));
end
yrt = mean(yrt_stack, 3);
yrts = average_shrink(yrt, 2);
for i = 11:15
    [est, history] = LSI_model_ADMM_Solver(yrts,psfs(:,:,i),para);
end






    