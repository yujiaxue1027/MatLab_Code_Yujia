%% try on cm2 psfs
rows = 256;
cols = 256;
layers = 8;

F3D = @(x) fftshift(fftn(ifftshift(x)));
Ft3D = @(x) fftshift(ifftn(ifftshift(x)));
VEC = @(x) x(:);
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
F2D = @(x) fftshift(fft2(ifftshift(x)));
Ft2D = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft2D(F2D(pad2d(obj)).*F2D(pad2d(psf)))));
deconv_tik = @(y,psf,miu) crop2d(real(Ft2D(F2D(pad2d(y)).*...
    conj(F2D(pad2d(psf)))./...
    ((abs(F2D(pad2d(psf)))).^2+miu))));
reversed = @(x) flip(flip(x,1),2);
my_xcorr2_pad = @(x,y) Ft2D(F2D(pad2d(x)).*conj(F2D(pad2d(y))));
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));

%% load image and synthesize measure
obj = im2double(rgb2gray(imread('target_example.jpg')));
obj = imresize(obj,[256,256]);
obj = obj - 0.08;
obj(obj<=0) = 0;
obj = obj./max(VEC(obj));
tmp = zeros(rows,cols,layers);
cols_per_layer = cols/layers;
for i = 1:layers
    tmp(:,(i-1)*cols_per_layer+1:i*cols_per_layer,i) = obj(:,(i-1)*cols_per_layer+1:i*cols_per_layer);
end
obj = tmp;

%% load psf and synthesize measurement
load psf_data_for_simulation.mat
% psfs = psfs_256(:,:,[18, 21, 24, 27, 30, 33, 36, 39]);
psfs = psfs_256(:,:,[10, 15, 20, 25, 30, 35, 40, 45]);
psfs_reverse = flip(psfs,3);
y = conv3d(obj,psfs_reverse);

%% deconv 3d
load default_parameters.mat
para.display_flag = 1;
para.tau_l1 = 1e-1; %3e-2
para.tau_tv = 2e-2;
% para.clip_max = 1;
para.maxiter = 512;  %1024
% [recon,~] = ADMM_LSI_deconv_3D(y_easy,psf_reverse,para);
[recon,~] = ADMM_LSI_deconv_3D(y,psfs_reverse,para);

% try single layer 2D deconv ADMM
para_single = para;
para.maxiter = 512;
para.tau_l1 = 0.05;
para.tau_tv = 1e-4;
[recon_single,~] = ADMM_LSI_deconv_l1tv(y,psfs_reverse(:,:,7),para);

