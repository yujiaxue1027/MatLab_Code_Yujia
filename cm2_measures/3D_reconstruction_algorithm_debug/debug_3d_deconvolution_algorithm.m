%% some debug code for 3d reconstruction algorithm derivation
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


%% generate some examples psfs
psfs = zeros(rows,cols,layers);
[x,y] = meshgrid(-128:127);
for i = 1:layers
    tmp = psfs(:,:,i);
    tmp(sqrt(x.^2 + y.^2)<=2*i) = 1;
    tmp = tmp./sum(VEC(tmp));
    psfs(:,:,i) = tmp;
end

psf2decompose = reshape(psfs,[rows*cols, layers]);
[U,S,V] = svds(psf2decompose,8);
psf_basis = reshape(U,[rows,cols,layers]);
psf_basis = abs(psf_basis);
for i = 1:layers
    tmp = psf_basis(:,:,i);
    tmp = tmp./sum(VEC(tmp));
    psf_basis(:,:,i) = tmp;
end
psfs = psf_basis;

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

%% compare slice by slice conv and sum AND 3d conv and slice
load caustics.mat
psf = psf(:,25+1:25+270,1:6:43);
npsf = zeros(256,256,8);
for i = 1:8
    tmp = psf(:,:,i);
    tmp = imresize(tmp,[256,256]);
    tmp = tmp./sum(tmp);
    npsf(:,:,i) = tmp;
end
psf = npsf;
psf_reverse = flip(psf,3);
y = conv3d(obj,psf_reverse);

clear i U S V x 

% load caustics.mat
% psf = psf(:,25+1:25+270,1:10:31);
% npsf = zeros(256,256,4);
% for i = 1:4
%     tmp = psf(:,:,i);
%     tmp = imresize(tmp,[256,256]);
%     tmp = tmp./sum(tmp);
%     npsf(:,:,i) = tmp;
% end
% psf = npsf;
% psf_reverse = flip(psf,3);
% obj_easy = zeros(256,256,4);
% obj_easy(:,:,1) = sum(obj(:,:,1:2),3);
% obj_easy(:,:,2) = sum(obj(:,:,3:4),3);
% obj_easy(:,:,3) = sum(obj(:,:,5:6),3);
% obj_easy(:,:,4) = sum(obj(:,:,7:8),3);
% y_easy = conv3d(obj_easy,psf_reverse);
% 
% clear i U S V x 

%% 
% psf_reverse = flip(psfs,3);
% y = conv3d(obj,psf_reverse);

%% easy case
% psf_reverse = flip(psfs(:,:,[1,3,5,7]),3);
% obj_easy = zeros(256,256,4);
% obj_easy(:,:,1) = sum(obj(:,:,1:2),3);
% obj_easy(:,:,2) = sum(obj(:,:,3:4),3);
% obj_easy(:,:,3) = sum(obj(:,:,5:6),3);
% obj_easy(:,:,4) = sum(obj(:,:,7:8),3);
% y_easy = conv3d(obj_easy,psf_reverse);

%% try 3d admm deconv
load default_parameters.mat
para.display_flag = 1;
para.tau_l1 = 3e-2;
para.tau_tv = 2e-2;
% para.clip_max = 1;
para.maxiter = 1024;
% [recon,~] = ADMM_LSI_deconv_3D(y_easy,psf_reverse,para);
[recon,~] = ADMM_LSI_deconv_3D(y,psf_reverse,para);

















