%% script to process 3D beads measurements
%% predefined functions
F3D = @(x) fftshift(fftn(ifftshift(x)));
Ft3D = @(x) fftshift(ifftn(ifftshift(x)));
VEC = @(x) x(:);
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
F2D = @(x) fftshift(fft2(ifftshift(x)));
Ft2D = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft2D(F2D(pad2d(obj)).*F2D(pad2d(psf)))));
deconv_tik = @(y,psf,miu) crop2d(real(Ft2D(F2D(pad2d(y)).*conj(F2D(pad2d(psf)))./((abs(F2D(pad2d(psf)))).^2+miu))));
reversed = @(x) flip(flip(x,1),2);
my_xcorr2_pad = @(x,y) Ft2D(F2D(pad2d(x)).*conj(F2D(pad2d(y))));
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));

%% read in images
rows = 1944;
cols = 2592;
y_0 = im2double(imread('15_25_T2.tif'));

y = bg_removal(y_0, 256); 
y = linear_normalize(y);

%% load in all PSFs
num_psf = 701;
psf_stack = zeros(rows,cols,num_psf);
for i = 1:num_psf
    tmp = im2double(imread(['psfs_p2/',num2str(i,'%.2d'),'.tif']));
    tmp = medfilt2(tmp,[3,3]);
    psf_stack(:,:,i) = tmp./sum(tmp(:)).*9;
    disp(i);
end

%% temp code to find which psf to use
for i = 1:10:701
    imwrite(crop2d(uint8(255*linear_normalize(deconv_tik(y,psf_stack(:,:,i),0.1)))),['temp/y_',num2str(i),'.png']);
end



%% export data to process on SCC
psfs1 = psf_stack(:,:,280:15:505);
psfs2 = psf_stack(:,:,200:20:500);
load default_parameters.mat
para.display_flag = 1;
para.tau_l1 = 1e-1; %3e-2
para.tau_tv = 2e-2;
% para.clip_max = 1;
para.maxiter = 512;  %1024
recon1 = ADMM_LSI_deconv_3D(y1,psfs1,para);

psfs = psf_stack(:,:,135:5:290);

% psfs = psf_stack(:,:,73:2:328);

load default_parameters.mat
para.display_flag = 1;
para.tau_l1 = 1e-1; %3e-2
para.tau_tv = 2e-2;
% para.clip_max = 1;
para.maxiter = 512;  %1024
recon1 = ADMM_LSI_deconv_3D(y1,psfs1,para);








