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


% psf_w = im2double(imread('psf_white.tif'));
% psf_w = psf_w./sum(psf_w(:));
% 
% y1 = im2double(imread('img1_negative_resolution_non_fluorescence.tif'));
% recon1 = LSI_model_ADMM_Solver(y1, psf_w, para_lsi);

% psf_fitc = im2double(imread('psf_fitc.tif'));
% psf_fitc = psf_fitc./sum(psf_fitc(:));
% psf_fitc_s = average_shrink(psf_fitc, 2);

y1 = im2double(imread('img2_fluorescence.tif'));
% y2 = im2double(imread('tissue1.tif'));
% y3 = im2double(imread('tissue2.tif'));
y1s = average_shrink(y1,2);
% y2s = average_shrink(y2,2);
% y3s = average_shrink(y3,2);

psfs = zeros(1944,2592,22);
for i = 1:22
    name = ['psfs/',num2str(i+51),'.tif'];
    tmp = im2double(imread(name));
    psfs(:,:,i) = tmp./sum(tmp(:));
end

tik_tau = 0.001;
recon_tik = zeros(1944,2592,22);
for i = 1:22
    tmp = deconv_tik(y1,psfs(:,:,i), 0.001);
    recon_tik(:,:,i) = clip(tmp,0,1000);
    close all;
    figure;
    imagesc(recon_tik(:,:,i));colormap gray;axis image;colorbar;
    disp(i)
end


recon1 = LSI_model_ADMM_Solver(y1, psf_fitc, para_lsi);
% recon2 = LSI_model_ADMM_Solver(y2, psf_fitc, para_lsi);
% recon3 = LSI_model_ADMM_Solver(y3, psf_fitc, para_lsi);

recon1s = LSI_model_ADMM_Solver(y1s, psf_fitc_s, para_lsi);
% recon2s = LSI_model_ADMM_Solver(y2s, psf_fitc_s, para_lsi);
% recon3s = LSI_model_ADMM_Solver(y3s, psf_fitc_s, para_lsi);