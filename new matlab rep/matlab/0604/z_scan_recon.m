%% z scan reconstruction
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

num_scans = 61; %61
src_folder = 'zscan0720/';
load resolution_target_result_1.mat
clear recon_lsv_large recon_noise recon_tik y yn
% para2.maxiter = 1; %
num_iter = para2.maxiter;
rows = 1944;
cols = 2592;
all_history = zeros(num_scans,rows,cols,num_iter);
all_recons = zeros(num_scans,rows,cols);
all_recons_tik = zeros(num_scans,rows,cols);
psf_tik = zeros(rows,cols);
for i = 1:size(bpsf_bank,3)
    psf_tik = psf_tik + bpsf_bank(:,:,i)*intp_mask(rows/2,cols/2,i);
end
psf_tik = psf_tik./sum(psf_tik(:));

for i = 1:num_scans
    input_img = im2double(imread([src_folder,num2str(i,'%.2d'),'.tif' ]));
    all_recons_tik(i,:,:) = deconv_tik(input_img,psf_tik,0.0001);
    [tmp_out, tmp_hist] = admm_lsv2(input_img,bpsf_bank,intp_mask,para2);
    all_recons(i,:,:) = tmp_out;
    all_history(i,:,:,:) = tmp_hist;
    close all;
    disp(i);
end

load resolution_target_result_4.mat
clear recon_sparse_n recon_sparse_nf y yn
% para2.maxiter = 1; %
num_iter = para2.maxiter;
rows = 1944;
cols = 2592;
% all_history_sparse = zeros(num_scans,rows,cols,num_iter);
all_recons_sparse = zeros(num_scans,rows,cols);

for i = 1:num_scans
    input_img = im2double(imread([src_folder,num2str(i,'%.2d'),'.tif' ]));
    [tmp_out, tmp_hist] = admm_lsv2(input_img,bpsf_bank_sparse,intp_mask_sparse,para2);
    all_recons_sparse(i,:,:) = tmp_out;
%     all_history_sparse(i,:,:,:) = tmp_hist;
    close all;
    disp(i);
end


