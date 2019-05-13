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

load resolution_targets_dataset.mat

num_scans = 41;
src_folder = 'tissue_diffuser/';

measures = zeros(1944, 2592, num_scans);
measures_s = zeros(972, 1296, num_scans);

for i = 1:num_scans
    img = im2double(imread([src_folder, num2str(i,'%.2d'), '.tif']));
    measures(:,:,i) = img;
    measures_s(:,:,i) = average_shrink(img,2);
    disp(i);
end

tik_s = zeros(972, 1296, num_scans);
for i = 1:num_scans
    tmp = measures_s(:,:,i);
    tmp = deconv_tik(tmp,psf_on_axis_s, 0.001);
    tik_s(:,:,i) = tmp;
    tmp = clip(tmp,0,100);
    tmp = tmp./max(tmp(:));
    imwrite(uint8(250*tmp), [src_folder,'tik_s_',num2str(i),'.png']);
    disp(i);
end

%%
% measure_s_2 = measures_s;
candidates_2 = 6:10;
recon_lsi_2 = zeros(972, 1296, length(candidates_2));
% recon_lsv = zeros(972, 1296, length(candidates));

for i = 1:length(candidates_2)
    tmp = measures_s_2(:,:,candidates_2(i));
    [tmplsi,~] = LSI_model_ADMM_Solver(tmp, psf_on_axis_s, para_lsi);
%     [tmplsv,~] = LSV_model_ADMM_Solver(tmp, bpsf_s, intp_mask_s, para_lsv);
    recon_lsi_2(:,:,i) = crop2d(tmplsi);
%     recon_lsv(:,:,i) = tmplsv;
    close all;
    disp(i);
end
save yfp2.mat recon_lsi_2 measures_s_2

% measure_s_3 = measures_s;
candidates_3 = 1:2:31;
recon_lsi_3 = zeros(972, 1296, length(candidates_3));
% recon_lsv = zeros(972, 1296, length(candidates));

for i = 1:length(candidates_3)
    tmp = measures_s_3(:,:,candidates_3(i));
    [tmplsi,~] = LSI_model_ADMM_Solver(tmp, psf_on_axis_s, para_lsi);
%     [tmplsv,~] = LSV_model_ADMM_Solver(tmp, bpsf_s, intp_mask_s, para_lsv);
    recon_lsi_3(:,:,i) = crop2d(tmplsi);
%     recon_lsv(:,:,i) = tmplsv;
    close all;
    disp(i);
end
save yfp3.mat recon_lsi_3 measures_s_3


% measure_s_4 = measures_s;
candidates_4 = 1:2:31;
recon_lsi_4 = zeros(972, 1296, length(candidates_4));
% recon_lsv = zeros(972, 1296, length(candidates));

for i = 1:length(candidates_4)
    tmp = measures_s_4(:,:,candidates_4(i));
    [tmplsi,~] = LSI_model_ADMM_Solver(tmp, psf_on_axis_s, para_lsi);
%     [tmplsv,~] = LSV_model_ADMM_Solver(tmp, bpsf_s, intp_mask_s, para_lsv);
    recon_lsi_4(:,:,i) = crop2d(tmplsi);
%     recon_lsv(:,:,i) = tmplsv;
    close all;
    disp(i);
end
save yfp4.mat recon_lsi_4 measures_s_4


