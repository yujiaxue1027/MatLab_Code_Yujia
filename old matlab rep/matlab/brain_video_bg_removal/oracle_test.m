%% load data
load('oracle_location_data.mat');

%% define func
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
myxcorr2 = @(x,y) crop2d(Ft(F(pad2d(x)).*F(pad2d(reversed(y)))));

%% display gt, measure
figure,imagesc(x_test),axis image off;colormap gray;title('gt');
figure,imagesc(y_test_nf),axis image off;colormap gray; title('noisy measure');

%% no cropping model
para1 = [];
para1.color ='gray';
para1.maxiter = 10;
para1.mu = 1;

[est_nf_nc,est_nf_nc_hist] = ...
    admm_oracle_location(y_test_nf, psf_test, mask_test, para1, 1);
figure,imagesc(est_nf_nc),axis image off;colormap gray; title('noise free, no cropping');
  
[est_n_nc,est_n_nc_hist] = ...
    admm_oracle_location(y_test, psf_test, mask_test, para1, 1);
figure,imagesc(est_n_nc),axis image off;colormap gray; title('noisy, no cropping');

%% cropping
para2 = [];
para2.color ='gray';
para2.maxiter = 10;
para2.mu1 = 1;
para2.mu2 = 0.5;

[est_nf_c,est_nf_c_hist] = ...
    admm_oracle_location_crop(y_test_nf, psf_test, mask_test, para2, 1);
figure,imagesc(crop2d(est_nf_c)),axis image off;colormap gray; title('noise free, cropping');
  
[est_n_c,est_n_c_hist] = ...
    admm_oracle_location_crop(y_test, psf_test, mask_test, para2, 1);
figure,imagesc(crop2d(est_n_c)),axis image off;colormap gray; title('noisy, cropping');