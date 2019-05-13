%% admm deconv with rescaled psf
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
load default_para.mat
load mbrain_1.mat
load('mousebrain.mat', 'img0');
psf0=psf;

new_dim = [680:10:750];
scale = [];
recons = zeros(1500,1500,length(new_dim));
for i = 1:length(new_dim)
scale = [scale,new_dim(i)/750];
psf2 = imresize(psf0,[new_dim(i),new_dim(i)]);
psf = padarray(psf2,0.5*(size(psf0)-size(psf2)));
save rescaled.mat psf y
pause(5);
[recon,~,~] = admm_recon3('rescaled.mat',para);
recons(:,:,i) = recon;
end