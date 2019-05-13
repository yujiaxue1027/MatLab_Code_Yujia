function deconvolution_tikhonov(mu)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
my_deconv_tik = @(y,psf,miu) crop2d(real(Ft(F(pad2d(y)).*conj(F(pad2d(psf)))./((abs(F(pad2d(psf)))).^2+miu))));
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));
load('data.mat','ys_bs','npsfs');
y = ys_bs(:,:,1);
psf = npsfs(:,:,1);
[yrow,ycol] = size(y);
[psfrow,psfcol] = size(psf);
disp([num2str(yrow),num2str(ycol),num2str(psfrow),num2str(psfcol),num2str(mu)]);
xhat = my_deconv_tik(y,psf,mu);
xhat = linear_normalize(xhat);
imwrite(uint8(255*xhat),['save/img_mu_',num2str(mu),'.png']);
save(['save/xhat_mu_',num2str(mu),'.mat'],'xhat');
end