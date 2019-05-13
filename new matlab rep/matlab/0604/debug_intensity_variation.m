%% debug intensity variation
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

obj = im2double(rgb2gray(imread('target.jpg')));
psf0 = im2double(imread('old_psf_6.tif'));
obj = make_the_same(obj,psf0);
psf0(psf0<=0.08) = 0;
psf0 = medfilt2(psf0, [3,3]);
psf0 = psf0./sum(psf0(:));

y0 = conv2d(obj, psf0);

psf1 = psf0;
ratio = 1.8;
psf1(670:1300,:) = psf1(670:1300,:).*ratio;
psf1(:,990:1600) = psf1(:,990:1600).*ratio;
psf1 = psf1./sum(psf1(:));

deconv_0 = deconv_tik(y0, psf0, 0.001);
deconv_1 = deconv_tik(y0, psf1, 0.001);

figure,imagesc(deconv_0),axis image off;colormap gray
figure,imagesc(deconv_1),axis image off;colormap gray
