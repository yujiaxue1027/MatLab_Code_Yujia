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
% my_xcorr2 = @(x,y) Ft(F(x).*conj(F(y)));
my_xcorr2_pad = @(x,y) Ft(F(pad2d(x)).*conj(F(pad2d(y))));
VEC = @(x) x(:);
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));

y = im2double(imread('xpolar_led_gfp2.tif'));
num_psf = 61;
rows = 1944;
cols = 2592;
psf_stack = zeros(rows,cols,num_psf);
for i = 1:61
    tmp = im2double(imread(['psf0322/',num2str(i,'%.2d'),'.tif']));
    psf_stack(:,:,i) = tmp./sum(tmp(:)).*9;
    psf_stack(:,:,i) = medfilt2(psf_stack(:,:,i),[3,3]);
end

for i = 10:1:20
    tmp = deconv_tik(y,psf_stack(:,:,i),0.001);
    figure,imagesc(clip(crop2d(tmp),0,100));axis image off;truesize,colormap(gray(256));
    title(i);
end
