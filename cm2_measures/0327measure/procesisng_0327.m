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

ys = zeros(1944,2592,3);
for i = 1:3
    ys(:,:,i) = im2double(imread(['y',num2str(i),'.tif']));
end
psfs = zeros(1944,2592,56);
for i = 1:56
    tmp = double(imread(['psf/',num2str(i,'%.2d'),'.tif']));
% tmp = medfilt2(tmp,[3,3]);
% psfs(:,:,i) = tmp./sum(tmp(:)).*9;
    psfs(:,:,i) = tmp;
end
psfs(psfs<=10) = 0;
for i = 1:56
    psfs(:,:,i) = psfs(:,:,i)./sum(VEC(psfs(:,:,i))).*9;
end

% psfs(psfs<=0.001) = 0;


% for i = 6:16
% imwrite(crop2d(uint8(255*linear_normalize(deconv_tik(ys(:,:,3),psfs(:,:,i),0.1)))),['z_',num2str(i),'.png']);
% end

ys = cat(3,ys,ys);
for i = 4:6
    ys(:,:,i) = bg_removal(ys(:,:,i), 128);
end
realpsf = zeros(size(ys));
realpsf(:,:,1) = psfs(:,:,23);
realpsf(:,:,2) = psfs(:,:,20);
realpsf(:,:,3) = psfs(:,:,11);
realpsf(:,:,4) = psfs(:,:,23);
realpsf(:,:,5) = psfs(:,:,20);
realpsf(:,:,6) = psfs(:,:,11);

psfs = realpsf;





