function output = inverse_cac(input,psf,length,reg)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
psf = padarray(psf,size(psf));
input = padarray(input,size(input));
f_mea = F(input);
f_psf = F(psf);
para = reg*max(abs(f_psf(:)));
f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
result = real(Ft(f_res));
result = mytruncate(result,-1,1);
output = result(round(size(result,1)/2)-length/2:round(size(result,1)/2)+length/2-1,...
    round(size(result,1)/2)-length/2:round(size(result,1)/2)+length/2-1);
end

