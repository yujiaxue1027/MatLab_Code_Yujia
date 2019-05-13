function output = Ft3D(input)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
output = zeros(size(input));
for i = 1:size(input,3)
    output(:,:,i) = Ft(input(:,:,i));
end

end

