function output = BT(input,Hs_conj)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
output = zeros(size(Hs_conj));
for i = 1:size(Hs_conj,3)
    output(:,:,i) = real(Ft(F(input).*Hs_conj(:,:,i)));
end

end

