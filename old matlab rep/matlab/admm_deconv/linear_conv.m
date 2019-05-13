function [output] = linear_conv(img,kernel)
    F = @(x) fftshift(fft2(ifftshift(x)));
    Ft = @(x) fftshift(ifft2(ifftshift(x)));
    output = Ft(F(padarray(img,size(img))).*F(padarray(kernel,size(kernel))));
    output = output(size(output,1)/3+1:2*size(output,1)/3,...
        size(output,2)/3+1:2*size(output,2)/3);
end