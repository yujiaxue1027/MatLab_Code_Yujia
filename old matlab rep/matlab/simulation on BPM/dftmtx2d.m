function output = dftmtx2d( M,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = zeros(M*N,M*N);
for p = 0:M*N-1
    for q = 0:M*N-1
        k = mod(p,M);
        l = (p-k)/M;
        m = mod(q,M);
        n = (q-m)/M;
        output(p+1,q+1) = exp(-1*1i*2*pi*(m*k/M+n*l/N));
    end
end
% my2ddft = @(x) fftshift(reshape(output*reshape(ifftshift(x),[M*N,1]),M,N));
% then F(x) = my2dfft(x) where F is fftshift(fft2(ifftshift()))

end

