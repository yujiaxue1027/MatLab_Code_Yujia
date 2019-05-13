%% write 2d dft of a 2d mat into a single mat * vec multiplication
%% determine the H mat in bpm
M = 64;
N = 64;
% x = reshape(1:M*N,M,N);
x = phantom('Modified Shepp-Logan',M);
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

% y0 = fft2(x);

% fm = zeros(M*N,M*N);
% for p = 0:M*N-1
%     for q = 0:M*N-1
%         k = mod(p,M);
%         l = (p-k)/M;
%         m = mod(q,M);
%         n = (q-m)/M;
%         fm(p+1,q+1) = exp(-1*i*2*pi*(m*k/M+n*l/N));
%     end
% end
y0 = F(x);
fm = dftmtx2d(M,N);
my2ddft = @(x) fftshift(reshape(fm*reshape(ifftshift(x),[M*N,1]),M,N));
myy = my2ddft(x);
ifm = conj(fm)/(M*N);
my2didft = @(x) fftshift(reshape(ifm*reshape(ifftshift(x),[M*N,1]),M,N));
myxhat = my2didft(myy);
xhat = Ft(y0);



xhat1 = abs(ifft2(y0));
xhat2 = abs(ifft2(myy));
fm = dftmtx(64*64);
figure,imagesc(x),colorbar;
figure,imagesc(xhat1),colorbar;
figure,imagesc(xhat2),colorbar;
figure,imagesc(xhat1-xhat2),colorbar;