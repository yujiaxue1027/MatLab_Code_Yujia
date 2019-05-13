function out = initCWD(data, PAR)
%% init for X
sz = [size(data.g,1),size(data.g,2)];  % size of the input image
col = size(data.g,3);  % number of color channels
k = PAR.nBase;

fprintf('Generating PSF coefficients ...');
t = tic;
h = reshape(data.h(:,:,1:k),[],k);
if isfield(data,'a')  % can be used  to force specific coefficients
    M = reshape(data.a',[],k);
else  % that are otherwise computed using least-square fit
    M = reshape(((h'*h)\(h'*data.D))',[],k);
end
fprintf(' done (%.3fs)\n', toc(t));

out.sz = sz;
out.col = col;
out.k = k;
out.M = M;
out.MtM = PAR.gamma/PAR.beta * sum(M.^2,2);
out.FICtC = fft2(circshift(padarray([0,-1,0;-1,4+PAR.gamma/PAR.beta,-1;0,-1,0],sz-3,0,'post'),[-1,-1]));
out.x = zeros(prod(sz),k,col);
out.a = out.x;

%% init for U
FH = zeros([sz, k]); % FFT of PSFs
FIHHt = zeros(sz);   % sum_i conj(FH_i)*FH_i
for i = 1:k
    FH(:,:,i) = fft2(circshift(padarray(data.h(:,:,i),sz-data.sPsf,0,'post'),-(data.sPsf-1)/2));
    FIHHt = FIHHt + abs(FH(:,:,i)).^2;
end
FIHHt = FIHHt + PAR.gamma/PAR.mu;

Htg = zeros([sz,k,col]);
for i = 1:k
    Htg(:,:,i,:) = real(ifft2(fft2(data.g).*repmat(conj(FH(:,:,i)),[1,1,col])));
end
Htg = PAR.mu/PAR.gamma * reshape(Htg,[],k,col);

out.tInit = toc(t(1));
out.FH = FH;
out.FIHHt = FIHHt; % repmat(FHtH,[1,1,sg(3)]);
out.Htg = Htg;
out.u = zeros(size(data.g));
out.b = zeros([sz,col,2]);
out.v = out.b;
out.eu = [];
out.t = [];

end
