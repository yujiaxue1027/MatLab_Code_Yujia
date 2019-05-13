function out = initRWD(data, PAR)
%% init for X
sz = [size(data.g2,1),size(data.g2,2)];  % size of the input image
col = size(data.g2,3);  % number of color channels
k = PAR.nBase;

fprintf('Generating PSF coefficients ...');
t = tic;
h = reshape(data.h2(:,:,1:k),[],k);
if isfield(data,'a')  % can be used  to force specific coefficients
    M = reshape(data.a',[],k);
else  % that are otherwise computed using least-square fit
    M = reshape(((h'*h)\(h'*data.D2))',[],k);
end
fprintf(' done (%.3fs)\n', toc(t));

fprintf('Computing M''M and M''g ...');
t(2) = tic;
MtM = repmat(reshape(M,[sz,k]),[1,1,1,k]);
MtM = MtM .* permute(MtM,[1,2,4,3]);
Mtg = repmat(reshape(M,[sz,k]),[1,1,1,col]) .* repmat(reshape(data.g2,[sz,1,col]),[1,1,k,1]);
% pad MtM and Mtg by zeros (for eliminating border artefacts)
r = (data.sPsf+1)/2 - 1;
Mtg = reshape(padarray(Mtg,r),[],k,col);
MtM = reshape(padarray(MtM,r),[],k,k);
fprintf(' done (%.3fs)\n', toc(t(2)));
su = sz + data.sPsf - 1;  % size of the resulting image
% mu/gamma*MtM + I -> MtM
MtM = PAR.mu/PAR.gamma * MtM;
for i = 1:k
    MtM(:,i,i) = MtM(:,i,i) + 1;
end

fprintf('Inverting M''M+I ...');
t(3) = tic;
MtM = cat(3, MtM, repmat(shiftdim(eye(k),-1),[size(MtM,1),1,1]));
for i = 1:k
    MtM(:,i,:) = MtM(:,i,:) ./ repmat(MtM(:,i,i),[1,1,2*k]);
    for j = i+1:k
        MtM(:,j,:) = MtM(:,j,:) - MtM(:,i,:).*repmat(MtM(:,j,i),[1,1,2*k]);
    end
end
for i = k:-1:1
    for j = 1:i-1
        MtM(:,j,:) = MtM(:,j,:) - MtM(:,i,:).*repmat(MtM(:,j,i),[1,1,2*k]);
    end
end
MtM = reshape(MtM(:,:,k+1:2*k),[],k);
fprintf(' done (%.3fs)\n', toc(t(3)));

Mtg = PAR.mu/PAR.gamma * Mtg; %reshape(,[],sg(3));

out.col = col;
out.su = su;
out.k = k;
out.iMtM = MtM;
out.Mtg = Mtg;
out.x = zeros(prod(su),k,col);
out.a = out.x;

%% init for U
FH = zeros([su, k]); % FFT of PSFs
FHtH = zeros(su);    % sum_i conj(FH_i)*FH_i
for i = 1:k
    FH(:,:,i) = fft2(data.h2(:,:,i),su(1),su(2));
    FHtH = FHtH + abs(FH(:,:,i)).^2;
end

% FD ... FFT of individual derivative operators (x, y, etc.)
FD = zeros([su,2]);
FD(:,:,1) = fft2([1 -1],su(1),su(2));
FD(:,:,2) = fft2([1;-1],su(1),su(2));
% FD(:,:,3) = fft2([1 0;0 -1] / sqrt(2),su(1),su(2));
% FD(:,:,4) = fft2([0 1;-1 0] / sqrt(2),su(1),su(2));
nDir = size(FD,3);
FDtD = repmat(sum(conj(FD).*FD,3),[1,1,col]);
FD = repmat(reshape(FD,[su,1,nDir]),[1,1,col,1]);

out.tInit = toc(t(1));
out.FH = FH;
out.FHtH = repmat(FHtH,[1,1,col]);
out.FD = FD;
out.FDtD = FDtD;
out.u = zeros([su,col]);
out.Hu = zeros([prod(su),k,col]);
out.b = zeros([su,col,nDir]);
out.v = out.b;
out.eu = [];
out.t = [];

end
