function data = iterU_RWD(data, PAR)

col = data.col;

%% solve for u
% H'(x+a) -> FHtx (in fourier domain)
FHtx = zeros([data.su,col]);
for i = 1:data.k
    img = reshape(data.x(:,i,:)+data.a(:,i,:),[data.su,col]);
    FHtx = FHtx + repmat(conj(data.FH(:,:,i)),[1,1,col]) .* fft2(img);%edgetaper(img,data.b{i}));
end

% (gamma*FHtx + beta*D'(v+b)) / (gamma*H'H + beta*D'D) -> Fu (in fourier domain)
c = PAR.beta/PAR.gamma;
b = FHtx + c*sum(conj(data.FD).*fft2(data.v+data.b),4);
Fu = b./(data.FHtH + c*data.FDtD);
for i = 1:data.k
    data.Hu(:,i,:) = reshape(real(ifft2(Fu.*repmat(data.FH(:,:,i),[1,1,col]))),[],1,col);
end
data.u = real(ifft2(Fu));

%% solve for v
% soft thresholding of Du-b -> v 
Dub = real(ifft2(data.FD.*repmat(Fu,[1,1,1,size(data.FD,4)]))) - data.b;
Dubn = repmat(sqrt(sum(sum(Dub.^2,4),3)),[1,1,col,size(Dub,4)]);
data.v = softThres(Dub,Dubn);

%% update b
data.b = data.v - Dub;

%% soft thresholding
    function R = softThres(x,n)
        t = PAR.alpha/PAR.beta;
        R = zeros(size(x));
        R(n>t) = x(n>t).*(1-t./n(n>t));
    end
end
