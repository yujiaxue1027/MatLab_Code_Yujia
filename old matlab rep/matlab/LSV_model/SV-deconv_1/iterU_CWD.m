function data = iterU_CWD(data, PAR)

col = data.col;

%% solve for u
% b = gamma/beta*M'(x+a) + C'(v+b)
b = PAR.gamma/PAR.beta * squeeze(sum((data.x + data.a) .* repmat(data.M,[1,1,col]),2));
b = b./repmat(data.MtM,[1,col]);
V = data.v + data.b;
Vx = padarray(V(:,1:end-1,:,1),[0,1]);
Vx = Vx(:,2:end,:) - Vx(:,1:end-1,:);
Vy = padarray(V(1:end-1,:,:,2),[1,0]);
Vy = Vy(2:end,:,:) - Vy(1:end-1,:,:);
b = b + reshape(Vx + Vy,[],col);
data.u = real(ifft2(fft2(reshape(b,data.sz)) ./ repmat(data.FICtC,[1,1,col])));

%% solve for v
% soft thresholding of Cu-b -> v
u = data.u;
Cu(:,:,:,1) = cat(2,u(:,1:end-1,:)-u(:,2:end,:),zeros(data.sz(1),1,col));
Cu(:,:,:,2) = cat(1,u(1:end-1,:,:)-u(2:end,:,:),zeros(1,data.sz(2),col));
Dub = Cu - data.b;
Dubn = repmat(sqrt(sum(sum(Dub.^2,4),3)),[1,1,col,size(Dub,4)]);
data.v = softThres(Dub,Dubn);

%% update b
data.b = data.b - Cu + data.v;

%% soft thresholding
    function R = softThres(x,n)
        t = PAR.alpha/PAR.beta;
        R = zeros(size(x));
        R(n>t) = x(n>t).*(1-t./n(n>t));
    end

end
