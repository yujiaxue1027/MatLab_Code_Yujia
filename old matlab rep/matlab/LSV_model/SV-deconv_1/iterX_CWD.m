function data = iterX_CWD(data)

k = data.k;
col = data.col;

% solve for x
Mu = repmat(data.M,[1,1,col]).*repmat(reshape(data.u,[],1,col),[1,k,1]);
b = data.Htg + Mu - data.a;
b = fft2(reshape(b,[data.sz,k,col]));
x = squeeze(sum(b .* repmat(data.FH,[1,1,1,col]),3));
x = x ./ repmat(data.FIHHt,[1,1,col]);
x = repmat(reshape(x,[data.sz,1,col]),[1,1,k,1]) .* repmat(conj(data.FH),[1,1,1,col]);
ti = real(ifft2(b - x));
data.x = reshape(ti(1:data.sz(1),1:data.sz(2),:,:),[],k,col);

% update a
data.a = data.a - Mu + data.x;

end
