function data = iterX_RWD(data)

k = data.k;

% solve for x
for i = 1:data.col
    data.x(:,:,i) = squeeze(sum(reshape(data.iMtM .* repmat(reshape(data.Mtg(:,:,i) + data.Hu(:,:,i)-data.a(:,:,i),[],1),[1,k]),[],k,k), 2));
end

% update a
data.a = data.a - data.Hu + data.x;

end
