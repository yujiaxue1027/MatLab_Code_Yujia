function h = genBase(D,k)

[u,~,~] = svd(D,'econ');
sh = sum(abs(u(:,1:k)));
h = u(:,1:k)./repmat(sh,size(D,1),1);
