function err = err2orig(u,data,off,idx)

u = u(1:size(data.g,1),1:size(data.g,2),:);
r = (data.sPsf - 1) / 2;
o = data.og(off*r(1)+1:end-off*r(1),off*r(2)+1:end-off*r(2),:);
u = u(off*r(1)+1:end-off*r(1),off*r(2)+1:end-off*r(2),:);
% m = repmat(data.msk(off*r(1)+1:end-off*r(1),off*r(2)+1:end-off*r(2)),[1,1,size(o,3)]);
if nargin<4
    idx = [1,size(u,2)];
end
u = u(:,idx(1):sum(idx)-1);
o = o(:,idx(1):sum(idx)-1);

err = sum((u(:)-o(:)).^2)/numel(o);
% err = sum(m(:).*(u(:)-o(:)).^2)/sum(m(:));
err = 10*log10(1/err);
