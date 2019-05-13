function h = cylinderblur(s,t,is)
vec = @(x) x(:);
% h = motionblur(s,k,w,t,is)
% Generate cylindrical blur of diameter s and shifted by t
% rendered in a window of size is


% default window size is [s, s]
if nargin < 3
    is = [s s];
end

[X Y] = meshgrid([0:is(2)]-is(2)/2 + t(2), [0:is(1)]-is(1)/2 + t(1));
I = [vec(X(1:end-1,1:end-1)),vec(X(1:end-1,2:end)),...
        vec(Y(1:end-1,1:end-1)),vec(Y(2:end,1:end-1))];
I2 = [vec(X(1:end-1,1:end-1)).^2+vec(Y(1:end-1,1:end-1)).^2,...
    vec(X(1:end-1,1:end-1)).^2+vec(Y(2:end,2:end)).^2,...
    vec(X(2:end,2:end)).^2+vec(Y(2:end,2:end)).^2,...
    vec(X(2:end,2:end)).^2+vec(Y(1:end-1,1:end-1)).^2];
par = [s/2]^2;
h = reshape(sum(I2 < par,2),size(X)-1);
indc = find(I(:,1)<=0 & I(:,2)>=0 & I(:,3)<=0 & I(:,4)>=0);
h(indc) = 1;
ind = find(h==1 | h==2 | h==3).';
for p = ind
    i = I(p,:);
    h(p) = dblquad(@afun, i(1),i(2),i(3),i(4),1e-1,[],par);
end
h(h==4) = 1;
h = h/sum(h(:));

function r = afun(x,y,p)
r = (x.^2+y.^2)<=p;
