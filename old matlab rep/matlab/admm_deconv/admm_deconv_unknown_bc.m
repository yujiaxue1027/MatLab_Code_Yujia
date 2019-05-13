%% admm deconv with unknown b.c.

% define FFT IFFT CONV
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
central_part = @(x) x(1+size(x,1)/3:2*size(x,1)/3,1+size(x,2)/3:2*size(x,2)/3);
F_pad = @(x) fftshift(fft2(ifftshift(padarray(x,size(x)))));
pad3 = @(x) padarray(x,size(x));
% load image
img0 = imread('cameraman.tif');
img0 = double(img0)./255;
% load motion blur kernel
kernel = imread('kernel1.png');
kernel = double(kernel);
factor = 0.125;
kernel = imresize(kernel,factor);
kernel = kernel./sum(kernel(:));
l_img = length(img0);
l_kernel = length(kernel);
kernel = padarray(kernel,repmat(l_img-l_kernel,1,2)./2);
% CONV
measure = circ_conv(img0,kernel);
figure,subplot(1,3,1),imagesc(img0),axis image,...
    colormap gray,title('original');
subplot(1,3,2),imagesc(kernel),axis image,...
    colormap gray, title('kernel');
subplot(1,3,3),imagesc(measure),axis image,...
    colormap gray, title('measure');
% 1 deconv without regularization without zero padding
result_no_reg_no_pad = real(Ft(F(measure)./(eps+F(kernel))));
figure,imagesc(result_no_reg_no_pad),colormap gray,...
    axis image off,title('deconv no reg no pad');
% 2 deconv without regularization with zero padding
result_no_reg_yes_pad = real(Ft(F_pad(measure)./(eps+F_pad(kernel))));
result_no_reg_yes_pad = central_part(result_no_reg_yes_pad);
figure,imagesc(result_no_reg_yes_pad),colormap gray,...
    axis image off,title('deconv no reg yes pad');
%% 3 deconv tikhonov with zero padding
para = 0.02;
max_spectrum = max(max(abs(F_pad(kernel))));
result_l2 = real(Ft(F_pad(measure).*conj(F_pad(kernel))...
    ./(abs(F_pad(kernel)).^2+(para*max_spectrum)^2)));
result_l2 = central_part(result_l2);
figure,imagesc(result_l2),colormap gray,...
    axis image off,title('deconv tikhonov');
%% 4 deconv l1 fista
tk = 1;
x0 = zeros(3*size(img0));
xkm1 = x0;
yk = x0;
xk = x0;
maxiter =100;
alpha = 1;% 0.03 for brain image
thres = 0.00000005;%0.01 for brain image
psf1 = pad3(kernel);
psf2 = flip(flip(pad3(kernel),1),2);
for i = 1:maxiter
    residual = central_part(circ_conv(xk,psf1))-measure;
    nresidual = pad3(residual);
    grad = circ_conv(nresidual,psf2);
    xk = tv_proximal(yk-alpha.*grad,thres);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    figure(51),imagesc((xk)),axis image off;colormap gray;title(i);
    mse = norm(residual(:));
    figure(52),plot(i,mse,'bo'),grid on,hold on;...
        title(['fista: data term, alpha: ',num2str(alpha),' thres: ',num2str(thres)])
end
%% 5 deconv l1 admm
tau = 0.1;
mu = 0.5;
kernel3x = pad3(kernel);
kernel3x_re = flip(flip(kernel3x,1),2);
xt = zeros(3*size(img0));
[u2t,y] = wavedec2(xt,1,'haar');
u1t = circ_conv(xt,kernel3x);
alpha1t = 0;
alpha2t = 0;
maxiter = 50;
mask = mu*ones(3*size(img0));
mask = mask + padarray(ones(size(img0)),size(img0));
for i = 1:maxiter
    tmp = circ_conv(u1t+alpha1t,kernel3x_re)+...
        waverec2(u2t+alpha2t,y,'haar');
    xtp1 = Ft(F(tmp)./(abs(F(kernel3x)).^2+1));
    
    u1tp1 = (pad3(measure)+mu*(circ_conv(xtp1,kernel3x)-alpha1t))./mask;

    [tmp,~] = wavedec2(xtp1,1,'haar');
    u2tp1 = wthresh(tmp-alpha2t,'s',tau/mu);
    
    alpha1tp1 = alpha1t - (circ_conv(xtp1,kernel3x)-u1tp1);
    [tmp,~] = wavedec2(xtp1,1,'haar');
    alpha2tp1 = alpha2t - (tmp-u2tp1);
    
    u2t = u2tp1;
    xt = xtp1;
    u1t = u1tp1;
    alpha1t = alpha1tp1;
    alpha2t = alpha2tp1;
    
    figure(53),imagesc(real(central_part(xt))),axis image off;colormap gray;title(i);
    residual = central_part(circ_conv(xt,kernel3x))-measure;
    mse = norm(residual(:));
    figure(54),plot(i,mse,'bo'),grid on,hold on;...
        title('admm');
end




