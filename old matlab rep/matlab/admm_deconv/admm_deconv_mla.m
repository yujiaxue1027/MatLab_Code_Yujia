%% admm deconv

% define FFT IFFT CONV
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
central_part = @(x) x(1+size(x,1)/3:2*size(x,1)/3,1+size(x,2)/3:2*size(x,2)/3);
F_pad = @(x) fftshift(fft2(ifftshift(padarray(x,size(x)))));
pad3 = @(x) padarray(x,size(x));
% load image
% load('mla_data.mat');
load('resolution_target_data.mat');
% load('diffuse_cam_data2.mat');
% load('yunzhe_9digits.mat');
% 5 deconv l1 admm
tau = 1.5;% 0.1 for brain image
% mu = 0.5 ;% 0.5 for brain image
mu1 = 0.5;
mu2 = 0.05;
kernel3x = pad3(kernel);
kernel3x_re = flip(flip(kernel3x,1),2);
xt = zeros(3*size(measure));
[u2t,y] = wavedec2(xt,1,'haar');
u1t = circ_conv(xt,kernel3x);
alpha1t = 0;
alpha2t = 0;
maxiter = 30;
mask = mu1*ones(3*size(measure));
mask = mask + padarray(ones(size(measure)),size(measure));
for i = 1:maxiter
    tmp = circ_conv(u1t+alpha1t,kernel3x_re)+...
        waverec2(u2t+alpha2t,y,'haar');
    xtp1 = Ft(F(tmp)./(abs(F(kernel3x)).^2+1));
    
    u1tp1 = (pad3(measure)+mu1*(circ_conv(xtp1,kernel3x)-alpha1t))./mask;

    [tmp,~] = wavedec2(xtp1,1,'haar');
    u2tp1 = wthresh(tmp-alpha2t,'s',tau/mu2);
    
    alpha1tp1 = alpha1t - (circ_conv(xtp1,kernel3x)-u1tp1);
    [tmp,~] = wavedec2(xtp1,1,'haar');
    alpha2tp1 = alpha2t - (tmp-u2tp1);
    
    u2t = u2tp1;
    xt = xtp1;
    u1t = u1tp1;
    alpha1t = alpha1tp1;
    alpha2t = alpha2tp1;
    
    figure(53),subplot(1,2,1),imagesc(real((xt))),axis image off;title(i);
    residual = central_part(circ_conv(xt,kernel3x))-measure;
    mse = norm(residual(:));
    subplot(1,2,2),plot(i,mse,'bo'),grid on,hold on;...
        title('admm');
    drawnow;
end




