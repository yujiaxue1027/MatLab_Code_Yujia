%% admm deconv with non-negative constraint
% define FFT IFFT CONV
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
central_part = @(x) x(1+size(x,1)/3:2*size(x,1)/3,1+size(x,2)/3:2*size(x,2)/3);
F_pad = @(x) fftshift(fft2(ifftshift(padarray(x,size(x)))));
pad3 = @(x) padarray(x,size(x));
% functions for ulugbek tv proximal 
average_h = @(x) (x+circshift(x,[0,1]))./sqrt(2);
difference_h = @(x) (x-circshift(x,[0,1]))./sqrt(2);
average_v = @(x) (x+circshift(x,[1,0]))./sqrt(2);
difference_v = @(x) (x-circshift(x,[1,0]))./sqrt(2);
average_hT = @(x) (x+circshift(x,[0,-1]))./sqrt(2);
difference_hT = @(x) (x-circshift(x,[0,-1]))./sqrt(2);
average_vT = @(x) (x+circshift(x,[-1,0]))./sqrt(2);
difference_vT = @(x) (x-circshift(x,[-1,0]))./sqrt(2);

% load image
% load('mla_data.mat');
% load('resolution_target_data.mat');
% load('new_res_target.mat');
load('new_res_target2.mat');
% load('diffuse_cam_data2.mat');
% load('yunzhe_9digits.mat');
kernel = kernel./sum(kernel(:));
max_val = 255;
% 5 deconv l1 admm
tau = 2.5;
mu_u = 0.5;
mu_v = 0.5;
mu_w = 0.5;
kernel3x = pad3(kernel);
kernel3x_re = flip(flip(kernel3x,1),2);
xt = zeros(3*size(measure));
ut = W_operator(xt);
vt = circ_conv(xt,kernel3x);
wt = xt;
gamma_ut = zeros(size(ut));
gamma_vt = zeros(size(wt));
gamma_wt = zeros(size(vt));
maxiter = 200;
mask = mu_v*ones(3*size(measure));
mask = mask + padarray(ones(size(measure)),size(measure));
history = [];
if exist('img0')
    img0 = (img0-mean(img0(:)))./std(img0(:));
end


for i = 1:maxiter
    % update u
    utp1 = softthres_d(W_operator(xt)+gamma_ut./mu_u,4*sqrt(2)*tau/mu_u);
    
    % update v
    vtp1 = (pad3(measure)+mu_v*circ_conv(xt,kernel3x)+gamma_vt)./mask;
    
    % update w
    wtp1 = min(max(0,xt+gamma_wt./mu_w),max_val);
    
    % update x
    tmp = WT_operator(mu_u*utp1/mu_v-gamma_ut/mu_v)+...
        circ_conv(vtp1-gamma_vt/mu_v,kernel3x_re)+...
        (mu_w*wtp1/mu_v-gamma_wt/mu_v);
    xtp1 = Ft(F(tmp)./((abs(F(kernel3x))).^2+(mu_u+mu_w)/mu_v));

    % update gammma
    gamma_utp1 = gamma_ut + mu_u*(W_operator(xtp1)-utp1);
    gamma_vtp1 = gamma_vt + mu_v*(circ_conv(xtp1,kernel3x)-vtp1);
    gamma_wtp1 = gamma_wt + mu_w*(xtp1-wtp1);
    
    
    
    % refreshing variables
    ut = utp1;
    vt = vtp1;
    wt = wtp1;
    xt = xtp1;
    gamma_ut = gamma_utp1;
    gamma_vt = gamma_vtp1;
    gamma_wt = gamma_wtp1;

    history(:,:,i) = xt;
    
    if exist('img0')
        xt_n = (xt-mean(xt(:)))./std(xt);
%         psnr = 20*log10(255/sqrt((sum((xt_n(:)-img0(:)).^2)/(768*768))));
        ncc_val = abs(sum(xt_n(:).*img0(:))/sqrt(sum(xt_n(:).*xt_n(:))*sum(img0(:).*img0(:))));
        figure(53),subplot(1,3,1),imagesc(real((xt))),colorbar,axis image off;title(i);
        subplot(1,3,2),plot(i,ncc_val,'ro'),grid on, hold on,...
            title('admm: ncc')
        residual = central_part(circ_conv(xt,kernel3x))-measure;
        mse = norm(residual(:));
        subplot(1,3,3),plot(i,mse,'bo'),grid on,hold on;...
            title('admm: norm of residual');
        drawnow;
    else
        figure(53),subplot(1,2,1),imagesc(real((xt))),colorbar,axis image off;title(i);
        
        residual = central_part(circ_conv(xt,kernel3x))-measure;
        mse = norm(residual(:));
        subplot(1,2,2),plot(i,mse,'bo'),grid on,hold on;...
            title('admm: norm of residual');
        drawnow; 
    end
end

%     tv1 = xt(:,2:end)-xt(:,1:end-1);
%     tv2 = xt(2:end,:)-xt(1:end-1,:);
%     tv = sqrt(sum(tv1(:).^2+tv2(:).^2));
%     subplot(1,3,3),plot(i,tv,'bo'),grid on,hold on;...
%         title('admm: tv of current guess');
%     tmp = circ_conv(u1t+alpha1t,kernel3x_re)+...
%          WT_operator(u2t+alpha2t)+...               
%         (u3t+alpha3t);
%     xtp1 = Ft(F(tmp)./(abs(F(kernel3x)).^2+2));
%     
%     u1tp1 = (pad3(measure)+mu1*(circ_conv(xtp1,kernel3x)-alpha1t))./mask;
% 
%     u2tp1 = W_operator(xtp1)-alpha2t;
%     u2tp1(:,:,1:2) = wthresh(u2tp1(:,:,1:2),'s',tau/mu2);
%     
%     u3tp1 = max(0,xtp1-alpha3t);
%     
%     alpha1tp1 = alpha1t - (circ_conv(xtp1,kernel3x)-u1tp1);
%     alpha2tp1 = alpha2t - (W_operator(xtp1)-u2tp1);
%     alpha3tp1 = alpha3t - (xtp1-u3tp1);
%     
%     u2t = u2tp1;
%     xt = xtp1;
%     u1t = u1tp1;
%     u3t = u3tp1;
%     alpha1t = alpha1tp1;
%     alpha2t = alpha2tp1;
%     alpha3t = alpha3tp1;



