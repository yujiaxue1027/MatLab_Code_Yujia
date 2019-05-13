%% implement fista to solve the crop-conv inverse problem
% load data
load('fista_data.mat');
% define fft and ifft
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% initialize fista
tk = 1;
x0 = zeros(3150,3150);
xkm1 = x0;
yk = x0;
xk = x0;
maxiter =100;
alpha = 0.003;
thres = 0.003;
psf2 = flip(flip(mla_psf,1),2);
psf2 = padarray(psf2,[575,575]);
mselist = [];
% if decrease_flag is 1, learning rate will decrease when residual is small
decrease_flag = 0;
back_length = 8;
decrease_thres = 1;
decrease_rate = 0.5;
%%main iteration here:
for i = 1:maxiter
    residual = forward_cac(xk,mla_psf,2000)-measure;
    nresidual = padarray(residual,[575,575]);% (3150 - 2000)/2= 575;
    grad = conv2(nresidual,psf2,'same');
    xk = myproximal(yk-alpha.*grad,thres);
%     xk = mytruncate(xk,0,1);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    figure(51),imagesc(xk),axis image off;colormap jet;title(i);
    mse = norm(residual(:));
    mselist = [mselist,mse];
    figure(52),plot(i,mse,'bo'),grid on,hold on;title(['data term, alpha: ',num2str(alpha),' thres: ',num2str(thres)])
    if decrease_flag
    if i>=40
        recentmse = mselist(i-back_length:i);
        p = polyfit([1:back_length+1],recentmse,1);
        if abs(p(1))<decrease_thres
            alpha = alpha*decrease_rate;
            decrease_thres = decrease_thres*decrease_rate;
            disp(['at ',num2str(i),'th iteration, decrease alpha, alpha: ',num2str(alpha)]);
        end
    end
    end
end