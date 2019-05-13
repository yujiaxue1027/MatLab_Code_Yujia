%% gradient of operation (convolution followed by cropping)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
img = zeros(512,512);
% psf = randn(512,512);
psf = double(imread('psf.tiff'));
psf = psf(25:2024,225:2224);
psf = imresize(psf,[512,512]);
psf = psf./(sum(psf(:)));
figure,imagesc(psf),axis image off;truesize,colormap gray,title('psf');
y = double(imread('obj.tiff'));
y = y(25:2024,225:2224);
y = imresize(y,[512,512]);
y = y./(max(y(:)));
figure,imagesc(y),axis image off; truesize,colormap gray;title('after cropping and convolution');
% after_conv = conv2(img,psf,'same');
% figure,imagesc(after_conv),axis image off; truesize,colormap gray;title('after convolution');
% y = after_conv;
% figure,imagesc(y),axis image off; truesize,colormap gray;title('after cropping and convolution');

%% deconv
input = padarray(y,size(y));
kernel = psf;
kernel = padarray(kernel,size(kernel));
figure,imagesc(input),axis image off; truesize,colormap gray;title('zero padded input');
figure,imagesc(kernel),axis image off; truesize,colormap gray;title('zero padded psf');
reg = 0.005;  %%decrease reg if you see copies, increase if if cannot tell the obj
f_mea = F(input);
f_psf = F(kernel);
para = reg*max(abs(f_psf(:)));
f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
result = real(Ft(f_res));
% result = mytruncate(result,0,1);
result = result(round(size(result,1)/2)-size(img,1)/2:round(size(result,1)/2)+size(img,1)/2-1,...
    round(size(result,1)/2)-size(img,1)/2:round(size(result,1)/2)+size(img,1)/2-1);
figure,imagesc(result),axis image off, colormap gray;colorbar;title('deconv result');

%% implement fista to solve the crop-conv inverse problem
%%initalize parameters for fista
tk = 1;
x0 = zeros(size(img));
xkm1 = x0;
yk = x0;
xk = x0;
maxiter = 500;%100
alpha = 1;%0.005
thres = 0.0001;%0.001
psf2 = flip(flip(psf,1),2);
mselist = [];
back_length = 8;
decrease_thres = 0.01;
decrease_rate = 0.2;
% psf2 = padarray(psf2,[575,575]);
%%main iteration here:
for i = 1:maxiter
    residual = conv2(xk,psf,'same')-y;
%     nresidual = padarray(residual,[575,575]);% (3150 - 2000)/2= 575;
    nresidual = residual;
    grad = conv2(nresidual,psf2,'same');
    xk = myproximal(yk-alpha.*grad,thres);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    figure(51),imagesc(xk),axis image off;truesize,colormap gray;title(i);
    mse = norm(residual(:));
    figure(52),plot(i,mse,'bo'),grid on,hold on;title(['data term, alpha: ',num2str(alpha),' thres: ',num2str(thres)])
    mselist = [mselist,mse];
    if i>=40
        recentmse = mselist(i-back_length:i);
        p = polyfit([1:back_length+1],recentmse,1);
        if abs(p(1))<decrease_thres
            alpha = alpha*decrease_rate;
            decrease_thres = decrease_thres*decrease_rate;
            disp(['at ',num2str(i),'th iteration, decrease alpha, alpha: ',num2str(alpha)]);
        end
    end
    pause(1);
end

