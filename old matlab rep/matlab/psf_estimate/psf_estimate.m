%% simulation: estimate kernel from measurement and known object
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
deconv_tik = @(y,psf,miu) crop2d(real(Ft(F(pad2d(y)).*...
    conj(F(pad2d(psf)))./...
    ((abs(F(pad2d(psf)))).^2+miu))));
reversed = @(x) flip(flip(x,1),2);
myxcorr2 = @(x,y) crop2d(Ft(F(pad2d(x)).*F(pad2d(reversed(y)))));

obj = generate_calibration_pattern(32,5);
obj = padarray(obj,size(obj))./255;
kernel = im2double(imread('k2.png'));
kernel = imresize(rgb2gray(kernel),[64,64]);
kernel = kernel./sum(kernel(:));
y0 = conv2(obj,kernel,'same');
figure;
subplot(1,3,1),imagesc(obj),axis image off;colormap gray;title('obj');
subplot(1,3,2),imagesc(kernel),axis image off;colormap gray;title('GT kernel');
subplot(1,3,3),imagesc(y0),axis image off; colormap gray;title('noise free measure');


%%
noise_level = 0.01;
y1 = y0 + noise_level.*randn(size(y0));

k_hat0 = deconv_tik(y0,obj,0.00001);
figure,imagesc(k_hat0),axis image off;colormap gray;...
    title('estimated kernel Tik from noise free measure');
k_hat1 = deconv_tik(y1,obj,100);
figure;
subplot(1,2,1),imagesc(y1),axis image off;colormap gray;title('noisy measure')
subplot(1,2,2),imagesc(k_hat1),axis image off;colormap gray;...
    title('estimated kernel Tik from noisy measure');

%%
y2_ = im2double((imread('14_23.tif')));
roi = [901:1300];
y2 = y2_(roi,roi);
figure,subplot(2,2,1),imagesc(y2),axis image off;colormap gray;
subplot(2,2,2),imagesc(obj),axis image off; colormap gray;
xcorr_obj = myxcorr2(obj,obj);
xcorr_y2 = myxcorr2(y2,y2);
subplot(2,2,3),imagesc(xcorr_y2),axis image off;colormap gray;
subplot(2,2,4),imagesc(xcorr_obj),axis image off;colormap gray;

peak_obj = zeros(size(xcorr_obj));
peak_y2 = zeros(size(xcorr_y2));
for i = 2:size(xcorr_obj,1)-1
    for j = 2:size(xcorr_obj,2)-1
        if xcorr_obj(i,j)>=xcorr_obj(i-1,j-1) &&...
                xcorr_obj(i,j)>=xcorr_obj(i-1,j) &&...
                xcorr_obj(i,j)>=xcorr_obj(i-1,j+1) &&...
                xcorr_obj(i,j)>=xcorr_obj(i,j-1) &&...
                xcorr_obj(i,j)>=xcorr_obj(i,j+1) &&...
                xcorr_obj(i,j)>=xcorr_obj(i+1,j-1) &&...
                xcorr_obj(i,j)>=xcorr_obj(i+1,j) &&...
                xcorr_obj(i,j)>=xcorr_obj(i+1,j+1) 
            peak_obj(i,j) = 1;
        end
    end
end
for i = 2:size(xcorr_y2,1)-1
    for j = 2:size(xcorr_y2,2)-1
        if xcorr_y2(i,j)>=xcorr_y2(i-1,j-1) &&...
                xcorr_y2(i,j)>=xcorr_y2(i-1,j) &&...
                xcorr_y2(i,j)>=xcorr_y2(i-1,j+1) &&...
                xcorr_y2(i,j)>=xcorr_y2(i,j-1) &&...
                xcorr_y2(i,j)>=xcorr_y2(i,j+1) &&...
                xcorr_y2(i,j)>=xcorr_y2(i+1,j-1) &&...
                xcorr_y2(i,j)>=xcorr_y2(i+1,j) &&...
                xcorr_y2(i,j)>=xcorr_y2(i+1,j+1) 
            peak_y2(i,j) = 1;
        end
    end
end
figure,imagesc(peak_obj),axis image off,colormap gray;
figure,imagesc(peak_y2),axis image off,colormap gray;

figure,imagesc(obj),axis image off;colormap parula;title('obj');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_obj = norm([x(1)-x(2),y(1)-y(2)]);
figure,imagesc(y2),axis image off;colormap parula;title('noise measure');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_y2 = norm([x(1)-x(2),y(1)-y(2)]);
scale_obj2y2 = diag_y2/diag_obj;


figure,imagesc(xcorr_obj),axis image off;colormap parula;title('xcorr obj');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_obj_x = norm([x(1)-x(2),y(1)-y(2)]);
figure,imagesc(xcorr_y2),axis image off;colormap parula;title('xcorr noise measure');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_y2_x = norm([x(1)-x(2),y(1)-y(2)]);
scale_obj2y2_x = diag_y2_x/diag_obj_x;

figure,imagesc(peak_obj),axis image off;colormap parula;title('peak pos obj xcorr');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_obj_xp = norm([x(1)-x(2),y(1)-y(2)]);
figure,imagesc(peak_y2),axis image off;colormap parula;title('peak pos measure xcorr');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
diag_y2_xp = norm([x(1)-x(2),y(1)-y(2)]);
scale_obj2y2_xp = diag_y2_xp/diag_obj_xp;

obj2 = imresize(obj,scale_obj2y2_xp);
y2_pad = make_the_same(y2,obj2);

%% tik psf estimate
k_hat_real = deconv_tik(y2_pad,obj2,10000);
figure;
subplot(1,2,1),imagesc(y2_pad),axis image off;colormap gray;title('noisy measure')
subplot(1,2,2),imagesc(k_hat_real),axis image off;colormap gray;...
    title('estimated kernel from real measure');
%% admm psf estimate
para= [];
para.tau1= 0.0001;
para.tau2= 0.001;
para.mu1= 0.0050;
para.mu4= 0.5000;
para.mu2= 0.0500;
para.mu3= 0.5000;
para.mu_inc= 1;
para.mu_dec= 1;
para.mu_tol= 2;
para.maxiter= 30;
para.color= 'gray';
[k_hat_admm,~] = admm_recon4(y2_pad,obj2,para);
figure,imagesc(k_hat_admm),axis image off;colormap gray;...
    title('estimated kernel from real measure using admm');

[k_hat_admm_simu,~] = admm_recon4(y1,obj,para);
figure,imagesc(k_hat_admm_simu),axis image off;colormap gray;...
    title('estimated kernel from simu measure using admm');

%% mla psf
obj3 = make_the_same(obj2,y2_);
[mla_hat_admm,~] = admm_recon4(y2_,obj3,para);








