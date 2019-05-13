%% psf estimate with real measured data
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

y = im2double(imread('exp0516/pattern.tif'));
psf30 = im2double(imread('exp0516/3.tif'));
psf60 = im2double(imread('exp0516/6.tif'));
bg = im2double(imread('exp0516/bg.tif'));
psf30 = abs(psf30-bg);
pattern = generate_calibration_pattern(32,5);
pattern = padarray(pattern,size(pattern))./255;

%% process
obj = pattern;
y2_ = y;
y2 = y2_(994-299:994+300,731-299:731+300);
y3 = y2_(1016-299:1016+300,1309-299:1309+300);
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

% figure,imagesc(obj),axis image off;colormap parula;title('obj');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% [x,y] = ginput(2);
% diag_obj = norm([x(1)-x(2),y(1)-y(2)]);
% figure,imagesc(y2),axis image off;colormap parula;title('noise measure');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% [x,y] = ginput(2);
% diag_y2 = norm([x(1)-x(2),y(1)-y(2)]);
% scale_obj2y2 = diag_y2/diag_obj;
% 
% 
% figure,imagesc(xcorr_obj),axis image off;colormap parula;title('xcorr obj');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% [x,y] = ginput(2);
% diag_obj_x = norm([x(1)-x(2),y(1)-y(2)]);
% figure,imagesc(xcorr_y2),axis image off;colormap parula;title('xcorr noise measure');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% [x,y] = ginput(2);
% diag_y2_x = norm([x(1)-x(2),y(1)-y(2)]);
% scale_obj2y2_x = diag_y2_x/diag_obj_x;
%% 
figure,imagesc(peak_obj),axis image off;colormap parula;title('peak pos obj xcorr');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
vec_obj = [x(2)-x(1),y(2)-y(1)];
diag_obj_xp = norm(vec_obj);
theta_obj = rad2deg(atan(vec_obj(2)/vec_obj(1)));
figure,imagesc(peak_y2),axis image off;colormap parula;title('peak pos measure xcorr');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(2);
vec_y = [x(2)-x(1),y(2)-y(1)];
diag_y2_xp = norm(vec_y);
theta_y = rad2deg(atan(vec_y(2)/vec_y(1)));
scale_obj2y2_xp = diag_y2_xp/diag_obj_xp;
to_rotate = theta_obj-theta_y;
obj2 = imresize(obj,scale_obj2y2_xp);
obj2 = imrotate(obj2,-1*to_rotate,'crop');
y2_pad = make_the_same(y2,obj2);
y3_pad = make_the_same(y3,obj2);
figure,imshow(obj2)
figure,imshow(y2_pad)
figure,imshow(y3_pad)
%% tik psf estimate
k_hat_real = deconv_tik(y2_pad,obj2./sum(obj2(:)),0.0001);
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
para.maxiter= 20;
para.color= 'gray';
[k_hat_admm,~] = admm_recon4(y2_pad,obj2,para);
figure,imagesc(k_hat_admm),axis image off;colormap gray;...
    title('estimated kernel from real measure using admm');
%%
[k_hat_admm2,~] = admm_recon4(y3_pad,obj2,para);
figure,imagesc(k_hat_admm2),axis image off;colormap gray;...
    title('estimated kernel from real measure using admm');



%% mla psf
obj3 = make_the_same(obj2,y2_);
[mla_hat_admm,~] = admm_recon4(y2_,obj3,para);




