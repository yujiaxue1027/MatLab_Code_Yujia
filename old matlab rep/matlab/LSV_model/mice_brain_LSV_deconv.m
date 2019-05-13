%% mice brain image simulation
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

%% load data
img0 = im2double(rgb2gray(imread('mice_brain.png')));
[rows,cols] = size(img0);

% filter1 = fspecial('gaussian',64,1);
% filter1 = filter1./sum(filter1(:));
% filter2 = fspecial('gaussian',64,16);
% filter2 = filter2./sum(filter2(:));
filter1 = imresize(im2double(rgb2gray(imread('k1.png'))),[64,64]);
filter1 = filter1./sum(filter1(:));
filter2 = imresize(im2double(rgb2gray(imread('k2.png'))),[64,64]);
filter2 = filter2./sum(filter2(:));

filters = [];
filters(:,:,1) = filter1;
filters(:,:,2) = filter2;

% masks = zeros(rows,cols,2);
% [xxx,yyy] = meshgrid(linspace(-0.5,0.5,size(img0,1)));
% masks(:,:,2) = sqrt(xxx.^2+yyy.^2).*1.3;
% masks(:,:,1) = 1 - masks(:,:,2);

masks = rand(rows,cols,2);
masks(:,:,2) = 1-masks(:,:,1);

offset = size(filters,1);
k_size = size(filters,1);
img1 = padarray(img0,[offset,offset]);
y = zeros(size(img1));
for i = 1:size(masks,1)
    for j = 1:size(masks,2)
            y(i+offset-k_size/2:i+offset+k_size/2-1,...
            j+offset-k_size/2:j+offset+k_size/2-1) = ...
            y(i+offset-k_size/2:i+offset+k_size/2-1,...
            j+offset-k_size/2:j+offset+k_size/2-1)+...
            img1(i+offset,j+offset).*(masks(i,j,1).*(filters(:,:,1))+...
            masks(i,j,2).*(filters(:,:,2)));
    end
    disp(i);
end
y = y(offset+1:offset+rows,offset+1:offset+cols);
y = y;
figure,imagesc(y),axis image off, colormap gray;

%% LSI deconv with basis PSFs
test_psf1 = filters(:,:,1);
test_psf2 = filters(:,:,2);
test_psf1 = padarray(test_psf1,0.5*(size(y)-size(test_psf1)));
test_psf2 = padarray(test_psf2,0.5*(size(y)-size(test_psf2)));
miu = 0.01;
recon_tik1 = deconv_tik(y,test_psf1,miu);
recon_tik2 = deconv_tik(y,test_psf2,miu);
figure;
subplot(1,2,1),imagesc(recon_tik1),axis image off,colormap gray;title('Tikhonov deconv 1');
subplot(1,2,2),imagesc(recon_tik2),axis image off,colormap gray;title('Tikhonov deconv 2');


%% LSV deconv using ADMM algorithm (TV constraint)
y1 = y + 0*randn(size(y));
y2 = y + 0.02*randn(size(y));
psfs = padarray(filters,0.5*[size(y1,1)-size(filters,1),...
    size(y1,2)-size(filters,2)]);
M = masks;
para = [];
para.tau = 0.0001;
para.mu1 = 0.5;
para.mu2 = 0.25;
para.color = 'gray';
para.maxiter = 30;
[recon_admm,~] = admm_lsv(y1,psfs,M,para);
para.tau = 0.001;
[recon_admm_noise,~] = admm_lsv(y2,psfs,M,para);
figure,imagesc(recon_admm),axis image off;colormap gray;truesize;
figure,imagesc(recon_admm_noise),axis image off;colormap gray;truesize;

%% LSV deconv using ADMM algorithm (tv l1 nonnega)
y3 = y2;
psfs = padarray(filters,0.5*[size(y2,1)-size(filters,1),...
    size(y2,2)-size(filters,2)]);
M = masks;
para2 = [];
para2.tau1 = 0.01;
para2.tau2 = 0.015;
para2.mu1 = 0.5;
para2.mu2 = 0.25;
para2.mu3 = 0.5;
para2.mu4 = 0.25;
para2.color = 'gray';
para2.maxiter = 30;
[recon_admm2,recon_admm_history2] = admm_lsv2(y2,psfs,M,para2);
figure,imagesc(recon_admm2),axis image off;colormap gray;truesize;
