%% LSV model deconvolution with ADMM
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

%% pre-define two basis PSFs and generated synthetic data
img0 = im2double(imread('barbara.png'));
img0 = cut_edge(img0,30);
[rows,cols] = size(img0);
filter1 = im2double(rgb2gray(imread('k1.png')));
filter1 = filter1./sum(filter1(:));
filter2 = im2double(rgb2gray(imread('k2.png')));
filter2 = filter2./sum(filter2(:));
filters(:,:,1) = filter1;
filters(:,:,2) = filter2;
masks = zeros(rows,cols,2);
masks(:,1:cols/2,1) = 1;
% masks(:,1:cols/2,2) = 0.1;
% masks(:,1+cols/2:end,1) = 0.3;
masks(:,1+cols/2:end,2) = 1;
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
end
y = y(offset+1:offset+rows,offset+1:offset+cols);
figure,imagesc(y),axis image off, colormap gray;

%% LSI deconv with basis PSFs
test_psf1 = 1*filters(:,:,1) + 0*filters(:,:,2);
test_psf2 = 0*filters(:,:,1) + 1*filters(:,:,2);
test_psf3 = squeeze(mean(filters,3));
test_psf1 = padarray(test_psf1,0.5*(size(y)-size(test_psf1)));
test_psf2 = padarray(test_psf2,0.5*(size(y)-size(test_psf2)));
test_psf3 = padarray(test_psf3,0.5*(size(y)-size(test_psf3)));
miu = 0.01;
recon_tik1 = deconv_tik(y,test_psf1,miu);
recon_tik2 = deconv_tik(y,test_psf2,miu);
recon_tik3 = deconv_tik(y,test_psf3,miu);
figure;
subplot(1,3,1),imagesc(recon_tik1),axis image off,colormap gray;title('Tikhonov deconv 1');
subplot(1,3,2),imagesc(recon_tik2),axis image off,colormap gray;title('Tikhonov deconv 2');
subplot(1,3,3),imagesc(recon_tik3),axis image off,colormap gray;title('Tikhonov deconv 3');
test_psf4 = test_psf1(:,129:384);
test_psf5 = test_psf2(:,129:384);
y4 = y(:,1:256);
y5 = y(:,257:end);
miu = 0.01;
recon_tik4 = deconv_tik(y4,(test_psf4),miu);
recon_tik5 = deconv_tik(y5,(test_psf5),miu);
figure;
subplot(1,2,1),imagesc(recon_tik4),axis image off,colormap gray;title('Tikhonov deconv 4');
subplot(1,2,2),imagesc(recon_tik5),axis image off,colormap gray;title('Tikhonov deconv 5');

%% LSV deconv using ADMM algorithm (TV constraint)
y = y;
psfs = padarray(filters,0.5*[size(y,1)-size(filters,1),...
    size(y,2)-size(filters,2)]);
M = masks;
para = [];
para.tau = 0.0001;
para.mu1 = 0.5;
para.mu2 = 0.05;
para.color = 'gray';
para.maxiter = 20;
[recon_admm,recon_admm_history] = admm_lsv(y,psfs,M,para);
figure,imagesc(recon_admm),axis image off;truesize,colormap gray;














