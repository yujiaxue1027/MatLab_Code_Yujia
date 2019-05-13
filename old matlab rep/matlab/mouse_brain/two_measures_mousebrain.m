%% MOUSE BRAIN: script on testing two measurements recon

%% part 1
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,size(x)/2);
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(x,psf) crop2d(Ft((F(pad2d(x)).*F(pad2d(psf)))));

% load('mousebrain.mat','img0');
img0 = (im2double(imread('newbrain2.bmp')));
img0 = img0(:,:,2);
% img0 = imresize(img0,[512,512]);
obj1 = img0./max(img0(:));
imagesc(obj1),axis image off;colormap gray;truesize;
obj1 = pad2d(obj1);

load('simu_psf.mat');
y1 = conv2d(obj1,psf);
figure,imagesc(y1),axis image off;truesize;colormap gray;
para= [];
para.tau1= 1.0000e-03;
para.tau2= 0.001000;
para.mu1= 0.0050;
para.mu4= 0.5000;
para.mu2= 0.0500;
para.mu3= 0.5000;
para.mu_inc= 1;
para.mu_dec= 1;
para.mu_tol= 2;
para.maxiter= 55;
para.color= 'gray';
[recon1,~] = admm_recon4(y1,psf,para);


%% part2
scale = 1.9;
obj2 = imresize(obj1,scale);
[rows,cols] = size(obj1);
[nrows,ncols] = size(obj2);
tmp1 = round((nrows-rows)/2);
tmp2 = nrows-rows-tmp1;
tmp3 = round((ncols-cols)/2);
tmp4 = ncols-cols-tmp3;
psf_ = padarray(padarray(psf,[tmp1,tmp3],'pre'),...
    [tmp2,tmp4],'post');
if mod(size(obj2,1),2)==1
    obj2 = obj2(1:end-1,:);
    psf_ = psf_(1:end-1,:);
end
if mod(size(obj2,2),2)==1
    obj2 = obj2(:,1:end-1);
    psf_ = psf_(:,1:end-1);
end
y2 = conv2d(obj2,psf_);
y2 = y2(tmp1+1:tmp1+rows,tmp3+1:tmp3+cols);
figure,imagesc(y2),axis image off,colormap gray;truesize;
para.maxiter = 40;
[recon2,~] = admm_recon4(y2,psf,para);

%% part3
psfs = zeros(512,512,2);
ys = zeros(512,512,2);
scale = 1.9;
obj2 = imresize(obj1,scale);
[rows,cols] = size(obj1);
[nrows,ncols] = size(obj2);
tmp1 = round((nrows-rows)/2);
tmp2 = nrows-rows-tmp1;
tmp3 = round((ncols-cols)/2);
tmp4 = ncols-cols-tmp3;
psf_ = padarray(padarray(psf,[tmp1,tmp3],'pre'),...
    [tmp2,tmp4],'post');
if mod(size(obj2,1),2)==1
    obj2 = obj2(1:end-1,:);
    psf_ = psf_(1:end-1,:);
end
if mod(size(obj2,2),2)==1
    obj2 = obj2(:,1:end-1);
    psf_ = psf_(:,1:end-1);
end
y2 = conv2d(obj2,psf_);
y2 = y2(tmp1+1:tmp1+rows,tmp3+1:tmp3+cols);
figure,imagesc(y2),axis image off,colormap gray;truesize;

ys(:,:,1) = y2;
psfs(:,:,1) = psf;

psf2 = circshift(psf,[50,0]);
psf_ = padarray(padarray(psf2,[tmp1,tmp3],'pre'),...
    [tmp2,tmp4],'post');
if mod(size(psf_,1),2)==1
    psf_ = psf_(1:end-1,:);
end
if mod(size(psf_,2),2)==1
    psf_ = psf_(:,1:end-1);
end
y2 = conv2d(obj2,psf_);
y2 = y2(tmp1+1:tmp1+rows,tmp3+1:tmp3+cols);
figure,imagesc(y2),axis image off,colormap gray;truesize;
ys(:,:,2) = y2;
psfs(:,:,2) = psf2;

[recon3,~] = admm_recon5(ys,psfs,para);



%% part4
psfs = zeros(512,512,2);
ys = zeros(512,512,2);
scale = 1.9;
obj2 = imresize(obj1,scale);
[rows,cols] = size(obj1);
[nrows,ncols] = size(obj2);
tmp1 = round((nrows-rows)/2);
tmp2 = nrows-rows-tmp1;
tmp3 = round((ncols-cols)/2);
tmp4 = ncols-cols-tmp3;
psf_ = padarray(padarray(psf,[tmp1,tmp3],'pre'),...
    [tmp2,tmp4],'post');
if mod(size(obj2,1),2)==1
    obj2 = obj2(1:end-1,:);
    psf_ = psf_(1:end-1,:);
end
if mod(size(obj2,2),2)==1
    obj2 = obj2(:,1:end-1);
    psf_ = psf_(:,1:end-1);
end
y2 = conv2d(obj2,psf_);
y2 = y2(tmp1+1:tmp1+rows,tmp3+1:tmp3+cols);
figure,imagesc(y2),axis image off,colormap gray;truesize;

ys(:,:,1) = y2;
psfs(:,:,1) = psf;

psf2 = imresize(psf,1.1);
[psf_r,psf_c] = size(psf);
[psf_r2,psf_c2] = size(psf2);
r_start = round((psf_r2-psf_r)/2);
c_start = round((psf_c2-psf_c)/2);
psf2 = psf2(r_start+1:r_start+psf_r,c_start+1:c_start+psf_c);

psf_ = padarray(padarray(psf2,[tmp1,tmp3],'pre'),...
    [tmp2,tmp4],'post');
if mod(size(psf_,1),2)==1
    psf_ = psf_(1:end-1,:);
end
if mod(size(psf_,2),2)==1
    psf_ = psf_(:,1:end-1);
end
y2 = conv2d(obj2,psf_);
y2 = y2(tmp1+1:tmp1+rows,tmp3+1:tmp3+cols);
figure,imagesc(y2),axis image off,colormap gray;truesize;
ys(:,:,2) = y2;
psfs(:,:,2) = psf2;

[recon4,~] = admm_recon5(ys,psfs,para);

% save two_measures_recon_mousebrain




