%% script to simulate LSV model and inverse problem
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
deconv2d = @(y,psf,miu) crop2d(real(Ft(F(pad2d(y)).*...
    conj(F(pad2d(psf)))./...
    ((F(pad2d(psf))).^2+miu))));
    
%% data preparation
img0 = im2double(imread('cameraman.tif'));
[rows,cols] = size(img0);
block_rows = 8;
block_cols = 8;
sub_rows = rows/block_rows;
sub_cols = cols/block_cols;
[x,y] = meshgrid(linspace(-0.5,0.5,256));
masks = zeros(256,256,block_rows,block_cols);
for i = 1:block_rows
    for j = 1:block_cols
        tmp = padarray(ones(sub_rows,sub_cols),...
            [rows-sub_rows,cols-sub_cols],'post');
        masks(:,:,i,j) = circshift(tmp,...
            [(i-1)*sub_rows,(j-1)*sub_cols]);
    end
end
clear i j tmp;
% kernels = zeros(size(masks));
% for i = 1:block_rows
%     for j = 1:block_cols
%         tmp = fspecial('gaussian',256,...
%             1*sqrt((i-4.5)^2+(j-4.5)^2));
%         tmp = tmp./sum(tmp(:));
%         kernels(:,:,i,j) = tmp;
%     end
% end
% clear i j tmp;
% y = zeros(rows,cols);
% for i = 1:block_rows
%     for j = 1:block_cols
%         y = y + conv2d(img0.*masks(:,:,i,j),kernels(:,:,i,j));
%     end
% end
% offset = 50;
% img1 = padarray(img0,[offset,offset]);
% measure = zeros(size(img1));
% kernel_size = 51;
% half_size = (kernel_size-1)/2;
% kernels = zeros(kernel_size,kernel_size,rows,cols);
% for i = 1:rows
%     for j = 1:cols
%         tmp = fspecial('gaussian',51,...
%             4*sqrt(x(i,j)^2+y(i,j)^2)+0.5);
%         tmp = tmp./sum(tmp(:));
%         kernels(:,:,i,j) = tmp;
%         measure(i+offset-half_size:i+offset+half_size,...
%             j+offset-half_size:j+offset+half_size) = ...
%             measure(i+offset-half_size:i+offset+half_size,...
%             j+offset-half_size:j+offset+half_size)+...
%             img1(i+offset,j+offset)*tmp;
%     end
%     disp(i);
% end
% clear i j tmp;
% measure = measure(offset+1:offset+rows,offset+1:offset+cols);
% figure,imagesc(measure),axis image off, colormap gray;

% %% inverse
% miu = 0.5;
% % recon1_tik = deconv2d(measure,...
% %     padarray(padarray(kernels(:,:,129,129),[102,102],'pre'),...
% %     [103,103],'post'),miu);
% % figure,imagesc(recon1_tik),axis image off, colormap gray;
% % title('LSI deconv with central PSF');
% recon2_tik = deconv2d(measure,...
%     padarray(padarray(mean(mean(kernels,4),3),[102,102],'pre'),...
%     [103,103],'post'),miu);
% figure,imagesc(recon2_tik),axis image off, colormap gray;
% title('LSI deconv with average PSF');

%% pre-define two basis kernels
k_size = 64;
[x,y] = meshgrid(linspace(-1,1,k_size));
k1 = exp(-200*(x.^2+y.^2));
k2 = exp(-400*((sqrt(x.^2+y.^2)-0.1).^2));
psf_mat = [k1(:),k2(:)];
[U,S,V] = svd(psf_mat);
bpsf1 = reshape(U(:,1),[k_size,k_size]);
bpsf2 = reshape(U(:,2),[k_size,k_size]);

offset = 50;
img1 = padarray(img0,[offset,offset]);
measure = zeros(size(img1));
[xx,yy] = meshgrid(linspace(-0.5,0.5,256));
for i = 1:rows
    for j = 1:cols
        a = 0.1+sqrt(xx(i,j)^2+yy(i,j)^2);
        tmp = a*bpsf1 + (1-a)*bpsf2;
        tmp = tmp./sum(tmp(:));
        measure(i+offset-k_size/2:i+offset+k_size/2-1,...
            j+offset-k_size/2:j+offset+k_size/2-1) = ...
            measure(i+offset-k_size/2:i+offset+k_size/2-1,...
            j+offset-k_size/2:j+offset+k_size/2-1)+...
            img1(i+offset,j+offset)*tmp;
    end
    disp(i);
end
clear i j tmp;
measure = measure(offset+1:offset+rows,offset+1:offset+cols);
figure,imagesc(measure),axis image off, colormap gray;

%% inverse
miu = 0.01;
bpsf1_deconv = padarray(bpsf1,[96,96]);
bpsf1_deconv = bpsf1_deconv./sum(bpsf1_deconv(:));
bpsf2_deconv = padarray(bpsf2,[96,96]);
bpsf2_deconv = bpsf2_deconv./sum(bpsf2_deconv(:));
recon1_tik = deconv2d(measure,bpsf1_deconv,miu);
figure,imagesc(recon1_tik),axis image off, colormap gray;
title('LSI deconv with basis psf 1');
recon2_tik = deconv2d(measure,bpsf2_deconv,miu);
figure,imagesc(recon2_tik),axis image off, colormap gray;
title('LSI deconv with basis psf 2');   
% recon3_tik = deconv2d(measure,0.9*bpsf2_deconv+0.1*bpsf1_deconv,miu);
% figure,imagesc(recon3_tik),axis image off, colormap gray;
% title('LSI deconv with basis psf 2'); 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        