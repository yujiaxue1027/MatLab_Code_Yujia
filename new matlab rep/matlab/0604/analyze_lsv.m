%% script to analyze LSV of mla

%% read in all data
path = '0606/';
num_scans = 21;
rows = 1944;
cols = 2592;
data = zeros(rows,cols,num_scans,num_scans);
for i = 1:num_scans
    for j = 1:num_scans
        name = [path,num2str(i,'%.2d'),'_',num2str(j,'%.2d'),'.tif'];
        data(:,:,i,j) = im2double(imread(name));
    end
    disp(i);
end
disp('finish loading images');
[src_h, src_v] = meshgrid(linspace(-6,6,num_scans));

%% start processing
k_size = 64;
psf_stack = zeros(k_size,k_size,num_scans,num_scans,3,3);
centers = zeros(2,num_scans,num_scans,3,3);% first row second col
window = 256;
center_pos = [961, 1300];
right_pos = [960, 1682];
down_pos = [1343, 1304];
drowdi = (down_pos(1) - center_pos(1))/((num_scans-1)/2);
drowdj = (right_pos(1) - center_pos(1))/((num_scans-1)/2);
dcoldi = (down_pos(2) - center_pos(2))/((num_scans-1)/2);
dcoldj = (right_pos(2) - center_pos(2))/((num_scans-1)/2);

down_pos2 = [1470, 1312];
right_pos2 = [954, 1806];
drowdii = down_pos2(1) - center_pos(1);
drowdjj = right_pos2(1) - center_pos(1);
dcoldii = down_pos2(2) - center_pos(2);
dcoldjj = right_pos2(2) - center_pos(2);

for ii = 1:3
    for jj = 1:3
        for i = 1:num_scans
            for j = 1:num_scans
                center_r = min(max(center_pos(1)...
                    +(i-11)*drowdi+(j-11)*drowdj...
                    +(ii-2)*drowdii+(jj-2)*drowdjj,1),rows);
                center_c = min(max(center_pos(2)...
                    +(i-11)*dcoldi+(j-11)*dcoldj...
                    +(ii-2)*dcoldii+(jj-2)*dcoldjj,1),cols);
                
                rs = round(min(max(center_r - window/2,1),rows));
                re = round(min(max(center_r + window/2-1,1),rows));
                cs = round(min(max(center_c - window/2,1),cols));
                ce = round(min(max(center_c + window/2-1,1),cols));
                patch = data(rs:re,cs:ce,i,j);
                [tmp_x,tmp_y] = meshgrid([cs:ce],[rs:re]);
                centers(1,i,j,ii,jj) = min(max(round(sum(sum(patch.*tmp_y))/sum(patch(:))),1),rows);
                centers(2,i,j,ii,jj) = min(max(round(sum(sum(patch.*tmp_x))/sum(patch(:))),1),cols);
                if (centers(1,i,j,ii,jj) - k_size/2)>0 && (centers(2,i,j,ii,jj) - k_size/2)>0 && ...
                        (centers(1,i,j,ii,jj) + k_size/2 -1)<=rows &&...
                        (centers(2,i,j,ii,jj) + k_size/2 -1)<=cols
                    psf_stack(:,:,i,j,ii,jj) = data(centers(1,i,j,ii,jj)-k_size/2:centers(1,i,j,ii,jj)+k_size/2-1,...
                        centers(2,i,j,ii,jj)-k_size/2:centers(2,i,j,ii,jj)+k_size/2-1,i,j);
                end
            end
        end
    end
    disp(ii);
end
disp('done');

%% display psf stack
cmax = 0.7*max(psf_stack(:));
cmin = min(psf_stack(:));
figure('rend','painters','pos',[50 50 1500 900]);
for i = 1:num_scans
    for j = 1:num_scans
        tmp_fig = zeros(num_scans);
        tmp_fig(i,j) = 1;
        subplot(3,4,5);
        imagesc(padarray(tmp_fig,[5,5])),axis image off; 
        title(['Horizontal: ',num2str(src_h(i,j)),'mm, Vertical: ',num2str(src_v(i,j)),'mm']);
        
        subplot(3,4,2);
        imagesc(psf_stack(:,:,i,j,1,1)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,3);
        imagesc(psf_stack(:,:,i,j,1,2)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,4);
        imagesc(psf_stack(:,:,i,j,1,3)),axis image off;colormap hot;%caxis([cmin,cmax]);
        
        subplot(3,4,6);
        imagesc(psf_stack(:,:,i,j,2,1)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,7);
        imagesc(psf_stack(:,:,i,j,2,2)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,8);
        imagesc(psf_stack(:,:,i,j,2,3)),axis image off;colormap hot;%caxis([cmin,cmax]);
        
        subplot(3,4,10);
        imagesc(psf_stack(:,:,i,j,3,1)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,11);
        imagesc(psf_stack(:,:,i,j,3,2)),axis image off;colormap hot;%caxis([cmin,cmax]);
        subplot(3,4,12);
        imagesc(psf_stack(:,:,i,j,3,3)),axis image off;colormap hot;%caxis([cmin,cmax]);
        
%         saveas(gcf,['figs/',num2str(i),'_',num2str(j),'.tif']);
        pause(0.01);
    end
end
%% save shift variant psf for visualization
for i = 1:3
    for j = 1:3
        img = zeros(rows,cols);
        for ii = 1:num_scans
            for jj=1:num_scans
                tmp_row = centers(1,ii,jj,i,j);
                tmp_col = centers(2,ii,jj,i,j);
                img(tmp_row-k_size/2:tmp_row+k_size/2-1, ...
                    tmp_col-k_size/2:tmp_col+k_size/2-1) ...
                    = psf_stack(:,:,ii,jj,i,j);
            end
        end
        [J,~] = gray2ind(img);
        imwrite(J,parula,['lsv_',num2str(i),num2str(j),'.tif']);
    end
end

img = reshape(data,rows,cols,[]);
img = mean(img,3);
imagesc(img);
img = img./max(img(:));
[J,~] = gray2ind(img);
imwrite(J*4,parula,['lsv_psf.tif']);
%% correlation analysis
correlation_map = zeros(num_scans,num_scans,3,3);
for i = 1:3
    for j = 1:3
        psf0 = psf_stack(:,:,11,11,i,j);
        psf0 = psf0./norm(psf0(:));
        for ii = 1:num_scans
            for jj=1:num_scans
                psf1 = psf_stack(:,:,ii,jj,i,j);
                psf1 = psf1./norm(psf1(:));
                correlation_map(ii,jj,i,j) = max(max(xcorr2(psf0,psf1)));
            end
        end
    end
end
[tgt_h, tgt_v] = meshgrid(linspace(-6,6,num_scans*50));
interpolated_xcorr_map = zeros(num_scans*50,num_scans*50,3,3);
for i = 1:3
    for j=1:3
        interpolated_xcorr_map(:,:,i,j) = griddata(src_h,src_v,correlation_map(:,:,i,j),tgt_h, tgt_v,'v4');
        [J,~] = gray2ind(interpolated_xcorr_map(:,:,i,j));
        imwrite(J,parula,['xcorr_map_',num2str(i),num2str(j),'.tif']);
    end
end
%% separation analysis
dist_h = zeros(num_scans);
dist_v = zeros(num_scans);
for i = 1:num_scans
    for j = 1:num_scans
        dist_h(i,j) = mean([
                          norm(centers(:,i,j,1,2)-centers(:,i,j,1,1)),...
                          norm(centers(:,i,j,1,3)-centers(:,i,j,1,2)),...
                          norm(centers(:,i,j,2,2)-centers(:,i,j,2,1)),...
                          norm(centers(:,i,j,2,3)-centers(:,i,j,2,2)),...
                          norm(centers(:,i,j,3,2)-centers(:,i,j,3,1)),...
                          norm(centers(:,i,j,3,3)-centers(:,i,j,3,2))]);
        dist_v(i,j) = mean([
                          norm(centers(:,i,j,2,1)-centers(:,i,j,1,1)),...
                          norm(centers(:,i,j,3,1)-centers(:,i,j,2,1)),...
                          norm(centers(:,i,j,2,2)-centers(:,i,j,1,2)),...
                          norm(centers(:,i,j,3,2)-centers(:,i,j,2,2)),...
                          norm(centers(:,i,j,2,3)-centers(:,i,j,1,3)),...
                          norm(centers(:,i,j,3,3)-centers(:,i,j,2,3))]);
    end
end
average_separation_horizontal = mean(dist_h(:));
average_separation_vertical = mean(dist_v(:));


%% PSF decomposition
% psf_mat = reshape(psf_stack,k_size*k_size,[]);
tmp_stack= permute(psf_stack,[3,4,1,2,5,6]);
psf_mat = zeros(num_scans,num_scans,k_size,k_size*3*3);
for i = 1:num_scans
    for j = 1:num_scans
        psf_mat(i,j,:,:) = reshape(squeeze(tmp_stack(i,j,:,:,:,:)),64,[]);
    end
    disp(i);
end
psf_mat = permute(psf_mat,[3,4,1,2]);
psf_mat = reshape(psf_mat,k_size*k_size*3*3,num_scans*num_scans);
[U,S,V] = svd(psf_mat);
singular_vals = diag(S);
num_bpsf = 20;
bpsf = zeros(k_size,k_size*3*3,num_bpsf);
for i = 1:num_bpsf
    bpsf(:,:,i) = reshape(U(:,i),[k_size,k_size*3*3]);
end
figure('rend','painters','pos',[50 50 1500 900]);
cmax = max(bpsf(:));
cmin = min(bpsf(:));
tmp_display = zeros(k_size*num_bpsf,k_size*3*3);
for i = 1:num_bpsf
    tmp_display((i-1)*k_size+1:i*k_size,:) = bpsf(:,:,i);
%     subplot(10,2,i);
%     imagesc(reshape(bpsf(:,:,i),k_size,k_size*3*3)),axis image;colormap hot;caxis([cmin,cmax])
%     title(['basis psf ',num2str(i)]);
end
imagesc(tmp_display'),axis image;colormap hot;caxis([cmin,cmax]);

%% compute coefficients
rows = 1944;
cols = 2592;
mask0 = zeros(num_scans,num_scans,num_bpsf);
decompose_mat = U(:,1:num_bpsf);
for i = 1:num_scans
    for j = 1:num_scans
        psf2decompose = reshape(squeeze(psf_stack(:,:,i,j,:,:)),k_size,k_size*3*3);
        psf2decompose = reshape(psf2decompose,[],1);
        mask0(i,j,:) = inv(decompose_mat'*decompose_mat)*...
                    decompose_mat'*psf2decompose; 
    end
    disp(i);
end

%% randomly display true psf and estimated psf
figure('rend','painters','pos',[50 50 1500 900]);
for trial = 1:4
    i = randi(num_scans);
    j = randi(num_scans);
    subplot(4,2,2*trial-1);
    tmp = psf_stack(:,:,i,j,:,:);
    cmin = min(tmp(:));
    cmax = max(tmp(:));
    imagesc(reshape(squeeze(psf_stack(:,:,i,j,:,:)),k_size,k_size*3*3)),...
        axis equal off,colormap hot;title('true psf');caxis([cmin,cmax]);
    subplot(4,2,2*trial);
    coeff = mask0(i,j,:);
    est = reshape(decompose_mat*coeff(:),[k_size,k_size*3*3]);
    imagesc(est),axis equal off,colormap hot;title('estimated psf');caxis([cmin,cmax]);
end
clear data;


%% interpolate grid to full fov
central_psf_co = centers(:,:,:,2,2);
row_co_31 = squeeze(central_psf_co(1,:,:));
col_co_31 = squeeze(central_psf_co(2,:,:));
[col_target,row_target] = meshgrid(1:cols,1:rows);
intp_mask = zeros(rows,cols,num_bpsf);
for i = 1:num_bpsf
    tmp = mask0(:,:,i);
    intp_mask(:,:,i) = griddata(row_co_31,col_co_31,tmp,row_target,col_target,'cubic');
    disp(i);
end
intp_mask(isnan(intp_mask)) = 0;

%% generate mla basis psfs
bpsf_bank = zeros(rows,cols,num_bpsf);
for k = 1:num_bpsf
    tmp = zeros(rows,cols);
    for i = 1:3
        for j = 1:3
            basis_patch = bpsf(:,64*((i+(j-1)*3)-1)+1:64*((i+(j-1)*3)-1)+64,k);
            basis_patch = padarray(basis_patch,0.5*[rows-k_size,cols-k_size]);
            tmp = tmp + circshift(basis_patch,[(i-2)*648-(j-2)*37,(j-2)*646+(i-2)*28]);
        end
    end
    bpsf_bank(:,:,k) = tmp;
    disp(k);
end



%% simulate with resolution target and mice brain image
img0 = im2double(rgb2gray(imread('target.jpg')));
img0 = imresize(img0,1.0);
img0 = padarray(img0,0.5*[rows-size(img0,1),cols-size(img0,2)]);
y = zeros(size(img0));
for i = 1:size(bpsf_bank,3)
    y = y + conv2d(img0.*intp_mask(:,:,i),bpsf_bank(:,:,i));
    disp(i);
end
figure,imagesc(y),axis image off;truesize,colormap gray;

%% tikhonov solution with central psf
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

psf_tik = zeros(rows,cols);
for i = 1:size(bpsf_bank,3)
    psf_tik = psf_tik + bpsf_bank(:,:,i)*intp_mask(rows/2,cols/2,i);
end
psf_tik = psf_tik./sum(psf_tik(:));
recon_tik = deconv_tik(y,psf_tik,0.00001);
figure,imagesc(recon_tik),axis image off;truesize,colormap gray;
%% lsi admm
para= [];
para.tau1= 0.001;
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
[recon_admm_lsi,lsi_hist] = admm_recon4(y,psf_tik,para);

%% lsv admm
para2 = [];
para2.tau1 = 0.01;
para2.tau2 = 0.015;
para2.mu1 = 0.5;
para2.mu2 = 0.25;
para2.mu3 = 0.5;
para2.mu4 = 0.25;
para2.color = 'gray';
para2.maxiter = 30;
[recon_lsv_large,lsv_large_hist] = admm_lsv2(y,bpsf_bank,intp_mask,para2);


%% shrink image for faster processing
y_small = average_shrink(y,2);
bpsf_bank_small = zeros(rows/2,cols/2,size(bpsf_bank,3));
masks_small = zeros(rows/2,cols/2,size(bpsf_bank,3));
for i = 1:size(bpsf_bank,3)
    bpsf_bank_small(:,:,i) = average_shrink(bpsf_bank(:,:,i),2);
    masks_small(:,:,i) = average_shrink(intp_mask(:,:,i),2);
    disp(i);
end
[recon_lsv_small,lsv_small_hist] = admm_lsv2(y_small,bpsf_bank_small,masks_small,para2);

%% mice brain image
img0 = im2double((imread('newbrain.png')));
img0 = imresize(img0(1:end-1,:,2),0.64);
img0 = padarray(img0,0.5*[rows-size(img0,1),cols-size(img0,2)]);
y2 = zeros(size(img0));
for i = 1:size(bpsf_bank,3)
    y2 = y2 + conv2(img0.*intp_mask(:,:,i),bpsf_bank(:,:,i),'same');
    disp(i);
end
figure,imagesc(y2),axis image off;truesize,colormap gray;

%% tikhonov solution with central psf
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

psf_tik = zeros(rows,cols);
for i = 1:size(bpsf_bank,3)
    psf_tik = psf_tik + bpsf_bank(:,:,i)*intp_mask(rows/2,cols/2,i);
end
psf_tik = psf_tik./sum(psf_tik(:));
recon_tik2 = deconv_tik(y2,psf_tik,0.00001);
figure,imagesc(recon_tik2),axis image off;truesize,colormap gray;
%% lsi admm
para= [];
para.tau1= 0.001;
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
[recon_admm_lsi2,lsi_hist2] = admm_recon4(y2,psf_tik,para);

%% lsv admm
para2 = [];
para2.tau1 = 0.01;
para2.tau2 = 0.015;
para2.mu1 = 0.5;
para2.mu2 = 0.25;
para2.mu3 = 0.5;
para2.mu4 = 0.25;
para2.color = 'gray';
para2.maxiter = 30;
[recon_lsv_large2,lsv_large_hist2] = admm_lsv2(y2,bpsf_bank,intp_mask,para2);


%% shrink image for faster processing
y_small2 = average_shrink(y2,2);
bpsf_bank_small = zeros(rows/2,cols/2,size(bpsf_bank,3));
masks_small = zeros(rows/2,cols/2,size(bpsf_bank,3));
for i = 1:size(bpsf_bank,3)
    bpsf_bank_small(:,:,i) = average_shrink(bpsf_bank(:,:,i),2);
    masks_small(:,:,i) = average_shrink(intp_mask(:,:,i),2);
    disp(i);
end
[recon_lsv_small2,lsv_small_hist2] = admm_lsv2(y_small2,bpsf_bank_small,masks_small,para2);

%% real measurement
% real measurement of resolution target
y_real = im2double(imread('real_medium.tif'));
[psf_shift_len,psf_shift_angle] = find_copy_separation(psf_tik);
[y_shift_len,y_shift_angle] = find_copy_separation(y_real);
mag = psf_shift_len/y_shift_len;
y_post = imresize(y_real,mag);
y_post = imrotate(y_post,psf_shift_angle-y_shift_angle,'crop');
y_post = make_the_same(y_post,psf_tik);
[post_len,post_ang] = find_copy_separation(y_post);
[recon_real,real_hist] = admm_lsv2(y_post,bpsf_bank,intp_mask,para2);




