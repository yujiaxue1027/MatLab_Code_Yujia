%% define useful functions
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
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
my_xcorr2 = @(x,y) Ft(F(x).*conj(F(y)));
VEC = @(x) x(:);

%% read in 26 * 26 grid of data
path = '0625R10S51/';
num_scans = 26;
rows = 1944;
cols = 2592;
data = zeros(rows,cols,num_scans,num_scans);
for i = 1:num_scans
    for j = 1:num_scans
        name = [path,num2str(i*2-1,'%.2d'),'_',num2str(j*2-1,'%.2d'),'.tif'];
        data(:,:,i,j) = im2double(imread(name));
    end
    disp(i);
end
disp('finish loading images');
[src_h, src_v] = meshgrid(linspace(-5,5,num_scans));


%% start processing
horizontal_start = 8;
horizontal_end = 23;
horizontal_start_pos = [916, 978]; % [row, col]
horizontal_end_pos = [909, 1866];

vertical_start = 7;
vertical_end = 21;
vertical_start_pos = [554,1272];
vertical_end_pos = [1390, 1277];

drowdi = (vertical_end_pos(1) - vertical_start_pos(1))/(vertical_end-vertical_start);
drowdj = (horizontal_end_pos(1) - horizontal_start_pos(1))/(horizontal_end-horizontal_start);
dcoldi = (vertical_end_pos(2) - vertical_start_pos(2))/(vertical_end-vertical_start);
dcoldj = (horizontal_end_pos(2) - horizontal_start_pos(2))/(horizontal_end-horizontal_start);

center_pos = [917,1277]; %% 13 13
pixel_pos = zeros(2,num_scans, num_scans);
for i = 1:num_scans
    for j = 1:num_scans
        pixel_pos(1,i,j) = center_pos(1) + (i-13)*drowdi + (j-13)*drowdj;
        pixel_pos(2,i,j) = center_pos(2) + (i-13)*dcoldi + (j-13)*dcoldj;
    end
end

%% magnification in pixel/mm
magnification_h = norm(horizontal_start_pos - horizontal_end_pos)/...
                  abs(src_h(13,horizontal_start) - src_h(13,horizontal_end));
magnification_v = norm(vertical_start_pos - vertical_end_pos)/...
                  abs(src_v(vertical_start,13) - src_v(vertical_end,13));              


%% without refine alignment
% for i = 1:num_scans
%     for j = 1:num_scans
%         tmp_img = data(:,:,i,j);
%         tmp_shift = round([center_pos(1)-pixel_pos(1,i,j),...
%                      center_pos(2)-pixel_pos(2,i,j)]);
%         tmp_img = pad2d(tmp_img);
%         tmp_img = circshift(tmp_img, tmp_shift);
%         tmp_img = crop2d(tmp_img);
%         data(:,:,i,j) = tmp_img;
%     end
%     disp(i);
% end
% disp('done');

%% with refine alignment
window_size = 256;
center_patch = data(center_pos(1)-window_size/2:center_pos(1)+window_size/2-1,...
                    center_pos(2)-window_size/2:center_pos(2)+window_size/2-1,13,13);
tot_shift_threshold = 20;
                
for i = 1:num_scans
    for j = 1:num_scans
        tmp_img = data(:,:,i,j);
        tmp_patch = tmp_img(pixel_pos(1,i,j)-window_size/2:pixel_pos(1,i,j)+window_size/2-1,...
                            pixel_pos(2,i,j)-window_size/2:pixel_pos(2,i,j)+window_size/2-1);
        
        
        tmp_xcorr2 = my_xcorr2(center_patch, tmp_patch);
        [tmp_row,tmp_col] = find(tmp_xcorr2==max(tmp_xcorr2(:)));
        tmp_row = tmp_row(1);
        tmp_col = tmp_col(1);
        tmp_row = tmp_row - (window_size/2+1);
        tmp_col = tmp_col - (window_size/2+1);
        
        tot_shift = sqrt(tmp_row^2 + tmp_col^2);
        if tot_shift <= tot_shift_threshold
            tmp_shift = round([center_pos(1)-pixel_pos(1,i,j)+tmp_row,...
                         center_pos(2)-pixel_pos(2,i,j)+tmp_col]);
            disp([num2str(i),' ',num2str(j),' refined']);
        else
            tmp_shift = round([center_pos(1)-pixel_pos(1,i,j),...
                         center_pos(2)-pixel_pos(2,i,j)]);
            disp([num2str(i),' ',num2str(j),' default']);
        end

        tmp_img = pad2d(tmp_img);
        tmp_img = circshift(tmp_img, tmp_shift);
        tmp_img = crop2d(tmp_img);
        data(:,:,i,j) = tmp_img;
    end
    disp(i);
end
disp('done');

%% optional compensate for intensity variation
tot_energy = zeros(26,26,3,4);
window_rows = [360,950,1550];
window_cols = [700, 1300, 1900, 2500];
window_size = 128;
for i = 1:3
    for j =1:4
        for ii = 1:26
            for jj = 1:26
                tmp_patch = data(window_rows(i)-window_size/2:window_rows(i)+window_size/2-1,...
                                 window_cols(j)-window_size/2:window_cols(j)+window_size/2-1,...
                                 ii,jj);
                tot_energy(ii,jj,i,j) = norm(VEC(tmp_patch));
            end
        end
    end
    disp(i);
end
tot_energy = tot_energy./max(tot_energy(:));

compensate_threshold = 0.03;
for i = 1:3
    for j =1:4
        for ii = 1:26
            for jj = 1:26
                tmp_patch = data(window_rows(i)-window_size/2:window_rows(i)+window_size/2-1,...
                                 window_cols(j)-window_size/2:window_cols(j)+window_size/2-1,...
                                 ii,jj);
                if tot_energy(ii,jj,i,j) >= 0.03
                    tmp_patch = tmp_patch./tot_energy(ii,jj,i,j);
                    data(window_rows(i)-window_size/2:window_rows(i)+window_size/2-1,...
                                 window_cols(j)-window_size/2:window_cols(j)+window_size/2-1,...
                                 ii,jj) = tmp_patch;
                end
            end
        end
    end
    disp(i);
end




%% PSF decomposition
data = reshape(data,[rows*cols, num_scans*num_scans]);
num_bpsf = 20;
[U,S,V] = svds(data,num_bpsf);
bpsf = zeros(rows,cols,num_bpsf);
for i = 1:num_bpsf
    bpsf(:,:,i) = reshape(U(:,i),[rows,cols]);
end

%% compute coefficients and compare the measured to estimated
sparse_coefficients = S*V';
sparse_coefficients = reshape(sparse_coefficients,[num_bpsf, num_scans, num_scans]);
clear data;

%% interpolate grid to full fov
[col_dense, row_dense] = meshgrid(-1*cols/2:cols/2-1, -1*rows/2:rows/2-1);
row_dense = row_dense./magnification_v;
col_dense = col_dense./magnification_h;
intp_mask = zeros(rows,cols, num_bpsf);
for i = 1:num_bpsf
    tmp = squeeze(sparse_coefficients(i,:,:));
    intp_mask(:,:,i) = griddata(src_v, src_h, tmp, row_dense, col_dense,'cubic');
    disp(i);
end
intp_mask(isnan(intp_mask)) = 0;

%% read in real measures
y1 = im2double(imread('y2.tif'));
y2 = im2double(imread('y.tif'));

%% tikhonov 
psf_tik = zeros(rows,cols);
for i = 1:size(bpsf,3)
    psf_tik = psf_tik + bpsf(:,:,i)*intp_mask(rows/2,cols/2,i);
end
psf_tik = psf_tik./sum(psf_tik(:));
y1_recon_tik = deconv_tik(y1, psf_tik, 0.0001);
figure,imagesc(col_dense(1,:), row_dense(:,1),y1_recon_tik),axis image; colormap gray;

% %% lsi admm
% para= [];
% para.tau1= 0.001;
% para.tau2= 0.001;
% para.mu1= 0.0050;
% para.mu4= 0.5000;
% para.mu2= 0.0500;
% para.mu3= 0.5000;
% para.mu_inc= 1;
% para.mu_dec= 1;
% para.mu_tol= 2;
% para.maxiter= 30;
% para.color= 'gray';
% [recon_admm_lsi,lsi_hist] = admm_recon4(y,psf_tik,para);

%% lsv admm
para2 = [];
para2.tau1 = 0.1;
para2.tau2 = 0.2;
para2.mu1 = 5;
para2.mu2 = 1;
para2.mu3 = 2.5;
para2.mu4 = 1;
para2.color = 'gray';
para2.maxiter = 30;
[recon_lsv,lsv_hist] = admm_lsv2(y1,bpsf,intp_mask,para2);

%% shrink
y1_s = average_shrink(y1,2);
bpsf_s = zeros(rows/2, cols/2, num_bpsf);
intp_mask_s = zeros(rows/2, cols/2, num_bpsf);
for i = 1:num_bpsf
    bpsf_s(:,:,i) = average_shrink(bpsf(:,:,i), 2);
    intp_mask_s(:,:,i) = average_shrink(intp_mask(:,:,i), 2);
    disp(i);
end
psf_tik_s = average_shrink(psf_tik, 2);
%% process small sized data
%% LSV
para_lsv = [];
para_lsv.tau1 = 0.001;
para_lsv.tau2 = 0.001;
para_lsv.mu1 = 0.005;
para_lsv.mu2 = 0.05;
para_lsv.mu3 = 0.5;
para_lsv.mu4 = 0.5;
para_lsv.color = 'gray';
para_lsv.maxiter = 30;
para_lsv.M_dagger_threshold = 200;
para_lsv.clip_min = 0;
para_lsv.clip_max = 1;
[recon_lsv_s, lsv_hist_s] = admm_lsv2(y1_s, bpsf_s, intp_mask_s, para_lsv);

%% LSI
para= [];
para.tau1= 0.001;
para.tau2= 0.001;
para.mu1= 0.0050;
para.mu2= 0.0500;
para.mu3= 0.5000;
para.mu4= 0.5000;
para.mu_inc= 1;
para.mu_dec= 1;
para.mu_tol= 2;
para.maxiter= 30;
para.color= 'gray';
[recon_lsi_s,lsi_hist_s] = admm_recon4(y1_s,psf_tik_s,para);






