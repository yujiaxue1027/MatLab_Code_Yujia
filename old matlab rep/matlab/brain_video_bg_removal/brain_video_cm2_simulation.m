%% Brain video CM2 simulation

%% func definition
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

%% load brain video
load('video.mat','video');
[rows,cols,num_frames] = size(video);

%% RPCA on brain video, for comparison
low_rank_frame = zeros(rows,cols);
sparse_data = zeros(size(video));
for i = 1:8
    for j = 1:8
        tmp = video(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:);
        tmp = reshape(tmp,174*128,num_frames);
        [A_hat, E_hat ,~] = inexact_alm_rpca(tmp,0.01, 1e-7, 100);
        low_rank_frame(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174) = reshape(A_hat(:,1),128,174);
        sparse_data(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:) = reshape(E_hat,128,174,num_frames);
    end
    disp([num2str(i),'/8 complete']);
end

%% synthesize CM2 measurements
load('sample_psf.mat','psf');
cm2_video = zeros(rows,cols,num_frames);
for i = 1:num_frames
    cm2_video(:,:,i) = conv2d(video(:,:,i),psf);
    disp(['frame: ',num2str(i),' done.'])
end

%% Without background removal
frame_of_interest = 64;
y = cm2_video(:,:,frame_of_interest);
para1 = [];
para1.tau1= 0.01;
para1.tau2= 0.01;
para1.mu1= 0.0050;
para1.mu4= 0.5000;
para1.mu2= 0.0500;
para1.mu3= 0.5000;
para1.mu_inc= 1;
para1.mu_dec= 1;
para1.mu_tol= 2;
para1.maxiter= 30;
para1.color= 'gray';
[est1,~] = admm_recon4(y,psf,para1);
figure,imagesc(est1),axis image off;colormap gray;
title('direct admm deconv with measured frame');


%% RPCA analysis
low_rank_frame_measure = zeros(rows,cols);
sparse_data_measure = zeros(size(cm2_video));
for i = 1:8
    for j = 1:8
        tmp = cm2_video(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:);
        tmp = reshape(tmp,174*128,num_frames);
        [A_hat, E_hat ,~] = inexact_alm_rpca(tmp,0.01, 1e-7, 100);
        low_rank_frame_measure(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174) = reshape(A_hat(:,1),128,174);
        sparse_data_measure(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:) = reshape(E_hat,128,174,num_frames);
    end
    disp([num2str(i),'/8 complete']);
end
figure,imagesc(low_rank_frame_measure),axis image off;colormap gray;
title('RPCA background estimation from measured video');
% figure,imagesc(conv2d(pad2d(low_rank_frame),psf2x)),axis image off;colormap gray;
figure,imagesc(sparse_data_measure(:,:,frame_of_interest)),axis image off;colormap gray;
title('sparse decomposition of one frame from measured video');
figure,imagesc(sparse_data(:,:,frame_of_interest)),axis image off;colormap gray;
title('sparse decomposition of original video (ground truth)')
figure,imagesc(conv2d(sparse_data(:,:,frame_of_interest),psf));
axis image off;colormap gray;
title('convolution of original sparse frame with mla psf')
est_neuron_video = zeros(size(sparse_data));
for i = 1:size(est_neuron_video,3)
    tmp = sparse_data_measure(:,:,i);
    est_neuron_video(:,:,i) = deconv_tik(tmp,psf,0.00001);
    disp([num2str(i),'/',num2str(num_frames),' frames done.']);
end
est_neuron_video(est_neuron_video<=0.02) = 0;

gt = abs(sparse_data);
gt = gt./max(gt(:));
my_est = est_neuron_video./max(est_neuron_video(:));
mean_neuron = mean(gt,3);
figure,imagesc(mean_neuron),axis image off,colormap jet;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(6);
x = round(x);
y = round(y);
figure;
for i = 1:6
    subplot(2,3,i),plot(squeeze(gt(y(i),x(i),:)),'b'),hold on;
    plot(squeeze(my_est(y(i),x(i),:)),'r');
end

%% compute standard deviation map
sd_orignal_video = std(video,0,3);
sd_orignal_cm = std(cm2_video,0,3);
sd_sparse_oracle_neuron = std(sparse_data,0,3);
sd_sparse_measure = std(sparse_data_measure,0,3);
%%
deconv_sd = deconv_tik(sd_sparse_measure,psf,0.0001);
deconv_sd(deconv_sd<=0) = 0;
mask = (deconv_sd>=0.0032);
figure,imshowpair(sd_sparse_oracle_neuron,mask,'montage');
est_neuron_video2 = zeros(size(sparse_data));
para_new = [];
para_new.mu = 1;
para_new.color = 'gray';
para_new.maxiter = 10;
num_frames_to_process = size(est_neuron_video2,3);
% num_frames_to_process = 20;
tic;
for i = 1:num_frames_to_process
    tmp = sparse_data_measure(:,:,i);
    [est_neuron_video2(:,:,i),~] = admm_oracle_location(tmp,psf,mask,para_new,0);
    disp([num2str(i),'/',num2str(num_frames_to_process),' frames done.']);
end
est_neuron_video2(est_neuron_video2<=0) = 0;
processed_time = toc;
disp(['Total processing time for ',num2str(num_frames_to_process),...
    ' frames: ',num2str(processed_time ),' secs']);

%% compare result
est_neuron_video2 = est_neuron_video2./max(est_neuron_video2(:));
figure,imagesc(sd_sparse_oracle_neuron),axis image off,colormap jet;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
[x,y] = ginput(6);
x = round(x);
y = round(y);
for i = 1:6
    subplot(2,3,i);
    plot(squeeze(gt(y(i),x(i),:)),'b');hold on;
    plot(squeeze(est_neuron_video2(y(i),x(i),:)),'r');
end

% para2 = para1;
% para2.tau1 = 0.0000000000001;
% para2.mu2 = 0.000000001;
% para2.tau2 = 0.1;
% para2.mu4 = 0.5;
% frame = sparse_data_measure(:,:,frame_of_interest);
% frame = frame./max(frame(:));
% [est2,~] = admm_recon4(frame,psf,para2);
% figure,imagesc(est2),axis image off;colormap gray;
% title('admm deconv of decomposed sparse frame');
% y = cm2_video(:,:,frame_of_interest);
% frame2 = y - dot(y(:),low_rank_frame_measure(:))/norm(y(:))/norm(low_rank_frame_measure(:))...
%     *low_rank_frame_measure;
% frame2 = frame2./max(frame2(:));
% [est3,~] = admm_recon4(frame2,psf,para2);
% figure,imagesc(est3),axis image off;colormap gray;
% title('admm deconv of (frame - bg)');

