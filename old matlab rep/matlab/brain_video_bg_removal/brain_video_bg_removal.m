%% process brain video background removal

%% mean frame
v = VideoReader('video.mp4');
frames = read(v,[1 Inf]);
frames = double(squeeze(frames(:,:,2,:)));
[rows,cols,num_frames] = size(frames);
mean_frame = mean(frames,3);
figure,imagesc(mean_frame),axis image off;colormap gray;colorbar;
new_frames = frames - repmat(mean_frame,[1,1,121]);

%% SVD principle frame
data = reshape(frames,rows*cols,num_frames);
[U,S,V] = svds(data,1);
principle_frame = reshape(U,rows,cols);
figure,imagesc(principle_frame),axis image off;colormap gray;colorbar;
new_frames2 = frames;
for i = 1:num_frames
    vec = frames(:,:,i);
    vec = vec(:);
    new_frames2(:,:,i) = frames(:,:,i) - dot(vec,U)*principle_frame;
end

%% SVD principle frame 3
data = reshape(frames,rows*cols,num_frames);
[U,S,V] = svds(data,3);
principle_frame1 = reshape(U(:,1),rows,cols);
principle_frame2 = reshape(U(:,2),rows,cols);
figure,imagesc(principle_frame1),axis image off;colormap gray;colorbar;
figure,imagesc(principle_frame2),axis image off;colormap gray;colorbar;
new_frames3 = frames;
for i = 1:num_frames
    vec = frames(:,:,i);
    vec = vec(:);
    new_frames3(:,:,i) = frames(:,:,i) - dot(vec,U(:,1))*principle_frame1...
        - dot(vec,U(:,2))*principle_frame2;
end
%% RPCA inexact alm
low_rank_frame = zeros(rows,cols);
sparse_data = zeros(size(frames));
for i = 1:8
    for j = 1:8
        tmp = frames(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:);
        tmp = reshape(tmp,174*128,num_frames);
        [A_hat, E_hat ,~] = inexact_alm_rpca(tmp,0.01, 1e-7, 100);
        low_rank_frame(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174) = reshape(A_hat(:,1),128,174);
        sparse_data(1+(i-1)*128:128+(i-1)*128,...
            1+(j-1)*174:174+(j-1)*174,:) = reshape(E_hat,128,174,num_frames);
    end
    disp(i);
end
figure,imagesc(low_rank_frame),axis image off,colormap gray

%% display
min_val = min(new_frames(:));
max_val = max(new_frames(:));
min_val_svd = min(new_frames2(:));
max_val_svd = max(new_frames2(:));
min_val_svd2 = min(new_frames3(:));
max_val_svd2 = max(new_frames3(:));
min_val_rpca = min(sparse_data(:));
max_val_rpca = max(sparse_data(:));

for i = 1:size(new_frames,3)
figure(1);
subplot(3,2,1),imagesc(abs(new_frames(:,:,i))),caxis([0,max_val]);colorbar;
colormap gray;axis image off;title(['frame: ',num2str(i)]);
subplot(3,2,2),imagesc((new_frames(:,:,i))),caxis([min_val,max_val]);colorbar;
colormap gray;axis image off;
subplot(3,2,3),imagesc(abs(new_frames2(:,:,i))),caxis([0,max_val_svd]);colorbar;
colormap gray;axis image off;
subplot(3,2,4),imagesc((new_frames2(:,:,i))),caxis([min_val_svd,max_val_svd]);colorbar;
colormap gray;axis image off;
subplot(3,2,5),imagesc(abs(sparse_data(:,:,i))),caxis([0,max_val_rpca]);colorbar;
colormap gray;axis image off;
subplot(3,2,6),imagesc((sparse_data(:,:,i))),caxis([min_val_rpca,max_val_rpca]);colorbar;
colormap gray;axis image off;
pause(0.1);
drawnow;
end