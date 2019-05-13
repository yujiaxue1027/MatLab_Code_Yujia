%% mla psf
rows = 1944;
cols = 2592;
num_imgs = 16;
imgs = zeros(rows,cols,num_imgs);
for i = 1:num_imgs
    imgs(:,:,i) = im2double(imread([num2str(i),'.tif']));
end

%% coarsely label center psf
pos = zeros(num_imgs,2);
for i = 1:num_imgs
    figure(1),imagesc(imgs(:,:,i)),colormap jet,axis image off;
    tmp = ginput(1);
    pos(i,1) = round(tmp(2));
    pos(i,2) = round(tmp(1));
end

%% compute center of central psf
cpsf = zeros(101,101,num_imgs);
for i = 1:num_imgs
    patch = imgs(pos(i,1)-50:pos(i,1)+50,...
                       pos(i,2)-50:pos(i,2)+50,...
                       i);
    [x,y] = meshgrid(-50:1:50);
    x = x + pos(i,2);
    y = y + pos(i,1);
    xbar = round(sum(x(:).*patch(:))./sum(patch(:)));
    ybar = round(sum(y(:).*patch(:))./sum(patch(:)));
    cpsf(:,:,i) = imgs(ybar-50:ybar+50,...
                       xbar-50:xbar+50,...
                       i);
end

%% calibrate relative location
rlt_locn = zeros(8,2);
for i = 1:8
    figure(2),imagesc(imgs(:,:,1)),axis image off;colormap jet;
    tmp = ginput(1);
    rlt_locn(i,1) = tmp(2);
    rlt_locn(i,2) = tmp(1);
end
rlt_locn = round(rlt_locn - repmat(pos(1,:),[8,1]));

%% compute center of other psfs
psfs = zeros(101,101,8,num_imgs);
for j = 1:8
    for i = 1:num_imgs
    patch = imgs(pos(i,1)-50+rlt_locn(j,1):pos(i,1)+50+rlt_locn(j,1),...
                 pos(i,2)-50+rlt_locn(j,2):pos(i,2)+50+rlt_locn(j,2),...
                       i);
    [x,y] = meshgrid(-50:1:50);
    x = x + pos(i,2)+rlt_locn(j,2);
    y = y + pos(i,1)+rlt_locn(j,1);
    xbar = round(sum(x(:).*patch(:))./sum(patch(:)));
    ybar = round(sum(y(:).*patch(:))./sum(patch(:)));
    psfs(:,:,j,i) = imgs(ybar-50:ybar+50,...
                       xbar-50:xbar+50,...
                       i);
    end
end

%% save results
data = zeros(101,101,3,3,num_imgs);
data(:,:,1,1,:) = psfs(:,:,1,:);
data(:,:,1,2,:) = psfs(:,:,2,:);
data(:,:,1,3,:) = psfs(:,:,3,:);
data(:,:,2,1,:) = psfs(:,:,4,:);
data(:,:,2,2,:) = cpsf;
data(:,:,2,3,:) = psfs(:,:,5,:);
data(:,:,3,1,:) = psfs(:,:,6,:);
data(:,:,3,2,:) = psfs(:,:,7,:);
data(:,:,3,3,:) = psfs(:,:,8,:);
border = 20;
for i = 1:num_imgs
    output = ones(3*101+4*border, 3*101+4*border);
    for ii = 0:2
        for jj = 0:2
            iistart = ii*(101+border)+border;
            jjstart = jj*(101+border)+border;
            output(iistart+1:iistart+101, jjstart+1:jjstart+101) = data(:,:,ii+1,jj+1,i);
        end
    end
    imwrite(uint8(255*output),'hot',['mla_psf_',num2str(i),'.tif']);
end









    
    