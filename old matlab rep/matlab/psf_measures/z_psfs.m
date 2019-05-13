%% z psfs
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
reversed = @(x) flip(flip(x,1),2);
myxcorr2 = @(x,y) crop2d(Ft(F(pad2d(x)).*F(pad2d(reversed(y)))));
path = 'z/';
num_imgs = 19;
rows = 1944;
cols = 2592;
data0 = zeros(rows,cols,num_imgs);
for i = 1:num_imgs
    name = [path,num2str(i),'.tif'];
    data0(:,:,i) = im2double(imread(name));
end
data = wthresh(data0(861:1060,1201:1400,:),'h',0.1);
ref = data(:,:,1);
data_2 = data0;
tmp = myxcorr2(ref,ref);
[~,tmp2] = max(tmp(:));
[refr,refc] = ind2sub(size(tmp),tmp2);
for i = 1:num_imgs
    img = data(:,:,i);
    R = myxcorr2(ref,img);
    [~,tmp2] = max(R(:));
    [r,c] = ind2sub(size(R),tmp2);
    new = circshift(data0(:,:,i),[r-refr,c-refc]);
    data_2(:,:,i) = new;
end

%% 
max_proj_xy = squeeze(max(data_2,[],3));   
max_proj_xy = medfilt2(max_proj_xy,[5,5]);
figure,imagesc(max_proj_xy),axis image off;colormap hot;
max_proj_xz = squeeze(max(data_2,[],1));
max_proj_xz = medfilt2(max_proj_xz,[5,5]);
max_proj_xz = imresize(max_proj_xz,[cols,20*num_imgs]);
figure,imagesc(max_proj_xz),axis image off;colormap hot;
max_proj_yz = squeeze(max(data_2,[],2)); 
max_proj_yz = medfilt2(max_proj_yz,[5,5]);   
max_proj_yz = imresize(max_proj_yz,[rows,20*num_imgs]);
figure,imagesc(max_proj_yz),axis image off;colormap hot;

%%
psf0 = data_2(:,:,1);
psf1 = data_2(:,:,10);
psf2 = data_2(:,:,19);
tmp = psf0 + psf1+ psf2;
figure,imagesc(tmp),axis image off;colormap hot;
psf0 = repmat(psf0,[1,1,3]);
psf0(:,:,2:3) = 0;
psf1 = repmat(psf1,[1,1,3]);
psf1(:,:,[1,3]) = 0;
psf2 = repmat(psf2,[1,1,3]);
psf2(:,:,1:2) = 0;
tmp2 = psf0 + psf1+ psf2;
figure,imagesc(tmp2),axis image off;truesize;

%%
data_3 = zeros(size(data_2));
data_3(:,:,1) = data_2(:,:,1);
data_3(:,:,10) = data_2(:,:,10);
data_3(:,:,19) = data_2(:,:,19);
max_proj_xz2 = squeeze(max(data_3,[],1));
% max_proj_xz2 = medfilt2(max_proj_xz2,[5,5]);
max_proj_xz2 = imresize(max_proj_xz2,[cols,20*num_imgs]);
figure,imagesc(max_proj_xz2),axis image off;colormap hot;

%%
cr = 933;
cc = 1286;
r = 32;
xsec_h = squeeze(data_2(cr,cc-r+1:cc+r,:));
xsec_v = squeeze(data_2(cr-r+1:cr+r,cc,:));
x_co = [-32:31].*2.2;
data2plot = xsec_h(:,[1,3,5,7,9,11,13,15,17,19]);
u = linspace(14,7,10);
f = 3.3;
v = f.*u./(u-f);
mag = v./u;
x_ref = x_co./mag(6);
for i = 1:6
    data2plot(:,i) = interp1(x_co./mag(i),data2plot(:,i),x_ref);
end


figure;
for i = 1:2:11
    plot(x_co,squeeze(xsec_h(:,i))),hold on;
end
grid on;
legend('- 1.0 mm','- 0.8 mm','- 0.6 mm','- 0.4 mm','- 0.2 mm','0.0 mm');


data2plot = data2plot./max(data2plot(:));
figure;
for i = 1:10
    plot(x_co,squeeze(data2plot(:,i))),hold on;
end
grid on;
legend('- 1.0 mm','- 0.8 mm','- 0.6 mm','- 0.4 mm','- 0.2 mm','0.0 mm'...
    ,'0.2 mm','0.4 mm','0.6 mm', '0.8 mm');

toexcel = data2plot./max(data2plot(:));
x_ref = x_ref(:);
   
    
    
    
    
    
    
    
    
    
    
    