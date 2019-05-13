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

load data0803.mat
load paras.mat
y1 = im2double(imread('y1.tif')); %2mm
y2 = im2double(imread('y2.tif')); %6mm
y3 = im2double(imread('y3.tif')); %8mm
y4 = im2double(imread('y4.tif')); %no support

y1s = average_shrink(y1,2);
y2s = average_shrink(y2,2);
y3s = average_shrink(y3,2);
y4s = average_shrink(y4,2);

[y1_lsi,y1_lsi_hist] = LSI_model_ADMM_Solver(y1s, psf_on_axis_s, para_lsi);
[y2_lsi,y2_lsi_hist] = LSI_model_ADMM_Solver(y2s, psf_on_axis_s, para_lsi);
[y3_lsi,y3_lsi_hist] = LSI_model_ADMM_Solver(y3s, psf_on_axis_s, para_lsi);
[y4_lsi,y4_lsi_hist] = LSI_model_ADMM_Solver(y4s, psf_on_axis_s, para_lsi);
y1_lsi = crop2d(y1_lsi);
y2_lsi = crop2d(y2_lsi);
y3_lsi = crop2d(y3_lsi);
y4_lsi = crop2d(y4_lsi);

para_lsv.M_dagger_thres = 200; % no clip on vmax
para_lsv.clip_max = 10;
[y1_lsv,y1_lsv_hist] = LSV_model_ADMM_Solver(y1s, bpsf_s, intp_mask_s, para_lsv);
para_lsv.M_dagger_thres = 2;
para_lsv.clip_max = 0.01;
[y2_lsv,y2_lsv_hist] = LSV_model_ADMM_Solver(y2s, bpsf_s, intp_mask_s, para_lsv);
[y3_lsv,y3_lsv_hist] = LSV_model_ADMM_Solver(y3s, bpsf_s, intp_mask_s, para_lsv);
para_lsv.M_dagger_thres = 2;
[y4_lsv,y4_lsv_hist] = LSV_model_ADMM_Solver(y4s, bpsf_s, intp_mask_s, para_lsv);

close all
pitch = 6.6;
x = [-648:1:647]*pitch*2;
y = [-486:1:485]*pitch*2;
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y1_lsi,0,10)),axis image; colormap gray,title('LSI, obj1');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y2_lsi,0,10)),axis image; colormap gray,title('LSI, obj2');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y3_lsi,0,10)),axis image; colormap gray,title('LSI, obj3');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y4_lsi,0,10)),axis image; colormap gray,title('LSI, obj4');

figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y1_lsv,0,10)),axis image; colormap gray,title('LSV, obj1');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y2_lsv,0,10)),axis image; colormap gray,title('LSV, obj2');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y3_lsv,0,10)),axis image; colormap gray,title('LSV, obj3');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(x,y,clip(y4_lsv,0,10)),axis image; colormap gray,title('LSV, obj4');


[y1_LSI,y1_LSI_hist] = LSI_model_ADMM_Solver(y1, psf_on_axis, para_lsi);
[y2_LSI,y2_LSI_hist] = LSI_model_ADMM_Solver(y2, psf_on_axis, para_lsi);
[y3_LSI,y3_LSI_hist] = LSI_model_ADMM_Solver(y3, psf_on_axis, para_lsi);
[y4_LSI,y4_LSI_hist] = LSI_model_ADMM_Solver(y4, psf_on_axis, para_lsi);
y1_LSI = crop2d(y1_LSI);
y2_LSI = crop2d(y2_LSI);
y3_LSI = crop2d(y3_LSI);
y4_LSI = crop2d(y4_LSI);

para_lsv.M_dagger_thres = 200; % no clip on vmax
para_lsv.clip_max = 10;
[y1_LSV,y1_LSV_hist] = LSV_model_ADMM_Solver(y1, bpsf, intp_mask, para_lsv);
para_lsv.M_dagger_thres = 2;
para_lsv.clip_max = 2;
[y2_LSV,y2_LSV_hist] = LSV_model_ADMM_Solver(y2, bpsf, intp_mask, para_lsv);
[y3_LSV,y3_LSV_hist] = LSV_model_ADMM_Solver(y3, bpsf, intp_mask, para_lsv);
para_lsv.M_dagger_thres = 2;
[y4_LSV,y4_LSV_hist] = LSV_model_ADMM_Solver(y4, bpsf, intp_mask, para_lsv);

close all;
pitch = 6.6;
X = [-1296:1:1295]*pitch;
Y = [-972:1:971]*pitch;
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y1_LSI,0,10)),axis image; colormap gray,title('LSI, obj1');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y2_LSI,0,10)),axis image; colormap gray,title('LSI, obj2');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y3_LSI,0,10)),axis image; colormap gray,title('LSI, obj3');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y4_LSI,0,10)),axis image; colormap gray,title('LSI, obj4');

figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y1_LSV,0,10)),axis image; colormap gray,title('LSV, obj1');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y2_LSV,0,10)),axis image; colormap gray,title('LSV, obj2');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y3_LSV,0,10)),axis image; colormap gray,title('LSV, obj3');
figure('rend','painters','pos',[50 50 1500 900]);imagesc(X,Y,clip(y4_LSV,0,10)),axis image; colormap gray,title('LSV, obj4');
