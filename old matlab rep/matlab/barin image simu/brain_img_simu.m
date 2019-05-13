%% reconstruction test from simulated psf and measurements
%%yujia xue 20171118
% 
% pitch = 2; %%um
obj = rgb2gray(imread('brain.png'));
scale = 1;
% obj = double(imresize((obj),scale));
obj = double(imresize((obj),[256,256]));
obj = obj./255;
% if size(obj,1)<=256
% obj = padarray(obj,[round((256-size(obj,1))/2),...
%     round((256-size(obj,1))/2)]);
% obj = imresize(obj,[256,256]);
% end

figure,imagesc(obj),axis image off,colormap gray,colorbar,title('true object');

%% single micro lens psf simulation
% NA = 
single_f_size = 25;
single_f_sigma = 1;
single_psf = fspecial('gaussian',[single_f_size,single_f_size],single_f_sigma);
single_psf = single_psf./(sum(single_psf(:)));
figure,surf(single_psf),shading flat,title('psf of single lenslet');

%% try convolve with single psf
sng_conv = myconv(obj,single_psf,'same');
figure,imagesc(sng_conv),axis image off, colormap gray,colorbar,title('conv with single lens psf');


%% copy psf to mla psf
spacing = 50;
copies = [3,3];
t_sng_psf = padarray(single_psf,[round((spacing-single_f_size)/2),...
    round((spacing-single_f_size)/2)]);
mla_psf = repmat(t_sng_psf,copies);
figure,surf(mla_psf),shading flat,title('mla psf');


%% conv obj with mla psf
mla_conv = myconv(obj,mla_psf,'same');
if size(obj,1)>256
    [rr,cc] = size(obj);
    rr = floor(rr./2);
    cc = floor(cc./2);
    mla_conv = mla_conv(rr-255:rr+256,cc-255:cc+256);
end
% mla_conv = mla_conv(128:128+255,128:128+255);
figure,imagesc(mla_conv),axis image off,colormap gray;colorbar,title('conv with mla psf');


% %% try deconv
% F = @(x) fftshift(fft2(ifftshift(x)));
% Ft = @(x) fftshift(ifft2(ifftshift(x)));
% f_mea = F(mla_conv);
% % figure,imagesc(log(abs(f_mea))),axis image off, colormap gray;colorbar;
% dscrp = size(mla_conv,1) - size(mla_psf,1);
% if mod(dscrp,2) == 0
%     mla_psf = padarray(mla_psf,[dscrp/2,dscrp/2]);
% else
%     mla_psf = padarray(mla_psf,[(dscrp+1)/2,(dscrp+1)/2]);
%     mla_psf = mla_psf(1:end-1,1:end-1);
% end 
% 
% 
% 
% f_psf = F(mla_psf);
% figure,imagesc(log(abs(f_psf))),axis image off, colormap gray;colorbar,title('TF');
% f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+0);
% f_res2 = (f_mea).*conj(f_psf)./(abs(f_psf).^2+5);
% f_res3 = (f_mea).*conj(f_psf)./(abs(f_psf).^2+1e11);
% 
% % figure,imagesc(log(abs(f_res))),axis image off, colormap gray;colorbar;
% result = abs(Ft(f_res));
% figure,imagesc(result),axis image off, colormap gray;colorbar,title('deconv result para 0');
% 
% result2 = abs(Ft(f_res2));
% figure,imagesc(result2),axis image off, colormap gray;colorbar,title('deconv result para 5');
% 
% result3 = abs(Ft(f_res3));
% figure,imagesc(result3),axis image off, colormap gray;colorbar,title('deconv result para 1e11');


%% add gausssian noise


F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

noise = 0.01*randn(size(mla_conv));
% noise = 0;
ps = sum((mla_conv(:)-mean(mla_conv(:))).^2);
pn = sum(noise(:).^2);
snr = 10*log10(ps/pn);
mla_conv_noise = mla_conv+noise;
figure,imagesc(mla_conv_noise);axis image off, colormap gray;colorbar;

mla_psf = padarray(mla_psf,size(mla_psf));
mla_conv = padarray(mla_conv_noise,size(mla_conv_noise));
dscrp = size(mla_conv,1) - size(mla_psf,1);
if mod(dscrp,2) == 0
    mla_psf = padarray(mla_psf,[dscrp/2,0]);
else
    mla_psf = padarray(mla_psf,[(dscrp+1)/2,0]);
    mla_psf = mla_psf(1:end-1,:);
end
dscrp = size(mla_conv,2) - size(mla_psf,2);
if mod(dscrp,2) == 0
    mla_psf = padarray(mla_psf,[0,dscrp/2]);
else
    mla_psf = padarray(mla_psf,[0,(dscrp+1)/2]);
    mla_psf = mla_psf(:,1:end-1);
end


f_mea = F(mla_conv);
f_psf = F(mla_psf);
para = 0.03*max(abs(f_psf(:)));
% para =0;
% f_mea_n = F(mla_conv_noise);
f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
% f_res = f_mea./(eps+f_psf);
% f_res2 = (f_mea_n).*conj(f_psf)./(abs(f_psf).^2+5);
% f_res3 = (f_mea_n).*conj(f_psf)./(abs(f_psf).^2+1e11);
% 
% figure,imagesc(log(abs(f_res))),axis image off, colormap gray;colorbar;
result = real(Ft(f_res));
result = mytruncate(result,0,255);
result = result(size(result,1)/3+1:2*size(result,1)/3,size(result,2)/3+1:2*size(result,2)/3);
figure,imagesc(result),axis image off, colormap gray;colorbar;
% 
% result2 = abs(Ft(f_res2));
% figure,imagesc(result2),axis image off, colormap gray;colorbar,title('deconv result para 5');
% 
% result3 = abs(Ft(f_res3));
% figure,imagesc(result3),axis image off, colormap gray;colorbar,title('deconv result para 1e11');
% 
% for para = 0:0.1:0.5
%     ftmp = (f_mea_n).*conj(f_psf)./(abs(f_psf).^2+para);
%     result = abs(Ft(ftmp));
%     figure,imagesc(result),axis image off, colormap gray;colorbar,title(para);
% end
%     
