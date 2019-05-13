%% generate psf of single lens, every thing in um
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.633;
object_pitch = 1;
sensor_pitch = 2.2;
single_lens_dia = 1000;
f = 3300; 
simu_area = 4747;   %% i want the spatial unit of the psf to be divisable by 2.2 um(which is sensor pitch)
simu_dim = round(simu_area/object_pitch);
aperture = zeros(simu_dim);
x = object_pitch.*[-(simu_dim-1)/2:1:(simu_dim-1)/2];
y = object_pitch.*[-(simu_dim-1)/2:1:(simu_dim-1)/2];
[x,y] = meshgrid(x,y);
aperture(sqrt(x.^2+y.^2)<=single_lens_dia/2) = 1;
f_aperture = F(aperture);
psf_up = abs(f_aperture);
fov = simu_area;
du = 1/fov;
u = du.*[-(simu_dim-1)/2:1:(simu_dim-1)/2];
v = du.*[-(simu_dim-1)/2:1:(simu_dim-1)/2];
[u,v] = meshgrid(u,v);
xprime = u.*lambda.*f;
yprime = v.*lambda.*f;
%%unit of f_aperture is 0.44um which is one fifth of sensor pitch, do 5*5
%%pixel binning here
pos_center = (simu_dim+1)/2;
psf_dim = 31;
binning = 5;
tmp_dim = binning*psf_dim;
idx_start = pos_center - (binning-1)/2 - binning*((psf_dim-1)/2);
tmp_area = psf_up(idx_start:idx_start+psf_dim*binning-1,idx_start:idx_start+psf_dim*binning-1);
psf = zeros(psf_dim,psf_dim);
for i = 1:psf_dim
    for j = 1:psf_dim
        psf(i,j) = mean(mean(tmp_area(binning*(i-1)+1:i*binning,binning*(j-1)+1:j*binning)));
    end
end
psf = psf./sum(psf(:));
psf(psf<=0.001) = 0;

%% mla psf lenslet spacing 1mm
lenslet_spacing = 1000; %%1mm =1000um
spacing = round(lenslet_spacing/sensor_pitch);
copies = [3,3];
t_sng_psf = padarray(psf,[round((spacing-psf_dim)/2),...
    round((spacing-psf_dim)/2)]);
mla_psf = repmat(t_sng_psf,copies);
img_dim = 2000;
mla_psf = padarray(mla_psf,[round((img_dim-size(mla_psf,1))/2),...
    round((img_dim-size(mla_psf,1))/2)]);
if size(mla_psf,1) ~= img_dim
    mla_psf = mla_psf(1:end-1,1:end-1);
end
figure,imagesc(mla_psf),colormap jet;colorbar,title('mla psf');



%% load brain image and generate mla measurement
obj = rgb2gray(imread('brain2.jpg'));
scale = 3; %%5 3 if the conv output do not exceed sensor area
obj = double(imresize((obj),scale));
obj = obj./max(obj(:));
obj = padarray(obj,[round((img_dim-size(obj,1))/2),...
    round((img_dim-size(obj,1))/2)]);
if size(obj,1) ~= img_dim
    obj = obj(1:end-1,1:end-1);
end
figure,imagesc(obj),colormap jet;colorbar;title('obj');
measure = myconv(obj,mla_psf,'same');
figure,imagesc(measure),colormap jet;colorbar;title('conv without noise')

%% adding noise
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
noise = 0.01*randn(size(measure));
ps = sum((measure(:)-mean(measure(:))).^2);
pn = sum(noise(:).^2);
snr = 10*log10(ps/(pn+eps));
measure_noise = measure+noise;
figure,imagesc(measure_noise);axis image off, colormap jet;colorbar;title('conv with noise');

%% deconvolution
reg = 0.01;  %%decrease reg if you see copies, increase if if cannot tell the obj
zero_padding_flag = 1;
if zero_padding_flag == 1
    mla_psf_n = padarray(mla_psf,size(mla_psf));
    measure_noise_n = padarray(measure_noise,size(measure_noise));
    f_mea = F(measure_noise_n);
    f_psf = F(mla_psf_n);
    para = reg*max(abs(f_psf(:)));
    f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
    result = real(Ft(f_res));
    result = result(size(result,1)/3+1:2*size(result,1)/3,size(result,2)/3+1:2*size(result,2)/3);
else
    f_mea = F(measure_noise);
    f_psf = F(mla_psf);
    para = reg*max(abs(f_psf(:)));
    f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
    result = real(Ft(f_res));
end
result = mytruncate(result,0,1);
figure,imagesc(result),axis image off, colormap jet;colorbar;title('result');





%% assume large object and small sensor, try deconv to retrieve large object
obj = rgb2gray(imread('brain2.jpg'));
scale = 9; %% 9*350>2000
obj = double(imresize((obj),scale));
obj = obj./max(obj(:));
figure,imagesc(obj),colormap jet;colorbar;title('large obj');
measure = myconv(obj,mla_psf,'same');
tmp_size = size(measure,1);
measure = measure(round(tmp_size/2)-img_dim/2:round(tmp_size/2)+img_dim/2-1,...
    round(tmp_size/2)-img_dim/2:round(tmp_size/2)+img_dim/2-1);
figure,imagesc(measure),colormap jet;colorbar;title('conv without noise');
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
noise = 0.01*randn(size(measure));
ps = sum((measure(:)-mean(measure(:))).^2);
pn = sum(noise(:).^2);
snr = 10*log10(ps/(pn+eps));
measure_noise = measure+noise;
figure,imagesc(measure_noise);axis image off, colormap jet;colorbar;title('conv with noise');
reg = 0.07;  %%decrease reg if you see copies, increase if if cannot tell the obj
mla_psf_n = padarray(mla_psf,size(mla_psf));
measure_noise_n = padarray(measure_noise,size(measure_noise));
f_mea = F(measure_noise_n);
f_psf = F(mla_psf_n);
para = reg*max(abs(f_psf(:)));
f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
result = real(Ft(f_res));
result = mytruncate(result,0,1);
result = result(round(size(result,1)/2)-size(obj,1)/2:round(size(result,1)/2)+size(obj,1)/2-1,...
    round(size(result,1)/2)-size(obj,1)/2:round(size(result,1)/2)+size(obj,1)/2-1);
figure,imagesc(result),axis image off, colormap jet;colorbar;title('result');
    
% 
% %% symmetric padding test
% mla_psf_n = padarray(mla_psf,size(mla_psf));
% tmat1 = flip(measure_noise,1);
% tmat2 = flip(measure_noise,2);
% tmat3 = zeros(size(measure_noise));
% measure_noise_n = [tmat3,tmat1,tmat3;tmat2,measure_noise,tmat2;tmat3,tmat1,tmat3];
% figure,imagesc(measure_noise_n),axis image off, colormap gray;colorbar;title('padded');
% f_mea = F(measure_noise_n);
% f_psf = F(mla_psf_n);
% para = reg*max(abs(f_psf(:)));
% f_res = (f_mea).*conj(f_psf)./(abs(f_psf).^2+para^2);
% result = real(Ft(f_res));
% result = mytruncate(result,0,1);
% result = result(round(size(result,1)/2)-size(obj,1)/2:round(size(result,1)/2)+size(obj,1)/2-1,...
%     round(size(result,1)/2)-size(obj,1)/2:round(size(result,1)/2)+size(obj,1)/2-1);
% figure,imagesc(result),axis image off, colormap gray;colorbar;title('result');
%     


%% implement fista to solve the crop-conv inverse problem
%%initalize parameters for fista
factor = 1;
nobj = imresize(obj,factor);
nmeasure_noise = imresize(measure_noise,factor);
nmla_psf = imresize(mla_psf,factor);
tk = 1;
x0 = zeros(size(nobj));
xkm1 = x0;
yk = x0;
xk = x0;
maxiter =100;%100
alpha = 0.003;%0.005
thres = 0.003;%0.0585
psf2 = flip(flip(nmla_psf,1),2);
psf2 = padarray(psf2,[575,575]);
result = zeros(3150,3150,maxiter);
mselist = [];
decrease_flag = 0;
back_length = 8;
decrease_thres = 1;
decrease_rate = 0.5;
%%main iteration here:
for i = 1:maxiter
    residual = forward_cac(xk,nmla_psf,2000)-nmeasure_noise;
    nresidual = padarray(residual,[575,575]);% (3150 - 2000)/2= 575;
    grad = conv2(nresidual,psf2,'same');
    xk = myproximal(yk-alpha.*grad,thres);
%     xk = mytruncate(xk,0,1);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    figure(51),imagesc(xk),axis image off;colormap jet;title(i);
    result(:,:,i) = xk;
    mse = norm(residual(:));
    mselist = [mselist,mse];
    figure(52),plot(i,mse,'bo'),grid on,hold on;title(['data term, alpha: ',num2str(alpha),' thres: ',num2str(thres)])
    if decrease_flag
    if i>=40
        recentmse = mselist(i-back_length:i);
        p = polyfit([1:back_length+1],recentmse,1);
        if abs(p(1))<decrease_thres
            alpha = alpha*decrease_rate;
            decrease_thres = decrease_thres*decrease_rate;
            disp(['at ',num2str(i),'th iteration, decrease alpha, alpha: ',num2str(alpha)]);
        end
    end
    end
    pause(1);
end










