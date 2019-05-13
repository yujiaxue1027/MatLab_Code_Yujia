function [output,history] = CM2_3D_ADMM_Solver(y,psfs,para)
%% Yujia's ADMM implementation to solve crop and conv model
% y = CHx, y: measure; C:cropping operator; H:conv operator; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|CHx-y|^2+tau1*|Psix|_1+nonnega_indicator{x}+tau2*|x|_1

%% load prepared data
[rows,cols,depths] = size(psfs);
mu1 = para.mu1;
mu2 = para.mu2;
mu3 = para.mu3;
tau_over_mu2 = para.l1_thres;
tau = tau_over_mu2*mu2;
maxiter = para.maxiter;
color = para.color;
clip_min = para.clip_min;
clip_max = para.clip_max;
%% ADMM parameters
history = zeros(2*rows,2*cols,depths, maxiter);

%% define operators
C = @(x) x(rows/2+1:rows/2+rows, cols/2+1:cols/2+cols);
CT  = @(x) padarray(x,[rows/2,cols/2]);
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
S = @(x) sum(x, 3);
ST = @(x) repmat(x,[1,1,depths]);

%% ADMM algorithm
STCTy = ST(CT(y));
xt = zeros(2*rows, 2*cols, depths);
vt = xt;
zt = xt;
wt = xt;
gamma1 = vt;
gamma2 = zt;
gamma3 = wt;

iter = 0;
f1 = figure('rend','painters','pos',[50 50 1500 900]);
f2 = figure('rend','painters','pos',[50 50 1500 900]);

while iter < maxiter
   iter = iter + 1;
   vtp = v_update_inv_operator(STCTy + mu1.*H(xt,psfs) + gamma1, mu1);
   ztp = wthresh(gamma2/mu2 + xt,'s',tau_over_mu2);
   wtp = clip(gamma3/mu3 + xt, clip_min,clip_max);
   
   
   
   tmp_part1 = mu1*HT(vtp-gamma1/mu1, psfs);
   tmp_part2 = mu2*(ztp-gamma2/mu2);
   tmp_part3 = mu3*(wtp - gamma3/mu3);
   xtp_numerator = tmp_part1 + tmp_part2 + tmp_part3;
   xtp = x_update_inv_operator(xtp_numerator, psfs, mu1,mu2,mu3);
   
   % update duals
   gamma1 = gamma1 + mu1*(H(xtp, psfs)-vtp);
   gamma2 = gamma2 + mu2*(xtp-ztp);
   gamma3 = gamma3 + mu3*(xtp-wtp);
  
   % update xt
   xt = xtp;
   output = xt;
   history(:,:,:,iter) = xt;
   
   % display and evaluate
   figure(f1);
   subplot(1,2,1),imagesc(S(xt)),colormap(color);axis image;colorbar;title(iter);

   residual = abs(C(S(H(xt,psfs)))-y);
   dterm = 0.5*sum(residual(:).^2);
   l1term = tau*sum(abs(xt(:)));
   cost = dterm+l1term;
   subplot(1,2,2),plot(iter,log10(cost),'bo'),grid on,hold on;...
            plot(iter,log10(dterm),'ro'),hold on;...
            plot(iter,log10(l1term),'mo'),hold on;...
            title('log axis: blue: cost; red: data fidelity; purple: l1');
   drawnow; 
   figure(f2);
   subplot(2,2,1),imagesc(S(tmp_part1)),colormap(color);axis image;colorbar;title(['data: ',num2str(iter)]);
   subplot(2,2,2),imagesc(S(tmp_part2)),colormap(color);axis image;colorbar;title(['l1: ',num2str(iter)]);
   subplot(2,2,3),imagesc(S(tmp_part3)),colormap(color);axis image;colorbar;title(['NN: ',num2str(iter)]);
   drawnow;
end
end


function OUTPUT = H(INPUT, PSFS)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
                1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
[ROWS,COLS,DEPTHS] = size(PSFS);
PSFS = padarray(PSFS,[ROWS/2,COLS/2,0]);
OUTPUT = zeros(ROWS*2, COLS*2, DEPTHS);
for ii = 1:DEPTHS
    OUTPUT(:,:,ii) = conv2d(INPUT(:,:,ii), PSFS(:,:,ii));
end
end

function OUTPUT = HT(INPUT, PSFS)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
                1+size(x,2)/4:size(x,2)/4*3);
conv2dT = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*conj(F(pad2d(psf))))));
[ROWS,COLS,DEPTHS] = size(PSFS);
PSFS = padarray(PSFS,[ROWS/2,COLS/2,0]);
OUTPUT = zeros(ROWS*2, COLS*2, DEPTHS);
for ii = 1:DEPTHS
    OUTPUT(:,:,ii) = conv2dT(INPUT(:,:,ii), PSFS(:,:,ii));
end
end

function OUTPUT = v_update_inv_operator(INPUT, MU)
part1 = INPUT./MU;
[ROWS,COLS,DEPTHS] = size(INPUT);
mask_CTC = ones(ROWS/2,COLS/2);
mask_CTC = padarray(mask_CTC,[ROWS/4,COLS/4]);
mask_CTC_3d = repmat(mask_CTC,[1,1,DEPTHS]);
part2 = INPUT.*mask_CTC_3d;
part2 = sum(part2,3);
denominator = ones(ROWS,COLS) + mask_CTC.*DEPTHS./MU;
part2 = part2./denominator;
part2 = repmat(part2,[1,1,DEPTHS])./(MU^2);
OUTPUT = part1 - part2;
end

function OUTPUT = x_update_inv_operator(INPUT, PSFS, MU1,MU2,MU3)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
[ROWS,COLS,DEPTHS] = size(PSFS);
PSFS = padarray(PSFS,[ROWS/2,COLS/2,0]);
OUTPUT = zeros(ROWS*2, COLS*2, DEPTHS);
for ii = 1:DEPTHS
    layer = INPUT(:,:,ii);
    tf = F(PSFS(:,:,ii));
    OUTPUT(:,:,ii) = real(Ft(F(layer)./(MU1.*tf.*conj(tf)+MU2+MU3)));
end
end