function [output,history] = CM2_2D_FISTA_Solver(y,psf,para)
% Yujia's implementation on solving the inverse problem y = CHx. 
% which is the 2D model of CM2, assuming LSI and cropping. 
% The inverse problem is solved using FISTA algorithm.
%% load prepared data
[rows,cols] = size(psf);
l1_thres = para.l1_thres;
num_iters = para.num_iters;
alpha = para.alpha;
color = para.color;
clip_min = para.clip_min;
clip_max = para.clip_max;


%% define operator
C = @(x) x(rows/2+1:rows/2+rows, cols/2+1:cols/2+cols);
CT  = @(x) padarray(x,[rows/2,cols/2]);
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

%% initialize FISTA
tk = 1;
x0 = zeros(2*rows,2*cols);
xkm1 = x0;
yk = x0;
xk = x0;
output = zeros(rows*2,cols*2);
history = zeros(rows*2,cols*2,num_iters);

f1 = figure('rend','painters','pos',[50 50 1500 900]);
%% main body of FISTA
for i = 1:num_iters
    residual = C(H(xk,psf)) - y;
    grad = HT(CT(residual),psf);
    xk = wthresh(yk-alpha.*grad,'s', l1_thres);
    xk = clip(xk,clip_min,clip_max);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    % update
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    output = xk;
    history(:,:,i) = xk;
    % plot
    figure(f1);
    imagesc(output),colormap(color),colorbar,title(i);
end
end

%% other functions
function OUTPUT = H(INPUT, PSF)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
                1+size(x,2)/4:size(x,2)/4*3);
conv2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*F(pad2d(psf)))));
[ROWS,COLS] = size(PSF);
PSF = padarray(PSF,[ROWS/2,COLS/2]);
OUTPUT = conv2d(INPUT, PSF);
end

function OUTPUT = HT(INPUT, PSF)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
                1+size(x,2)/4:size(x,2)/4*3);
convT2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*conj(F(pad2d(psf))))));
[ROWS,COLS] = size(PSF);
PSF = padarray(PSF,[ROWS/2,COLS/2]);
OUTPUT = convT2d(INPUT, PSF);
end



