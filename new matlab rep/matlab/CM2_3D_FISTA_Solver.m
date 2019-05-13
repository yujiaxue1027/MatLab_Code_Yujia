function [output,history] = CM2_3D_FISTA_Solver(y,psfs,para)
% Yujia's implementation on solving the inverse problem y = CSHx. 
% which is the 3D model of CM2, assuming each layer is LSI and cropping
% at the end. The inverse problem is solved using FISTA algorithm.
%% load prepared data
[rows,cols,depths] = size(psfs);
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
S = @(x) sum(x, 3);
ST = @(x) repmat(x,[1,1,depths]);

%% initialize FISTA
tk = 1;
x0 = zeros(2*rows,2*cols,depth);
xkm1 = x0;
yk = x0;
xk = x0;
output = zeros(rows*2,cols*2,depth);
history = zeros(rows*2,cols*2,depth,num_iters);

f1 = figure('rend','painters','pos',[50 50 1500 900]);
f2 = figure('rend','painters','pos',[50 50 1500 900]);
num_subfig = ceil(sqrt(depth));
%% main body of FISTA
for i = 1:num_iters
    residual = C(S(H(xk,psfs))) - y;
    grad = HT(ST(CT(residual)),psfs);
    xk = wthresh(yk-alpha.*grad,'s', l1_thres);
    xk = clip(xk,clip_min,clip_max);
    tkp1 = (1+sqrt(1+4*tk*tk))/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    % update
    yk = ykp1;
    tk = tkp1;
    xkm1 = xk;
    output = xk;
    history(:,:,:,i) = xk;
    % plot
    figure(f1);
    imagesc(S(output)),colormap(color),colorbar,title(i);
    figure(f2);
    for j = 1:depths
        subplot(num_subfig,num_subfig,j);
        imagesc(output(:,:,j)),colormap(color),colorbar;
    end 
end
end

%% other functions
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
for i = 1:DEPTHS
    OUTPUT(:,:,i) = conv2d(INPUT(:,:,i), PSFS(:,:,i));
end
end

function OUTPUT = HT(INPUT, PSFS)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
                1+size(x,2)/4:size(x,2)/4*3);
convT2d = @(obj,psf) crop2d(real(Ft(F(pad2d(obj)).*conj(F(pad2d(psf))))));
[ROWS,COLS,DEPTHS] = size(PSFS);
PSFS = padarray(PSFS,[ROWS/2,COLS/2,0]);
OUTPUT = zeros(ROWS*2, COLS*2, DEPTHS);
for i = 1:DEPTHS
    OUTPUT(:,:,i) = convT2d(INPUT(:,:,i), PSFS(:,:,i));
end
end



