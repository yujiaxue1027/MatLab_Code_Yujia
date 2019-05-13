%%script on simulation of BPM (beam propagation method)
%%detailed explanation of the simulation parameters should be the same as
%Ulugbek's paper, unit: all in nanometer

%%construct the object layer
%%two layer first is a rectangle window, the second is lena
img_row = 512;
img_col = 512;
dz0 = 50;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.633;

% o_slices_r(:,:,1) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,2) = double(imread('lena.png'))./255;
% o_slices_r(:,:,3) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,4) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,5) = double(imread('lena.png'))./255;
% o_slices_r(:,:,6) = double(imread('mask2.png'))./255;
% o_slices_t = 1-o_slices_r;

scene = 4;

switch scene
    case 1
        %% first layer random reflection; second layer mirror-like
        Nslices = 2;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices]);
        xigma = 0.1;
        o_slices_t(:,:,1) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
        r2 = 0.95; %reflectance on the second layer
        o_slices_t(:,:,2) = (1-r2)*ones(img_row,img_col);
    case 2
        %%first layer resolution target, second layer random reflection
        Nslices = 2;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices]);
        o_slices_t(:,:,1) = 1 - imread('resolution512.png')./255;
        xigma = 0.1;
        o_slices_t(:,:,2) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
    case 3
        %%first 2 layers are resolution target, third layer random reflection
        Nslices = 3;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices]);
        o_slices_t(:,:,1) = 1 - imread('resolution512.png')./255;
        o_slices_t(:,:,2) = 1 - imread('resolution5122.png')./255;
        xigma = 0.1;
        o_slices_t(:,:,3) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
    case 4
        %%first 2 layers are resolution target, third layer mirror
        Nslices = 3;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices]);
        o_slices_t(:,:,1) = 1 - imread('resolution512.png')./255;
        o_slices_t(:,:,2) = 1 - imread('resolution5122.png')./255;
        r2 = 0.95;
        o_slices_t(:,:,3) = (1-r2)*ones(img_row,img_col);
end





%%assume r + t = 1
o_slices_r = 1-o_slices_t;

%%some setup parameters
pitch = 1;
FOV = (pitch*img_row);
du = 1/FOV;
Umax = du * img_row/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);

k2 = pi*lambda*(u.^2+v.^2);
i0 = exp(1i*2*pi*(x*0+y*0));

[ phi, psi ] = Fwd_Prop_MultiSlice_v2( i0, o_slices_t, k2, dz);
%%plot incident and exit field at each slice
for i = 1:Nslices
figure,subplot(1,2,1),imagesc(abs(phi(:,:,i))),title(['incident field of ',num2str(i),' layer']);
subplot(1,2,2),imagesc(abs(psi(:,:,i))),title(['exit field of ',num2str(i),' layer']);
colormap gray, truesize
end

%%compute back scattered field by single backscattering approximation and
%%beam propagation
backfield = zeros(img_row,img_col,Nslices); %%reflected field at each layer
%% reflected field at the first layer is first layer's incident field times reflectance of the first layer
backfield(:,:,1) = phi(:,:,1).*o_slices_r(:,:,1);
Prop = @(f,h) Ft(F(f).*h);
%%for each layer, first times its reflectance, then apply beam propagation
for i = 2:Nslices
    %%reflected here
    tmp = phi(:,:,i).*o_slices_r(:,:,i);
    %%propagate for a distance between ith and i-1th layer
    H = exp(1i*k2*dz(i-1));
    ti0 = Prop(tmp,H);
    toslices = o_slices_t(:,:,i-1:-1:1);
    tdz = dz(i-2:-1:1);
    [ tphi, tpsi ] = Fwd_Prop_MultiSlice_v2( ti0, toslices, k2, tdz);
    tback = tpsi(:,:,end);
    backfield(:,:,i) = tback;
end
%%the total back scattered field should be the sum of all backscattered
%%field from each layer
back_scattered_field = sum(backfield,3);

%%define a final propagation distance from the first layer to the detection
%%plane
zn = 50;
Hn = exp(1i*k2*zn);
nfield = Prop(back_scattered_field,Hn);
I = abs(nfield.^2);
figure,imagesc(I),title('intensity measurement of back scattered field'),colormap gray, truesize, axis image off;
    







