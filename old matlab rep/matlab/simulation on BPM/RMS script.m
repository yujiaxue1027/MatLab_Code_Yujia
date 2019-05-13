%%script on simulation of BPM (beam propagation method)
%%detailed explanation of the simulation parameters should be the same as
%Ulugbek's paper, unit: all in nanometer

%%construct the object layer
%%two layer first is a rectangle window, the second is lena

%%new on 1005 test generate forward model 
img_row = 256;
img_col = 256;
dz0 = 0.144;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.633;
pitch = 0.144;
scene = 5;

switch scene
    case 1
        %% first layer random reflection; second layer mirror-like
        Nslices = 2;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices-1]);
        xigma = 0.1;
        o_slices_t(:,:,1) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
        r2 = 0.95; %reflectance on the second layer
        o_slices_t(:,:,2) = (1-r2)*ones(img_row,img_col);
    case 2
        %%first layer resolution target, second layer random reflection
        Nslices = 2;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices-1]);
        o_slices_t(:,:,1) = 1 - imresize(double(imread('resolution512.png')),[img_row,img_col])./255;
        xigma = 0.1;
        o_slices_t(:,:,2) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
    case 3
        %%first 2 layers are resolution target, third layer random reflection
        Nslices = 3;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices-1]);
        o_slices_t(:,:,1) = 1 - imread('resolution512.png')./255;
        o_slices_t(:,:,2) = 1 - imread('resolution5122.png')./255;
        xigma = 0.1;
        o_slices_t(:,:,3) = mytruncate((xigma*randn(img_row,img_col)+0.5),0,1);
    case 4
        %%first 2 layers are resolution target, third layer mirror
        Nslices = 3;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices-1]);
        o_slices_t(:,:,1) = 1 - imread('resolution512.png')./255;
        o_slices_t(:,:,2) = 1 - imread('resolution5122.png')./255;
        r2 = 0.95;
        o_slices_t(:,:,3) = (1-r2)*ones(img_row,img_col);
   case 5
        %%256*256*256 layers a bead in a cube space
        Nslices = 256;
        n0 = 1.518;
        n1 = 1.548;
        o_slices_r = zeros(img_row,img_col,Nslices);
        o_slices_t = zeros(img_row,img_col,Nslices);
        n = n0*ones(img_row,img_col,Nslices);
        dz = repmat(dz0,[1,Nslices-1]);
        tcenter = [(img_col+1)/2*pitch,(1+img_row)/2*pitch,(Nslices+1)/2*dz0];
        radius = 3.5;
        for j = 1:img_col
            for i = 1:img_row
                for k = 1:Nslices
                    tpos = [j*pitch,i*pitch,k*dz0];
                    tdist = norm(tpos-tcenter);
                    if tdist <= radius
                        n(i,j,k) = n1;
                    end
                end
            end
        end
        %% my derivation on oct 31 night
        beta = 2*pi/lambda*(n1-n0)*dz0;
        eibeta = exp(1i*beta);
        eimbeta = exp(-1i*beta);
%         o_slices_t = eibeta*(1+eibeta*(eibeta-eimbeta)/(4*n0*n0).*(n-n0).*(n-n0));
        o_slices_t = exp(1i*2*pi/lambda.*(n-n0)*dz0);
        o_slices_r = eibeta*(eibeta-eimbeta)/2/n0.*(n-n0);
        
        %% ulugbek's object function
%         o_slices_t_ulug = exp(1i*2*pi/lambda.*(n-n0)*dz0);
        
%         %% 1
%         o_slices_t = 2*n0./(n+n0);
%         o_slices_r = 2*n0./(n0+n)-1;
%         %% 2
%         o_slices_t = 4*n0*n./(n+n0).^2;
%         o_slices_r = 2*n0./(n0+n)-1;
%         %% 3
%         o_slices_t = 4*n0*n./(n+n0).^2;
%         o_slices_r = 2*n0./(n0+n)-1+2*n0./(n+n0)*(2*n./(n+n0)-1);

end
% o_slices_r = o_slices_t-1;




%%assume r + t = 1


%%some setup parameters

FOV = (pitch*img_row);
du = 1/FOV;
Umax = du * img_row/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);

k2 = n0*pi*lambda*(u.^2+v.^2);
%%define i0 as the field that is dz0 ahead of the first layer, which can
%%be thought as 0th slice
thetax = 0;
thetay = 0;
i0 = exp(1i*2*pi*(x*thetax+y*thetay));
zinitial = dz0;
Prop = @(f,h) Ft(F(f).*h);
Hinitial = exp(1i*k2*zinitial);
i1 = Prop(i0,Hinitial);
[totalfield, backfield, transmission] = Fwd_BackScattering_MultiSlice( i1, o_slices_t,o_slices_r, k2, dz );
[totalfield_u, backfield_u, transmission_u] = Fwd_BackScattering_MultiSlice( i1, o_slices_t_ulug,o_slices_r, k2, dz );

% zfinal = dz0;
% Hfinal = exp(1i*k2*zfinal);
% finalfield = Prop(totalfield,Hfinal);
figure(1),imagesc(abs(totalfield)),title('intensity measurement of back scattered field'),colormap gray, truesize, axis image off,colorbar;
figure(2),imagesc(abs(transmission)),title('intensity measurement of transmission field'),colormap gray, truesize, axis image off,colorbar;
saveas(figure(1),'bpm_reflc3.png');
saveas(figure(2),'bpm_trans3.png');








