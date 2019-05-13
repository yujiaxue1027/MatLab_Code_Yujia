%%simple bpm simulation

img_row = 512;
img_col = 512;
Nslices = 10;
% o_slices_r = zeros(img_row,img_col,Nslices);
% o_slices_t = zeros(img_row,img_col,Nslices);
% o_slices_r = ones(img_row,img_col,Nslices);
o_slices_t = ones(img_row,img_col,Nslices);
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
%define dz , dz has a dimension of 1 by Nslices-1, uniform interval of dz0
dz0 = 50;
dz = repmat(dz0,[1,Nslices-1]);
lambda = 0.643;

%%defnie first layer object
firstlayer = zeros(img_row,img_col);
center = [round((img_row+1)/2),round((img_col+1)/2)];
for i = 1:img_row
    for j = 1:img_col
        dia = 21;
        pos = [i,j];
        if norm(pos-center) <= dia/2
            firstlayer(i,j)  = 1;
        end
    end
end


o_slices_t(:,:,1) = firstlayer;

% o_slices_r(:,:,1) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,2) = double(imread('lena.png'))./255;
% o_slices_r(:,:,3) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,4) = double(imread('mask2.png'))./255;
% o_slices_r(:,:,5) = double(imread('lena.png'))./255;
% o_slices_r(:,:,6) = double(imread('mask2.png'))./255;
% o_slices_t = 1-o_slices_r;
% o_slices_r = 1-o_slices_t;

pitch = 1;
FOV = (pitch*img_row);
du = 1/FOV;
Umax = du * img_row/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);

k2 = pi*lambda*(u.^2+v.^2);
i0 = exp(1i*2*pi*(x*0+y*0));

[ phi, psi ] = Fwd_Prop_MultiSlice_v2( i0, o_slices_t, k2, dz);
for i = 1:Nslices
figure,subplot(1,2,1),imagesc(abs(phi(:,:,i))),title(['incident field of ',num2str(i),' layer']);
subplot(1,2,2),imagesc(abs(psi(:,:,i))),title(['exit field of ',num2str(i),' layer']);
colormap gray, truesize
end
% 
% backfield = zeros(img_row,img_col,Nslices);
% backfield(:,:,1) = phi(:,:,1).*o_slices_r(:,:,1);
% Prop = @(f,h) Ft(F(f).*h);
% for i = 2:Nslices
%     tmp = phi(:,:,i).*o_slices_r(:,:,i);
%     H = exp(1i*k2*dz(i-1));
%     % propagate from neiboring slices
%     ti0 = Prop(tmp,H);
%     toslices = o_slices_t(:,:,i-1:-1:1);
%     tdz = dz(i-2:-1:1);
%     [ tphi, tpsi ] = Fwd_Prop_MultiSlice_v2( ti0, toslices, k2, tdz);
%     tback = tpsi(:,:,end);
%     backfield(:,:,i) = tback;
% end
% back_scattered_field = sum(backfield,3);
% 
% zn = 50;
% Hn = exp(1i*k2*zn);
% nfield = Prop(back_scattered_field,Hn);
I = abs(psi(:,:,end).^2);
figure,imagesc(I),title('intensity measurement of transmission field'),colormap gray, truesize, axis image off;
figure,imagesc(abs(psi(:,:,end))),title('amplitude of transmission field'),colormap gray, truesize, axis image off;


theta = asin(1.22*lambda/dia);
L = dz0*(Nslices-1);
r_theo_um = tan(theta)*L;
disp(['theoretical radius: ',num2str(r_theo_um),' um']);
r_theo_pixel = r_theo_um/pitch;
disp(['which is: ',num2str(r_theo_pixel),' pixels']);
plot(abs(squeeze(psi(257,:,end))))







