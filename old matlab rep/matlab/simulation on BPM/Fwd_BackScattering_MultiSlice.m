function [totalfield, backfield, transmission] = Fwd_BackScattering_MultiSlice( s1, o_slices_t,o_slices_r, k2, dz )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
[img_row,img_col,Nslices] = size(o_slices_t);
[ phi, psi ] = Fwd_Prop_MultiSlice_v2( s1, o_slices_t, k2, dz);
transmission = psi(:,:,end);
%%compute back scattered field by single backscattering approximation and
%%beam propagation
backfield = zeros(img_row,img_col,Nslices); %%reflected field at each layer
%% reflected field at the first layer is first layer's incident field times reflectance of the first layer
backfield(:,:,1) = phi(:,:,1).*o_slices_r(:,:,1);
%%for each layer, first times its reflectance, then apply beam propagation
Prop = @(f,h) Ft(F(f).*h);

for i = 2:Nslices
    %%reflected here
    tmp = phi(:,:,i).*o_slices_r(:,:,i);
    %%propagate for a distance between ith and i-1th layer
    H = exp(1i*k2*dz(i-1));
    ti0 = Prop(tmp,H);
    toslices = o_slices_t(:,:,i-1:-1:1);
    tdz = dz(i-2:-1:1);
    [ ~, tpsi ] = Fwd_Prop_MultiSlice_v2( ti0, toslices, k2, tdz);
    tback = tpsi(:,:,end);
    backfield(:,:,i) = tback;
end
%%the total back scattered field should be the sum of all backscattered
%%field from each layer
totalfield = sum(backfield,3);

end

