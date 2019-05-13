function output = conv3d(obj,psf)
    F3D = @(x) fftshift(fftn(ifftshift(x)));
    Ft3D = @(x) fftshift(ifftn(ifftshift(x)));
    output = crop3d(real(Ft3D(F3D(pad3d(obj)).*F3D(pad3d(psf)))));
end