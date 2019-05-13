%% script on mla psf analysis when the pixels downsample the psf
F = @(x) fftshift(fft(ifftshift(x)));
Ft = @(x) fftshift(ifft(ifftshift(x)));
psf_single = fspecial('gaussian',255,20);
psf_single = psf_single(:,128);
figure,plot(psf_single),title('single psf');
psf1 = zeros(3001,1);
psf1(1501-127:1501+127) = psf_single;
psf2 = zeros(3001,1);
psf2(501-127:501+127) = psf_single;
psf2(1501-127:1501+127) = psf_single;
psf2(2501-127:2501+127) = psf_single;
figure,plot(psf1),title('psf1');
figure,plot(psf2),title('psf2');
psf1_ds = binning_1d(psf1,32);
psf2_ds = binning_1d(psf2,32);
psf1_ds = psf1_ds./max(psf1_ds);
psf2_ds = psf2_ds./max(psf2_ds);
figure,plot(psf1_ds),title('psf1 downsampled');
figure,plot(psf2_ds),title('psf2 downsampled');
tf1 = abs(F(psf1_ds));
tf2 = abs(F(psf2_ds));
figure,plot(tf1),title('tf1');
figure,plot(tf2),title('tf2');




