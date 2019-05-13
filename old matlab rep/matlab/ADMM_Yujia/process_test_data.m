img = imread('measure_0220.tif');
psf = imread('psf_0220.tif');
psf = double(psf);
psf = psf./sum(psf(:));
img = double(img)./250.0;

y = img;
save test_0220_1 y psf