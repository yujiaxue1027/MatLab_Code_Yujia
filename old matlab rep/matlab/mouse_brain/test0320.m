psf = im2double(imread('0320/psf.tif'));
bg = im2double(imread('0320/bg.tif'));
y1 = im2double(imread('0320/y1.tif'));
y2 = im2double(imread('0320/y2.tif'));
psf = psf-bg;
psf(psf<=0) = 0;

psf = average_shrink(psf,4);
y1 = average_shrink(y1,4);
y2 = average_shrink(y2,4);
y1 = y1./max(y1(:));
y2 = y2./max(y2(:));
para= [];
para.tau1= 1.0000e-03;
para.tau2= 0.1000;
para.mu1= 0.0050;
para.mu4= 0.5000;
para.mu2= 0.0500;
para.mu3= 0.5000;
para.mu_inc= 1;
para.mu_dec= 1;
para.mu_tol= 2;
para.maxiter= 55;
para.color= 'gray';

[recon1,~,~] = admm_recon4(y1,psf,para);
[recon2,~,~] = admm_recon4(y2,psf,para);
