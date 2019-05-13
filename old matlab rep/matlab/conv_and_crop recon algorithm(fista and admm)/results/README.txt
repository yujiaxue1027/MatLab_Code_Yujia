The reconstruction result from ADMM has scaling problem. In order to compute the PSNR for both FISTA and ADMM, I normalize groundtruth, recon_admm and recon_fista by first subtract the mean and divided by its standard deviation.

Both algorithm use the same input and same psf.

PSNR:
FISTA: 22.9db
ADMM: 27.8db

I compute PSNR as: psnr(gt,recon,max,M,N) = 20*log10(max/(sum((gt(:)-recon(:)).^2)/(M*N)))

gt: groundtruth, recon: reconstruction, max: max of gt - min of gt, M: rows, N: cols