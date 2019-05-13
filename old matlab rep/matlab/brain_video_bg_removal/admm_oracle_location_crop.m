function [output,history] = admm_oracle_location_crop(y,psf,mask,para,plot_flag)
%% assume LSI forward model, know oracle location of individual neurons

%% load prepared data (groundtruth, y, psf)
[rows,cols] = size(psf);
mu1 = para.mu1;
mu2 = para.mu2;
maxiter = para.maxiter;
color = para.color;



%% ADMM parameters
history = zeros(rows*2,cols*2,maxiter);

%% define operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);

%% ADMM algorithm
psf = psf./norm(psf(:));
Hs = F(pad2d(psf));
Hs_conj = conj(Hs);
H = @(x) real(Ft(F(x).*Hs));
HT = @(x) real(Ft(F(x).*Hs_conj));
HTH = abs(Hs.*Hs_conj);

mask = pad2d(mask);
xt = zeros(2*rows,2*cols);
gamma1 = zeros(2*rows,2*cols);
gamma2 = zeros(2*rows,2*cols);
CTy = pad2d(y);
Hxtp = zeros(2*rows,2*cols);

x_mask = 1./(mu1*HTH + mu2);  %Denominator of x update (diagonized in FFT domain)
v_mask = 1./(pad2d(ones(size(y))) + mu1);


iter = 0;
if plot_flag
    figure;
end
while iter < maxiter
   iter = iter + 1;
   Hxt = Hxtp;
   vtp = v_mask.*(mu1*(gamma1/mu1 + Hxt) + CTy);
   
   wtp = max(gamma2/mu2 + xt,0);
   wtp(mask~=1) = 0;
   
   xtp_numerator = mu1*HT(vtp-gamma1/mu1)+mu2*(wtp-gamma2/mu2);

   xtp = real(Ft(F(xtp_numerator).*x_mask));
   
   Hxtp = H(xtp);
   % update duals
   gamma1 = gamma1+mu1*(Hxtp-vtp);
   gamma2 = gamma2+mu2*(xtp-wtp);
   
   % update xt
   xt = xtp;  
   
   
   history(:,:,iter) = xt;
   % display and evaluate
   if plot_flag
        subplot(1,2,1),imagesc(xt),colormap(color);axis image off;colorbar;title(iter);
        residual = (crop2d(H(xt))-y);
        dterm = 0.5*sum(residual(:).^2);
        subplot(1,2,2),plot(iter,(dterm),'ro'),hold on;...
                    title('red: data fidelity');
        drawnow; 
   end
  
end
output = xt;
end