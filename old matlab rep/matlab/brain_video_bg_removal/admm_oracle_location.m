function [output,history] = admm_oracle_location(y,psf,mask,para,plot_flag)
%% assume LSI forward model, know oracle location of individual neurons

%% load prepared data (groundtruth, y, psf)
[rows,cols] = size(psf);
mu = para.mu;
maxiter = para.maxiter;
color = para.color;



%% ADMM parameters
history = zeros(rows,cols,maxiter);

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
H = @(x) real(crop2d(Ft(F(pad2d(x)).*Hs)));
HT = @(x) real(crop2d(Ft(F(pad2d(x)).*Hs_conj)));
HTH = abs(Hs.*Hs_conj);

xt = zeros(rows,cols);
gamma = zeros(rows,cols);
HTy = HT(y);

x_mask = 1./(HTH + mu);  %Denominator of x update (diagonized in FFT domain)

iter = 0;
if plot_flag
    figure;
end
while iter < maxiter
   iter = iter + 1;
   wtp = max(gamma/mu + xt,0);
   wtp(mask~=1) = 0;
   xtp_numerator = HTy + mu*wtp - gamma;
   xtp = real(crop2d(Ft(F(pad2d(xtp_numerator)).*x_mask)));
   
   % update duals
   gamma = gamma+mu*(xtp-wtp);
   
   % update xt
   xt = xtp;  
   
   
   history(:,:,iter) = xt;
   % display and evaluate
   if plot_flag
        subplot(1,2,1),imagesc(xt),colormap(color);axis image off;colorbar;title(iter);
        residual = (H(xt)-y);
        dterm = 0.5*sum(residual(:).^2);
        subplot(1,2,2),plot(iter,(dterm),'ro'),hold on;...
                    title('red: data fidelity');
        drawnow; 
   end
  
end
output = xt;
end