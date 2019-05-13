function output = admm_recon2(data_file_name,tau1,tau2,mu1,mu2,mu3,mu4,maxiter,color)
%% Yujia's ADMM implementation to solve crop and conv model
% y = CHx, y: measure; C:cropping operator; H:conv operator; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|CHx-y|^2+tau1*|Psix|_1+nonnega_indicator{x}+tau2*|x|_1

%% load prepared data (groundtruth, y, psf)
load(data_file_name);
[rows,cols] = size(psf);


%% ADMM parameters


%% define operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
C = @(x) x(1+rows/2:rows+rows/2,1+cols/2:cols+cols/2);
CT = @(x) padarray(x,[rows/2,cols/2]);

%% ADMM algorithm
psf = psf./norm(psf(:));
Hs = F(CT(psf));
Hs_conj = conj(Hs);
H = @(x) real(Ft(F(x).*Hs));
HT = @(x) real(Ft(F(x).*Hs_conj));
HTH = abs(Hs.*Hs_conj);

xt = zeros(2*rows,2*cols);
gamma1 = zeros(2*rows,2*cols);
gamma3 = zeros(2*rows,2*cols);
gamma4 = zeros(2*rows,2*cols);
CTy = CT(y);

PsiTPsi = generate_laplacian(rows,cols);
gamma2_1 = zeros(2*rows-1,2*cols);
gamma2_2 = zeros(2*rows,2*cols-1);
PsiT = @(P1,P2) cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + ...
    cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));
Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2));
[ut1, ut2] = Psi(zeros(2*rows, 2*cols));
Psixt1 = ut1;
Psixt2 = ut2;

x_mask = 1./(mu1*HTH + mu2*PsiTPsi + mu3 + mu4);  %Denominator of x update (diagonized in FFT domain)
v_mask = 1./(CT(ones(size(y))) + mu1); %Denominator of v update (itself is diagonized)

iter = 0;
Hxtp = zeros(2*rows,2*cols);
figure;
while iter < maxiter
   iter = iter + 1;
   Hxt = Hxtp;
   vtp = v_mask.*(mu1*(gamma1/mu1 + Hxt) + CTy);
   wtp = max(gamma3/mu3 + xt,0);
   ztp = wthresh(gamma4/mu4 + xt,'s',tau2/mu4);
   
   [ut1,ut2] = soft_thres_2d(Psixt1+gamma2_1/mu2,Psixt2+gamma2_2/mu2,tau1/mu2);
   xtp_numerator = mu1*HT(vtp-gamma1/mu1)+...
       mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2)+...
       mu3*(wtp-gamma3/mu3)+...
       mu4*(ztp - gamma4/mu4);
   xtp = real(Ft(F(xtp_numerator).*x_mask));
   
   % update duals
   Hxtp = H(xtp);
   gamma1 = gamma1+mu1*(Hxtp-vtp);
   
   [Psixt1,Psixt2] = Psi(xtp);
   gamma2_1 = gamma2_1 + mu2*(Psixt1-ut1);
   gamma2_2 = gamma2_2 + mu2*(Psixt2-ut2);
   
   gamma3 = gamma3 + mu3*(xtp-wtp);

   gamma4 = gamma4 + mu4*(xtp-ztp);
   
   xt = xtp;
   
   % display and evaluate
   subplot(1,2,1),imagesc(xt),colormap(color);axis image;colorbar;title(iter);
   residual = C(H(xt))-y;
   mse = norm(residual(:));
   subplot(1,2,2),plot(iter,mse,'bo'),grid on,hold on;...
            title('admm: norm of residual');
   drawnow; 
end
output = xt;

end