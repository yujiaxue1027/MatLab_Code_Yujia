function [output,costlist,out_paras] = admm_recon6(y,psf,para,rs,re,cs,ce)
%% Yujia's ADMM implementation to solve crop and conv model
% y = CHx, y: measure; C:cropping operator; H:conv operator; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|CHx-y|^2+tau1*|Psix|_1+nonnega_indicator{x}+tau2*|x|_1

%% load prepared data (groundtruth, y, psf)
[rows,cols] = size(psf);
tau1 = para.tau1;
tau2 = para.tau2;
mu1 = para.mu1;
mu2 = para.mu2;
mu3 = para.mu3;
mu4 = para.mu4;
maxiter = para.maxiter;
color = para.color;
mu_inc = para.mu_inc;
mu_dec = para.mu_dec;
mu_tol = para.mu_tol;



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
   wtp(1:rs,:)=0;
   wtp(re:end,:)=0;
   wtp(:,1:cs)=0;
   wtp(:,ce:end)=0;
   ztp = wthresh(gamma4/mu4 + xt,'s',tau2/mu4);
   
   [ut1,ut2] = soft_thres_2d(Psixt1+gamma2_1/mu2,Psixt2+gamma2_2/mu2,tau1/mu2);
   xtp_numerator = mu1*HT(vtp-gamma1/mu1)+...
       mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2)+...
       mu3*(wtp-gamma3/mu3)+...
       mu4*(ztp - gamma4/mu4);
   val1 = norm(norm(mu1*HT(vtp-gamma1/mu1)));
   val2 = norm(norm(mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2)));
   val3 = norm(norm(mu3*(wtp-gamma3/mu3)));
   val4 = norm(norm(mu4*(ztp - gamma4/mu4)));
   disp(['val1: ',num2str(val1),'  val2: ',num2str(val2),'  val3: ',num2str(val3),'  val4: ',num2str(val4)]);
   xtp = real(Ft(F(xtp_numerator).*x_mask));
   
   % update duals
   Hxtp = H(xtp);
   gamma1 = gamma1+mu1*(Hxtp-vtp);
   
   [Psixt1,Psixt2] = Psi(xtp);
   gamma2_1 = gamma2_1 + mu2*(Psixt1-ut1);
   gamma2_2 = gamma2_2 + mu2*(Psixt2-ut2);
   
   gamma3 = gamma3 + mu3*(xtp-wtp);

   gamma4 = gamma4 + mu4*(xtp-ztp);
   

   
   % update mu1 mu2 mu3 mu4
   
   mu1_p_r = norm(Hxtp(:)-vtp(:));
   mu1_d_r = mu1*norm(Hxtp(:)-Hxt(:));
   mu1 = ADMM_update_mu(mu1,mu_tol,mu_inc,mu_dec,mu1_p_r,mu1_d_r);
   
   mu2_p_r_1 = norm(Psixt1(:)-ut1(:));
   mu2_p_r_2 = norm(Psixt2(:)-ut2(:)); 
   mu2_p_r = sqrt(mu2_p_r_1^2 + mu2_p_r_2^2);
   [last_Psixt1,last_Psixt2] = Psi(xt);
   mu2_d_r_1 = norm(Psixt1(:)-last_Psixt1(:));
   mu2_d_r_2 = norm(Psixt2(:)-last_Psixt2(:));  
   mu2_d_r = mu2*sqrt(mu2_d_r_1^2 + mu2_d_r_2^2);
   mu2 = ADMM_update_mu(mu2,mu_tol,mu_inc,mu_dec,mu2_p_r,mu2_d_r);

   mu3_p_r = norm(xtp(:)-wtp(:));
   mu3_d_r = mu3*norm(xtp(:)-xt(:));
   mu3 = ADMM_update_mu(mu3,mu_tol,mu_inc,mu_dec,mu3_p_r,mu3_d_r);

   mu4_p_r = norm(xtp(:)-ztp(:));
   mu4_d_r = mu4*norm(xtp(:)-xt(:));   
   mu4 = ADMM_update_mu(mu4,mu_tol,mu_inc,mu_dec,mu4_p_r,mu4_d_r);

   % update masks
   x_mask = 1./(mu1*HTH + mu2*PsiTPsi + mu3 + mu4);  %Denominator of x update (diagonized in FFT domain)
   v_mask = 1./(CT(ones(size(y))) + mu1); %Denominator of v update (itself is diagonized)

   % update xt
   xt = xtp;  
   
   % display and evaluate
   if mu_inc~=1 || mu_dec~=1
          subplot(1,3,1),imagesc(xt),colormap(color);axis image;colorbar;title(iter);
       residual = abs(C(H(xt))-y);
       dterm = 0.5*sum(residual(:).^2);
       [tv_x,tv_y] = Psi(xt);
       tvterm = tau1*sqrt(sum(tv_x(:).^2)+sum(tv_y(:).^2));
       l1term = tau2*sum(abs(xt(:)));
       cost = dterm+tvterm+l1term;
       subplot(1,3,2),plot(iter,log10(cost),'bo'),grid on,hold on;...
                plot(iter,log10(dterm),'ro'),hold on;...
                plot(iter,log10(tvterm),'go'),hold on;...
                plot(iter,log10(l1term),'mo'),hold on;...
                title('log axis: blue: cost; red: data fidelity; green: tv; purple: l1');
       subplot(1,3,3),plot(iter,log10(mu1),'bo'),grid on,hold on;...
           plot(iter,log10(mu2),'ro'),hold on;...
           plot(iter,log10(mu3),'go'),hold on;...
           plot(iter,log10(mu4),'mo'),hold on;...
           title('log axis: blue: mu1; red: mu2; green: mu3; purple: mu4');
       drawnow; 
   else
       subplot(1,2,1),imagesc(xt),colormap(color);axis image;colorbar;title(iter);
   residual = abs(C(H(xt))-y);
   dterm = 0.5*sum(residual(:).^2);
   [tv_x,tv_y] = Psi(xt);
   tvterm = tau1*sqrt(sum(tv_x(:).^2)+sum(tv_y(:).^2));
   l1term = tau2*sum(abs(xt(:)));
   cost = dterm+tvterm+l1term;
   subplot(1,2,2),plot(iter,log10(cost),'bo'),grid on,hold on;...
            plot(iter,log10(dterm),'ro'),hold on;...
            plot(iter,log10(tvterm),'go'),hold on;...
            plot(iter,log10(l1term),'mo'),hold on;...
            title('log axis: blue: cost; red: data fidelity; green: tv; purple: l1');
   drawnow; 
   end
   % print mu1 mu2 mu3 mu4 every 10 iters
   if mod(iter,10) == 1
       disp(['mu1: ',num2str(mu1),' mu2: ',num2str(mu2),...
           ' mu3: ',num2str(mu3),' mu4: ',num2str(mu4)]);
   end
end
output = xt;
costlist.tot = cost;
costlist.tv = tvterm;
costlist.l1 = l1term;
out_paras.mu1 = mu1;
out_paras.mu2 = mu2;
out_paras.mu3 = mu3;
out_paras.mu4 = mu4;
out_paras.tau1 = tau1;
out_paras.tau2 = tau2;
end