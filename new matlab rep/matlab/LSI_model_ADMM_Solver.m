function [output,history] = LSI_model_ADMM_Solver(y,psf,para, varargin)
%% Yujia's ADMM implementation to solve crop and conv model
% y = CHx, y: measure; C:cropping operator; H:conv operator; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|CHx-y|^2+tau1*|Psix|_1+nonnega_indicator{x}+tau2*|x|_1

%% load prepared data
[rows,cols] = size(psf);
mu1 = para.mu1;
mu2 = para.mu2;
mu3 = para.mu3;
mu4 = para.mu4;
tau1_over_mu2 = para.tv_thres;
tau2_over_mu4 = para.l1_thres;
tau1 = tau1_over_mu2*mu2;
tau2 = tau2_over_mu4*mu4;
maxiter = para.maxiter;
color = para.color;
clip_min = para.clip_min;
clip_max = para.clip_max;
if length(varargin) == 3
    custom_display_region_flag = true;
    display_row_start = varargin{1};
    display_col_start = varargin{2};
    display_width = varargin{3};
else
    custom_display_region_flag = false;
end
%% ADMM parameters
history = zeros(2*rows,2*cols,maxiter);

%% define operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
C = @(x) x(1+rows/2:rows+rows/2,1+cols/2:cols+cols/2);
CT = @(x) padarray(x,[rows/2,cols/2]);
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
VEC = @(x) x(:);

%% ADMM algorithm
psf = psf./norm(psf(:));
Hs = F(CT(psf));
Hs_conj = conj(Hs);
H = @(x) real(Ft(F(x).*Hs));
HT = @(x) real(Ft(F(x).*Hs_conj));
HTH =abs(Hs.*Hs_conj);

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
f1 = figure('rend','painters','pos',[50 50 1500 900]);
f2 = figure('rend','painters','pos',[50 50 1500 900]);

while iter < maxiter
   iter = iter + 1;
   Hxt = Hxtp;
   vtp = v_mask.*(mu1*(gamma1/mu1 + Hxt) + CTy);
   wtp = clip(xt+gamma3/mu3, clip_min, clip_max);
   ztp = wthresh(gamma4/mu4 + xt,'s',tau2/mu4);
   
   [ut1,ut2] = soft_thres_2d(Psixt1+gamma2_1/mu2,Psixt2+gamma2_2/mu2,tau1/mu2);
   tmp_part1 = mu1*HT(vtp-gamma1/mu1);
   tmp_part2 = mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2);
   tmp_part3 = mu3*(wtp-gamma3/mu3);
   tmp_part4 = mu4*(ztp - gamma4/mu4);
   xtp_numerator = tmp_part1+...
       tmp_part2+...
       tmp_part3+...
       tmp_part4;
%    val1 = norm(norm(mu1*HT(vtp-gamma1/mu1)));
%    val2 = norm(norm(mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2)));
%    val3 = norm(norm(mu3*(wtp-gamma3/mu3)));
%    val4 = norm(norm(mu4*(ztp - gamma4/mu4)));
%    disp(['val1: ',num2str(val1),'  val2: ',num2str(val2),'  val3: ',num2str(val3),'  val4: ',num2str(val4)]);
   xtp = real(Ft(F(xtp_numerator).*x_mask));
   tmp_part1todisplay = (Ft(F(tmp_part1).*x_mask));
   tmp_part2todisplay = (Ft(F(tmp_part2).*x_mask));
   tmp_part3todisplay = (Ft(F(tmp_part3).*x_mask));
   tmp_part4todisplay = (Ft(F(tmp_part4).*x_mask));
   
   % update duals
   Hxtp = H(xtp);
   gamma1 = gamma1+mu1*(Hxtp-vtp);
   
   [Psixt1,Psixt2] = Psi(xtp);
   gamma2_1 = gamma2_1 + mu2*(Psixt1-ut1);
   gamma2_2 = gamma2_2 + mu2*(Psixt2-ut2);
   
   gamma3 = gamma3 + mu3*(xtp-wtp);

   gamma4 = gamma4 + mu4*(xtp-ztp);
  
   % update xt
   xt = xtp;  
   
   history(:,:,iter) = xt;
   
   % display and evaluate
   figure(f1);
   if custom_display_region_flag
       subplot(1,2,1),imagesc(xt(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),...
           colormap(color);axis image;colorbar;title(iter);
   else
       subplot(1,2,1),imagesc(xt),colormap(color);axis image;colorbar;title(iter);
   end
   residual = abs(C(H(xt))-y);
   dterm = 0.5*sum(residual(:).^2);
   [tv_x,tv_y] = Psi(xt);
   tv_x = [tv_x; zeros(1,2*cols)];
   tv_y = [tv_y, zeros(2*rows,1)];
   tvterm = tau1*sum(sqrt(tv_x(:).^2 + tv_y(:).^2));
   l1term = tau2*sum(abs(xt(:)));
   cost = dterm+tvterm+l1term;
   subplot(1,2,2),plot(iter,log10(cost),'bo'),grid on,hold on;...
            plot(iter,log10(dterm),'ro'),hold on;...
            plot(iter,log10(tvterm),'go'),hold on;...
            plot(iter,log10(l1term),'mo'),hold on;...
            title('log axis: blue: cost; red: data fidelity; green: tv; purple: l1');
   drawnow; 
   figure(f2);
   if custom_display_region_flag
       subplot(2,2,1),imagesc(tmp_part1todisplay(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['data: ',num2str(iter)]);
       subplot(2,2,2),imagesc(tmp_part2todisplay(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['tv: ',num2str(iter)]);
       subplot(2,2,3),imagesc(tmp_part3todisplay(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['NN: ',num2str(iter)]);
       subplot(2,2,4),imagesc(tmp_part4todisplay(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['l1: ',num2str(iter)]);
   else
       subplot(2,2,1),imagesc(tmp_part1todisplay),colormap(color);axis image;colorbar;title(['data: ',num2str(iter)]);
       subplot(2,2,2),imagesc(tmp_part2todisplay),colormap(color);axis image;colorbar;title(['tv: ',num2str(iter)]);
       subplot(2,2,3),imagesc(tmp_part3todisplay),colormap(color);axis image;colorbar;title(['NN: ',num2str(iter)]);
       subplot(2,2,4),imagesc(tmp_part4todisplay),colormap(color);axis image;colorbar;title(['l1: ',num2str(iter)]); 
   end
   drawnow;
   
   % print losses
   loss_cvy = C(vtp) - y;
   loss_cvy = 0.5 * sum(VEC(loss_cvy.^2));
   loss_tau1u = tau1 * ( sum(abs(VEC(ut1))) + sum(abs(VEC(ut2))));
   loss_tau2z = tau2 * sum(abs(VEC(ztp)));
   loss_mu1 = (Hxtp - vtp + gamma1/mu1);
   loss_mu1 = 0.5 * mu1 * sum(VEC(loss_mu1.^2));
   loss_mu2_1 = Psixt1 - ut1 + gamma2_1/mu2;
   loss_mu2_2 = Psixt2 - ut2 + gamma2_2/mu2;
   loss_mu2 = 0.5 * mu2 * (sum(VEC(loss_mu2_1.^2)) + sum(VEC(loss_mu2_2.^2)));
   loss_mu3 = xtp - wtp + gamma3/mu3;
   loss_mu3 = 0.5 * mu3 * sum(VEC(loss_mu3.^2));
   loss_mu4 = xtp - ztp + gamma4/mu4;
   loss_mu4 = 0.5 * mu4 * sum(VEC(loss_mu4.^2));
   loss_tot = loss_cvy + loss_tau1u + loss_tau2z + loss_mu1 + loss_mu2 + loss_mu3 + loss_mu4; 
   disp(['iter: ',num2str(iter),', cvy: ',num2str(round(loss_cvy,3)), ', tau1u: ',num2str(round(loss_tau1u,3))...
       ', tau2z: ', num2str(round(loss_tau2z,3)), ', mu1: ',  num2str(round(loss_mu1,3)), ', mu2: ',  num2str(round(loss_mu2,3))...
       , ', mu3: ',  num2str(round(loss_mu3,3)), ', mu4: ',  num2str(round(loss_mu4,3))...
       , ', total: ',  num2str(round(loss_tot,3))]);
   
end
output = xt;
end







