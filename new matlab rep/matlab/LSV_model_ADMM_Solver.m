function [output,history] = LSV_model_ADMM_Solver(y,psfs,M,para, varargin)
%% Yujia's ADMM implementation for LSV deconvolution
% y = BMx, y: measure; B:[B1 B2 ... Bk] convolution matrix;
% M:[M1 M2 ... Mk]^T mask matrix; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|Bv-y|^2+tau1*|U|_1+indicator(w)+tau2*|z|_1...
% 0.5*mu1*|Mx-v+gamma1/mu1|^2 + 0.5*mu2*|Psix-u+gamma2/mu2|^2+...
% 0.5*mu3*|x-w+gamma3/mu3|^2 + 0.5*mu4*|x-z+gamma4/mu3|^2

%% ADMM parameters
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
[rows,cols, num_psfs] = size(psfs);
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
M_dagger_threshold = para.M_dagger_thres;
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

%% record history results
history = zeros(rows,cols,maxiter);

%% define operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));


%% ADMM algorithm
for i = 1:num_psfs
    psfs(:,:,i) = psfs(:,:,i)./norm(psfs(:,:,i));
end
TFs = F3D(psfs);
TFs_conj = conj(TFs);
B = @(x) real(sum(Ft3D(F3D(x).*TFs),3));
FBBT = sum(abs(TFs.*TFs_conj),3);

BTy = BT(y,TFs_conj);
xt = zeros(rows,cols);
gamma1 = zeros(rows,cols,num_psfs);
gamma3 = zeros(rows,cols);
gamma4 = zeros(rows,cols);

PsiTPsi = generate_laplacian_lsv_model(rows,cols);
gamma2_1 = zeros(rows-1,cols);
gamma2_2 = zeros(rows,cols-1);
PsiT = @(P1,P2) cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + ...
    cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));
Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2));
[ut1, ut2] = Psi(xt);
Psixt1 = ut1;
Psixt2 = ut2;
x_update_denominator = (mu1 + mu2*PsiTPsi + mu3 + mu4);

iter = 0;
f1 = figure('rend','painters','pos',[50 50 1500 900]);
f2 = figure('rend','painters','pos',[50 50 1500 900]);
while iter < maxiter
   iter = iter + 1;
   v_update_numerator = BTy + mu1*M.*repmat(xt,[1,1,num_psfs]) + gamma1;
   sum_x_tf = sum(F3D(v_update_numerator).*TFs,3);
   tmp_denominator = mu1 + FBBT;
   vtp =  Ft3D(F3D(v_update_numerator)-sum_x_tf.*TFs_conj./tmp_denominator)./mu1;  
   
   [ut1,ut2] = soft_thres_2d(Psixt1+gamma2_1/mu2,Psixt2+gamma2_2/mu2,tau1_over_mu2);
   
   wtp = clip(xt+gamma3/mu3, clip_min, clip_max);

   ztp = wthresh(xt+gamma4/mu4,'s',tau2_over_mu4);
   
   tmp_part1 = mu1*M_dagger(vtp-gamma1/mu1,M, M_dagger_threshold);
   tmp_part1(isnan(tmp_part1)) = 0;
   tmp_part2 = mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2);
   tmp_part3 = mu3*(wtp-gamma3/mu3);
   tmp_part4 = mu4*(ztp-gamma4/mu4);
   xtp_numerator = tmp_part1+...
       tmp_part2+...
       tmp_part3+tmp_part4;
   xtp = real(Ft(F(xtp_numerator)./x_update_denominator));
   
   % update duals
   gamma1 = gamma1+mu1*(M.*repmat(xtp,[1,1,num_psfs])-vtp);
   
   [Psixt1,Psixt2] = Psi(xtp);
   gamma2_1 = gamma2_1 + mu2*(Psixt1-ut1);
   gamma2_2 = gamma2_2 + mu2*(Psixt2-ut2);
   
   gamma3 = gamma3+mu3*(xtp-wtp);
   gamma4 = gamma4+mu4*(xtp-ztp);

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
   
   Mxhat = M.*repmat(xt,[1,1,num_psfs]);
   yhat = B(Mxhat);
   residual = abs(yhat-y);
   dterm = 0.5*sum(residual(:).^2);
   [tv_x,tv_y] = Psi(xt);
   tv_x = [tv_x; zeros(1,cols)];
   tv_y = [tv_y, zeros(rows,1)];
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
       subplot(2,2,1),imagesc(tmp_part1(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['data: ',num2str(iter)]);
       subplot(2,2,2),imagesc(tmp_part2(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['tv: ',num2str(iter)]);
       subplot(2,2,3),imagesc(tmp_part3(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['NN: ',num2str(iter)]);
       subplot(2,2,4),imagesc(tmp_part4(display_row_start:display_row_start+display_width-1,...
           display_col_start:display_col_start+display_width-1)),colormap(color);axis image;colorbar;title(['l1: ',num2str(iter)]);
   else
       subplot(2,2,1),imagesc(tmp_part1),colormap(color);axis image;colorbar;title(['data: ',num2str(iter)]);
       subplot(2,2,2),imagesc(tmp_part2),colormap(color);axis image;colorbar;title(['tv: ',num2str(iter)]);
       subplot(2,2,3),imagesc(tmp_part3),colormap(color);axis image;colorbar;title(['NN: ',num2str(iter)]);
       subplot(2,2,4),imagesc(tmp_part4),colormap(color);axis image;colorbar;title(['l1: ',num2str(iter)]);
   end
   drawnow;
   
   
%    disp(['D term: ',num2str(dterm)]);
%    disp(['TV: ',num2str(tvterm)]);
end
output = xt;
end