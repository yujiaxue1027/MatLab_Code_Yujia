function [output,history] = admm_lsv(y,psfs,M,para)
%% Yujia's ADMM implementation for LSV deconvolution
% y = BMx, y: measure; B:[B1 B2 ... Bk] convolution matrix;
% M:[M1 M2 ... Mk]^T mask matrix; x: obj
% apply admm to solve: xhat = argmin
% 0.5*|Bv-y|^2+tau*|U|_1+0.5*mu1*|Mx-v+gamma1/mu1|^2 + 
% 0.5*mu2*|Psix-u+gamma2/mu2|^2

%% ADMM parameters
[rows,cols, num_psfs] = size(psfs);
tau = para.tau;
mu1 = para.mu1;
mu2 = para.mu2;
maxiter = para.maxiter;
color = para.color;

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

xt = zeros(rows,cols);
gamma1 = zeros(rows,cols,num_psfs);
BTy = BT(y,TFs_conj);

PsiTPsi = generate_laplacian_lsv_model(rows,cols);
gamma2_1 = zeros(rows-1,cols);
gamma2_2 = zeros(rows,cols-1);
PsiT = @(P1,P2) cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + ...
    cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));
Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2));
[ut1, ut2] = Psi(zeros(rows, cols));
Psixt1 = ut1;
Psixt2 = ut2;
x_update_denominator = (mu1 + mu2*PsiTPsi);

iter = 0;
figure;
while iter < maxiter
   iter = iter + 1;
   v_update_numerator = BTy + mu1*M.*repmat(xt,[1,1,num_psfs]) + gamma1;
   sum_x_tf = sum(F3D(v_update_numerator).*TFs,3);
   tmp_denominator = mu1 + FBBT;
   vtp =  Ft3D(F3D(v_update_numerator)-sum_x_tf.*TFs_conj./tmp_denominator)./mu1;  
   
   [ut1,ut2] = soft_thres_2d(Psixt1+gamma2_1/mu2,Psixt2+gamma2_2/mu2,tau/mu2);
   
   xtp_numerator = mu1*M_dagger(vtp-gamma1/mu1,M)+...
       mu2*PsiT(ut1-gamma2_1/mu2,ut2-gamma2_2/mu2);
   xtp = real(Ft(F(xtp_numerator)./x_update_denominator));
   
   % update duals
   gamma1 = gamma1+mu1*(M.*repmat(xtp,[1,1,num_psfs])-vtp);
   
   [Psixt1,Psixt2] = Psi(xtp);
   gamma2_1 = gamma2_1 + mu2*(Psixt1-ut1);
   gamma2_2 = gamma2_2 + mu2*(Psixt2-ut2);
   

   % update xt
   xt = xtp;  
   
   
   history(:,:,iter) = xt;
   % display and evaluate

   subplot(1,2,1),imagesc(xt),colormap(color);axis image;colorbar;title(iter);
   Mxhat = M.*repmat(xt,[1,1,num_psfs]);
   yhat = B(Mxhat);
   residual = abs(yhat-y);
   dterm = 0.5*sum(residual(:).^2);
   [tv_x,tv_y] = Psi(xt);
   tvterm = tau*sqrt(sum(tv_x(:).^2)+sum(tv_y(:).^2));
   cost = dterm+tvterm;
   subplot(1,2,2),plot(iter,log10(cost),'bo'),grid on,hold on;...
            plot(iter,log10(dterm),'ro'),hold on;...
            plot(iter,log10(tvterm),'go'),hold on;...
            title('log axis: blue: cost; red: data fidelity; green: tv');
   drawnow; 
   
end
output = xt;
end