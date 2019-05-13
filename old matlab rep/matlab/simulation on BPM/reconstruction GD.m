% %%reconstruction algorithm script
% 
% %%first start with basic gradient descent
% x = sin([0:pi/100:99*pi/100]');
% %%define the model as h , linear
% h = zeros(50,100);
% for i = 1:50
%     h(i,2*i-1:2*i) = 1/2;
% end
% 
% %%------------------
% y = h*x;
% alpha = 0.02;
% delta = 0.001;
% 
% mean_x = mean(x(:));
% x_0 = repmat(mean_x,size(x));
% x_n = x_0;
% for i = 1:1000
%     y_n = h*x_n;
%     e = y_n-y;
%     gradient = real(h'*e);
%     x_np1 = x_n - alpha*gradient;
%     x_n = x_np1;
%     en = h*x_n - y;
%     mseerror = sum(abs(en.^2));
%     figure(1),plot(x,'c'),hold on, plot(x_n,'r'),title(['iter = ',num2str(i),'alpha = ',num2str(alpha),'mse is: ',num2str(mseerror)]);
%     figure(2),plot(i,mseerror,'r+'),title('mse error'),hold on;
%     pause(0.03)
%     if mseerror <= delta
%         break;
%     end
% end
% 
% 
% 
% %%test gd on quadratic function
% x0 = sin([0.01:0.01:10]./10*2*pi)';
% y0 = x0.^2+2*x0+3;
% f = @(x) x.^2+2*x+3;
% x_ini = repmat(mean(y0),size(x0));
% alpha = 0.02;
% delta = 0.001;
% xn = x_ini;
% maxiter = 1000;
% for i = 1:maxiter
%     yn = f(xn);
%     e = yn-y0;
%     gradient = real((2*diag(xn)+2*eye(size(x0,1)))*(e));
%     xnp1 = xn - alpha*gradient;
%     xn = xnp1;
%     en = f(xn) - y0;
%     mseerror = sum(abs(en.^2));
%     figure(1),plot(x0,'c'),hold on, plot(xn,'r'),title(['iter = ',num2str(i),'alpha = ',num2str(alpha),'mse is: ',num2str(mseerror)]);
%     figure(2),plot(i,mseerror,'r+'),title('mse error'),hold on;
%     pause(0.03)
%     if mseerror <= delta
%         break;
%     end
% end











% %%test with phantom, ulugbek paper, gd only, no fista
% img_dim = 64;
% img = phantom('Modified Shepp-Logan',img_dim);
% img = img./max(img(:));
% x_mean = mean(img(:));
% Nslices = 1;
% dim1 = [img_dim,img_dim,Nslices];
% dim2 = [img_dim*img_dim*Nslices,1];
% xobj = zeros(dim1);
% xobj(:,:,1) = img;
% x0 = reshape(xobj,dim2);
% %%setup bpm and generate y0
% img_row = img_dim;
% img_col = img_dim;
% dz0 = 10;
% F = @(x) fftshift(fft2(ifftshift(x)));
% Ft = @(x) fftshift(ifft2(ifftshift(x)));
% lambda = 0.643;
% pitch = 1;
% Nslices = 1;
% o_slices_t = xobj;
% FOV = (pitch*img_row);
% du = 1/FOV;
% Umax = du * img_row/2;
% [u,v] = meshgrid(-Umax:du:Umax-du);
% [x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);
% k2 = pi*lambda*(u.^2+v.^2);
% %%define i0 as the field that is dz0 ahead of the first layer, which can
% %%be thought as 0th slice
% i0 = exp(1i*2*pi*(x*0+y*0));
% zinitial = dz0;
% Prop = @(f,h) Ft(F(f).*h);
% Hinitial = exp(1i*k2*zinitial);
% i1 = Prop(i0,Hinitial);
% %%phi incident, psi output
% % [ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, o_slices_t, k2, dz);
% yobj = i1.*o_slices_t;
% dimy = [64*64,1];
% y0 = reshape(yobj,dimy);
% fm = dftmtx2d(img_row,img_col);
% ifm = conj(fm)./(img_row*img_col);
% phasemat = exp(1i*k2*dz0);
% Hmat = ifm*(diag(reshape(phasemat,dimy))*fm);
% dimX = img_row*img_col*Nslices;
% dimY = img_row*img_col;
% %%initial x is mean of x, run forward model, get current output, current
% %%residual
% outdir = ['single_layer_',datestr(now, 'yy_mm_DD_HH_MM')];
% mkdir(outdir);
% x_ini = repmat(x_mean,size(x0));
% % x_ini = reshape(repmat(abs(yobj),[1,1,2]),size(x0));
% xn = x_ini;
% % [ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xn,size(o_slices_t)), k2, dz);
% psi = i1.*reshape(xn,size(o_slices_t));
% yn = reshape(psi(:,:,end),dimy);
% residualn = yn-y0;
% alpha = 0.1;
% delta = 0.01;
% maxiter =100;
% %%main optimization GD 
% for iter = 1 : maxiter
%     %%update partial derivative of s
%     gradienth = real(diag(reshape(i1,dim2))*residualn);
% %     xn = xn - alpha*gradienth;
%     xn = xn - alpha*gradienth;
% %     [ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xn,size(o_slices_t)), k2, dz);
%     yn = reshape(diag(reshape(i1,dim2))*xn,dimy);
%     residualn = yn-y0;
%     mseerror = sum(abs(residualn.^2));
%     tobj = reshape(xn,dim1);
%     thologram = reshape(yn,img_row,img_col);
%     for i = 1:Nslices
%         figure(1),subplot(2,Nslices,i),imagesc(abs(tobj(:,:,i))),colormap gray,title(['estimation of ',num2str(i),'th slice, after ',num2str(iter),'th iteration']);
%         figure(1),subplot(2,Nslices,i+Nslices),imagesc(abs(xobj(:,:,i))),colormap gray,title(['ground truth of ',num2str(i),'th slice']);
%     end
%     figure(2),subplot(1,2,1),imagesc(abs(yobj)),colormap gray,title('real measurement');
%     figure(2),subplot(1,2,2),imagesc(abs(thologram)),colormap gray,title('output from current estimation');
%     saveas(figure(1),[outdir,'/object_iter_',num2str(iter),'.png']);    
%     saveas(figure(2),[outdir,'/hologram_iter_',num2str(iter),'.png']);
%     figure(3),plot(iter,mseerror,'r+'),title('mse error'),hold on;
%     disp(iter);
%     disp(['mse error is: ',num2str(mseerror)]);
%     disp(['mean of gradient is: ',num2str(mean(gradienth(:)))]);
%     input('next');
%     if mseerror <= delta
%         break;
%     end
% 
% end






%%test with phantom, ulugbek paper, gd only, no fista
img_dim = 64;
% img = phantom('Modified Shepp-Logan',img_dim);
% img = double(imresize(imread('lena.png'),[img_dim,img_dim]));
% imga = imresize(double(rgb2gray(imread('a.png'))),[img_dim,img_dim]);
% imgb = imresize(double(rgb2gray(imread('b.png'))),[img_dim,img_dim]);
% imga = imga./max(imga(:));
% imgb = imgb./max(imgb(:));
img = 1-phantom('Modified Shepp-Logan',img_dim);
img = 0.3+img./2; %%x in [0.3 0.8]
x_mean = mean(img(:));
Nslices = 2;
dim1 = [img_dim,img_dim,Nslices];
dim2 = [img_dim*img_dim*Nslices,1];
xobj = zeros(dim1);
xobj(:,:,1) = img;
xobj(:,:,2) = img;
% xobj(:,:,2) = img;
% xobj(:,:,2) = ones(img_dim,img_dim);
x0 = reshape(xobj,dim2);
%%setup bpm and generate y0
img_row = img_dim;
img_col = img_dim;
dz0 =2;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.643;
pitch = 0.5;
o_slices_t = xobj;
dz = repmat(dz0,[1,Nslices-1]);
FOV = (pitch*img_row);
du = 1/FOV;
Umax = du * img_row/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);
k2 = pi*lambda*(u.^2+v.^2);
%%define i0 as the field that is dz0 ahead of the first layer, which can
%%be thought as 0th slice
i0 = exp(1i*2*pi*(x*sin(2*pi/180)/lambda+y*sin(0*pi/180)/lambda));
zinitial = dz0;
Prop = @(f,h) Ft(F(f).*h);
Hinitial = exp(1i*k2*zinitial);
i1 = Prop(i0,Hinitial);
%%phi incident, psi output
[ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, o_slices_t, k2, dz);
yobj = psi(:,:,end);
dimy = [img_dim*img_dim,1];
y0 = reshape(yobj,dimy);
fm = dftmtx2d(img_row,img_col);
ifm = conj(fm)./(img_row*img_col);
phasemat = exp(1i*k2*dz0);
Hmat = ifm*(diag(reshape(phasemat,dimy))*fm);
dimX = img_row*img_col*Nslices;
dimY = img_row*img_col;

%%initialize pds0h is partial derivative of S0 w.r.t. x which is zero
%%the partial derivative of s w.r.t. x at each slice is pdsh
pds0h = zeros(dimX,dimY);
pdsh = zeros(dimX,dimY,Nslices);
%%initial x is mean of x, run forward model, get current output, current
%%residual
outdir = datestr(now, 'yy_mm_DD_HH_MM');
mkdir(outdir);
x_ini = repmat(x_mean,size(x0));
% x_ini = reshape(repmat(abs(yobj),[1,1,Nslices]),size(x0));
% x_ini = x0 + 0.01*randn(size(x0));
xn = x_ini;
[ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xn,size(o_slices_t)), k2, dz);
yn = reshape(psi(:,:,end),dimy);
residualn = yn-y0;
alpha = 0.01;
delta = 0.0001;
maxiter =100;
%%main optimization GD 
for iter = 1 : maxiter
    %%update partial derivative of s
    tic;
    skm1 = reshape(i0,dimY,1);
    pdsm1h = pds0h;
    for i = 1:Nslices
        pdpkh = zeros(dimX,dimY);
        pdpkh((i-1)*dimY+1:(i-1)*dimY+dimY,:) = eye(dimY);
        tmp = pdsm1h*Hmat'*diag(conj(xn((i-1)*dimY+1:(i-1)*dimY+dimY)))+...
            pdpkh*diag(conj(Hmat*skm1));

        skm1 = reshape(psi(:,:,i),dimY,1);
        pdsh(:,:,i) = tmp;
        pdsm1h = tmp;
    end
    toc
    gradienth = real(pdsh(:,:,end)*residualn);
    xnn = mytruncate(xn- alpha*gradienth,0.2,0.9);
    updateratio = norm(xnn-xn)/norm(xn);
    xn = xnn;
    [ phi, psi ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xn,size(o_slices_t)), k2, dz);
    yn = reshape(psi(:,:,end),dimy);
    residualn = yn-y0;
    mseerror = sum(abs(residualn.^2));
    tobj = reshape(xn,dim1);
    thologram = reshape(yn,img_row,img_col);
    for i = 1:Nslices
        figure(1),subplot(2,Nslices,i),imagesc(abs(tobj(:,:,i)),[0,1]),colormap gray,colorbar, title(['estimation of ',num2str(i),'th slice, after ',num2str(iter),'th iteration']);
        figure(1),subplot(2,Nslices,i+Nslices),imagesc(abs(xobj(:,:,i)),[0,1]),colormap gray,colorbar, title(['ground truth of ',num2str(i),'th slice']);
    end
    figure(2),subplot(1,2,1),imagesc(abs(yobj)),colorbar,colormap gray,title('real measurement');
    figure(2),subplot(1,2,2),imagesc(abs(thologram)),colorbar,colormap gray,title('output from current estimation');
    saveas(figure(1),[outdir,'/object_iter_',num2str(iter),'.png']);    
    saveas(figure(2),[outdir,'/hologram_iter_',num2str(iter),'.png']);
    figure(3),plot(iter,mseerror,'r+'),title('mse error'),hold on;
    disp(['iter: ',num2str(iter)]);
    disp(['mse error is: ',num2str(mseerror)]);
    disp(['mean of gradient is: ',num2str(mean(gradienth(:)))]);
    disp(['update ratio is: ',num2str(updateratio)]);
%     input('next');
    if updateratio <= delta
        break;
    end
end



