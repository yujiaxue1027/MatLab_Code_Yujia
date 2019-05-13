%%reproduce dr ulugbek algorithm with only gradient descend
%%test with phantom, ulugbek paper, gd only, no fista
img_dim = 64;
% img = phantom('Modified Shepp-Logan',img_dim);
% img = double(imresize(imread('lena.png'),[img_dim,img_dim]));
% imga = imresize(double(rgb2gray(imread('a.png'))),[img_dim,img_dim]);
% imgb = imresize(double(rgb2gray(imread('b.png'))),[img_dim,img_dim]);
% imga = imga./max(imga(:));
% imgb = imgb./max(imgb(:));
% 
% %%phantom test
% img = 1-phantom('Modified Shepp-Logan',img_dim);
% img = 0.3+img./2; %%x in [0.3 0.8]
% x_mean = mean(img(:));
% Nslices = 2;
% dim1 = [img_dim,img_dim,Nslices];
% dim2 = [img_dim*img_dim*Nslices,1];
% xobj = zeros(dim1);
% xobj(:,:,1) = img;
% lena = double(imread('lena.png'));
% lena = lena./max(lena(:));
% lena = 0.3 + lena./2;
% xobj(:,:,2) = double(imresize((lena),[img_dim,img_dim]));

%%a few beads
img_dim = 64;
Nslices =3;
img1 = double(rgb2gray(imread('layer1.png')));
img1 = imresize(img1./255,[img_dim,img_dim]);
img2 = double(rgb2gray(imread('layer2.png')));
img2 = imresize(img2./255,[img_dim,img_dim]);
img3 = double(rgb2gray(imread('layer3.png')));
img3 = imresize(img3./255,[img_dim,img_dim]);
dim1 = [img_dim,img_dim,Nslices];
dim2 = [img_dim*img_dim*Nslices,1];
xobj = zeros(dim1);
xobj(:,:,1) = img1;
xobj(:,:,2) = img2;
xobj(:,:,3) = img3;
x_mean = mean(xobj(:));




% xobj(:,:,2) = img;
% xobj(:,:,2) = ones(img_dim,img_dim);
x0 = reshape(xobj,dim2);
%%setup bpm and generate y0
img_row = img_dim;
img_col = img_dim;
dz0 = 5;
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
thetax = 4;
thetay = 4;
i0 = exp(1i*2*pi*(x*sin(thetax*pi/180)/lambda+y*sin(thetay*pi/180)/lambda));
zinitial = dz0;
Prop = @(f,h) Ft(F(f).*h);
Hinitial = exp(1i*k2*zinitial);
i1 = Prop(i0,Hinitial);
%%phi incident, psi output
[ ~, psi ] = Fwd_Prop_MultiSlice_v2( i1, o_slices_t, k2, dz);
yobj = psi(:,:,end);
dimy = [img_dim*img_dim,1];
yk = reshape(yobj,dimy);
fm = dftmtx2d(img_row,img_col);
ifm = conj(fm)./(img_row*img_col);
phasemat = exp(1i*k2*dz0);
Hmat = ifm*(diag(reshape(phasemat,dimy))*fm);
dimX = img_row*img_col*Nslices;
dimY = img_row*img_col;
outdir = ['ulugbek_',datestr(now, 'yy_mm_DD_HH_MM')];
mkdir(outdir);
x_ini = repmat(x_mean,size(x0));
% x_ini = reshape(repmat(abs(yobj),[1,1,Nslices]),size(x0));
% x_ini = x0 + 0.01*randn(size(x0));
xhat = x_ini;
[ ~, yhat ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xhat,size(o_slices_t)), k2, dz);
yhatk = reshape(yhat(:,:,end),dimy);
alpha = 0.005;
delta = 0.000001;
maxiter =500;

rm = yhatk - yk;
sm = 0;
myestimate = zeros(dimX,maxiter);

for iter = 1 : maxiter
    
    for i = Nslices:-1:2
        parPparX = zeros(dimX,dimY);
        parPparX((i-1)*dimY+1:(i-1)*dimY+dimY,:) = eye(dimY);
        smm1 = sm + parPparX*diag(conj(Hmat*reshape(yhat(:,:,i-1),dimY,1)))*rm;
        rmm = Hmat'*diag(conj(xhat((i-1)*dimY+1:(i-1)*dimY+dimY)))*rm;
        sm = smm1;
        rm = rmm;
    end
    %%the reason to separate i = 1 is there is no yhat0 , only i0 the
    %%illumination
    i=1;
    parPparX = zeros(dimX,dimY);
    parPparX((i-1)*dimY+1:(i-1)*dimY+dimY,:) = eye(dimY);
    smm1 = sm + parPparX*diag(conj(Hmat*reshape(i0,dimY,1)))*rm;
    rmm = Hmat'*diag(conj(xhat((i-1)*dimY+1:(i-1)*dimY+dimY)))*rm;
    sm = smm1;
    rm = rmm;
    
    gradienth = real(smm1);
    xhatnew = mytruncate(xhat- alpha*gradienth,0.2,0.9);
    updateratio = norm(xhatnew-xhat)/norm(xhat);
    xhat = xhatnew;
    myestimate(:,iter) = xhat;
    [~,yhat ] = Fwd_Prop_MultiSlice_v2( i1, reshape(xhat,size(o_slices_t)), k2, dz);
    yhatk = reshape(yhat(:,:,end),dimy);
    rm = yhatk - yk;
	sm = 0;
    mseerror = sum(abs(rm.^2));
    tobj = reshape(xhat,dim1);
    thologram = reshape(yhatk,img_row,img_col);
    for i = 1:Nslices
        figure(1),subplot(2,Nslices,i),imagesc(abs(tobj(:,:,i))),colormap gray,colorbar, title(['estimation of ',num2str(i),'th slice, after ',num2str(iter),'th iteration']);
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

   

