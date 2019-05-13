%%new script on GD only reconstuction with multiple measurements
%% randomly select illumination and measurements and see if it converges faster. i hope so.
%% Yujia Xue
img_dim = 64;
% img = 1-phantom('Modified Shepp-Logan',img_dim);
% img = 0.3+img./2; %%x in [0.3 0.8]
% imglena = double(imresize(imread('lena.png'),[img_dim,img_dim]));
% imglena = imglena./max(imglena(:));
% imglena = 0.3+imglena./2; %% range in 0.3 - 0.8
% x_mean = mean(img(:));
% x_mean_lena = mean(imglena(:));
img1 = imresize(double(rgb2gray(imread('layer1.png'))),[img_dim,img_dim]);
img1 = img1./max(img1(:));
img2 = imresize(double(rgb2gray(imread('layer2.png'))),[img_dim,img_dim]);
img2 = img2./max(img2(:));
img3 = imresize(double(rgb2gray(imread('layer3.png'))),[img_dim,img_dim]);
img3 = img3./max(img3(:));
Nslices = 3;
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
dz0 = 10;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.643;
pitch = 1;
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
thetax = [-9:3:9];
thetay = [-9:3:9];
measurements = zeros(img_row,img_col,length(thetax),length(thetay));
Prop = @(f,h) Ft(F(f).*h);
i0 = zeros(img_row,img_col,length(thetax),length(thetay));
i1 = zeros(img_row,img_col,length(thetax),length(thetay));
for i = 1:length(thetax)
    for j = 1:length(thetay)
        i0(:,:,i,j) = exp(1i*2*pi*(x*sin(thetax(i)*pi/180)/lambda+y*sin(thetay(j)*pi/180)/lambda));
        zinitial = dz0;
        Hinitial = exp(1i*k2*zinitial);
        i1(:,:,i,j) = Prop(i0(:,:,i,j),Hinitial);
        %%phi incident, psi output
        [ ~, psi ] = Fwd_Prop_MultiSlice_v2( i1(:,:,i,j), o_slices_t, k2, dz);
        measurements(:,:,i,j) = psi(:,:,end);
    end
end


dimy = [img_dim*img_dim,1];
fm = dftmtx2d(img_row,img_col);
ifm = conj(fm)./(img_row*img_col);
phasemat = exp(1i*k2*dz0);
Hmat = ifm*(diag(reshape(phasemat,dimy))*fm);
dimX = img_row*img_col*Nslices;
dimY = img_row*img_col;
outdir = datestr(now, 'yy_mm_DD_HH_MM');
outdir = ['randomselect_',outdir];
mkdir(outdir);
x_ini = repmat(x_mean,size(x0));
xn = x_ini;
pds0h = zeros(dimX,dimY);%% partial derivative of s0 which is illumination field so that not depend on object

maxiter =200;

alpha = 0.2;
delta = 0.00001;
%%store all history estimate of x
history_estimate = [];
%%main optimization GD 
for iter = 1 : maxiter
    tic;
    titer = mod(iter-1,length(thetax)*length(thetay))+1;
    iteri = floor((titer-1)/length(thetay))+1;
    iterj = mod(titer-1,length(thetay))+1;
%     iteri = round(rand(1)*length(thetax)+0.5);
%     iterj = round(rand(1)*length(thetay)+0.5);
    currenti0 = i0(:,:,iteri,iterj);
    currenti1 = i1(:,:,iteri,iterj);
    currentyobj = measurements(:,:,iteri,iterj);
    currenty0 = reshape(currentyobj,dimy);
    [ ~, psi ] = Fwd_Prop_MultiSlice_v2( currenti1, reshape(xn,size(o_slices_t)), k2, dz);
    yn = reshape(psi(:,:,end),dimy);
    residualn = yn-currenty0;
    pdsh = zeros(dimX,dimY,Nslices);
    %%update partial derivative of s
    skm1 = reshape(currenti0,dimY,1);
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
    gradienth = real(pdsh(:,:,end)*residualn);
    xnn = mytruncate(xn- alpha*gradienth,0.28,0.82);
    updateratio = norm(xnn-xn)/norm(xn);
    xn = xnn;
    history_estimate(:,size(history_estimate,2)+1) = xn;
    [ ~, psi ] = Fwd_Prop_MultiSlice_v2( currenti1, reshape(xn,size(o_slices_t)), k2, dz);
    yn = reshape(psi(:,:,end),dimy);
    residualn = yn-currenty0;
    mseerror = sum(abs(residualn.^2));
    tobj = reshape(xn,dim1);
    thologram = reshape(yn,img_row,img_col);
    for i = 1:Nslices
        figure(1),subplot(2,Nslices,i),imagesc(abs(tobj(:,:,i)),[0,1]),axis image, colormap gray,colorbar, title(['estimation of ',num2str(i),'th slice, after ',num2str(iter),'th iteration']);
        figure(1),subplot(2,Nslices,i+Nslices),imagesc(abs(xobj(:,:,i)),[0,1]),axis image,colormap gray,colorbar, title(['ground truth of ',num2str(i),'th slice']);
    end
    figure(2),subplot(1,2,1),imagesc(abs(currentyobj)),axis image,colorbar,colormap gray,title('real measurement');
    figure(2),subplot(1,2,2),imagesc(abs(thologram)),axis image,colorbar,colormap gray,title('output from current estimation');
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
    toc
end



