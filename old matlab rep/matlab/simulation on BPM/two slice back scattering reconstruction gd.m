%%two slices backscattering reconstruction with GD
img_dim = 64;
img = 1-phantom('Modified Shepp-Logan',img_dim);
img = 0.3+img./2; %%x in [0.3 0.8]
x_mean = mean(img(:));
Nslices = 2;
dim1 = [img_dim,img_dim,Nslices];
% dim2 = [img_dim*img_dim*Nslices,1];
xobj = zeros(dim1);
xobj(:,:,1) = img;
xobj(:,:,2) = img;
img_row = img_dim;
img_col = img_dim;
dz0 =10;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
lambda = 0.643;
pitch = 1;
o_slices_t = xobj;
o_slices_r = 1 - o_slices_t;
dz = repmat(dz0,[1,Nslices-1]);
FOV = (pitch*img_row);
du = 1/FOV;
Umax = du * img_row/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-img_row/2:img_row/2-1]*pitch);
k2 = pi*lambda*(u.^2+v.^2);

thetax = -6:3:6;
thetay = -6:3:6;
zinitial = dz0;
Prop = @(f,h) Ft(F(f).*h);
Hinitial = exp(1i*k2*zinitial);
s0 = zeros(img_dim,img_dim,length(thetax),length(thetay));
s1 = zeros(img_dim,img_dim,length(thetax),length(thetay));
y0 = zeros(img_dim,img_dim,length(thetax),length(thetay));
for i = 1:length(thetax)
    for j = 1:length(thetay)
        s0(:,:,i,j) = exp(1i*2*pi*(x*sin(thetax(i)*pi/180)/lambda+y*sin(thetay(j)*pi/180)/lambda));
        s1(:,:,i,j) = Prop(s0(:,:,i,j),Hinitial);
        [back_scattered_field, backfield] = Fwd_BackScattering_MultiSlice( s1(:,:,i,j), o_slices_t, k2, dz );
%         figure,imagesc(abs(back_scattered_field)),axis image, colormap gray, colorbar,title('total back field');
        y0(:,:,i,j) = back_scattered_field;
    end
end

outdir = datestr(now, 'yy_mm_DD_HH_MM');
outdir = ['two slices multi illum_',outdir];
mkdir(outdir);
x_ini = repmat(x_mean,size(o_slices_t));



history = [];
xn = x_ini;
xnr = 1 - xn;
H = exp(1i*k2*dz0);
alpha = 0.2;
delta = 0.00001;
maxiter = 200;
for iter = 1:maxiter
    %gradient
    tic;
%     titer = mod(iter-1,length(thetax)*length(thetay))+1;
%     iteri = floor((titer-1)/length(thetay))+1;
%     iterj = mod(titer-1,length(thetay))+1;
    iteri = round(rand(1)*length(thetax)+0.5);
    iterj = round(rand(1)*length(thetay)+0.5);
    [yn,~] = Fwd_BackScattering_MultiSlice( s1(:,:,iteri,iterj), xn, k2, dz );
    residualn = yn-y0(:,:,iteri,iterj);
    parBparx1 = zeros(img_dim,img_dim,img_dim,img_dim);
    parBparx2 = zeros(img_dim,img_dim,img_dim,img_dim);
    for i = 1:img_dim
        for j = 1:img_dim
            tparR = zeros(img_dim,img_dim);
            tparR(i,j)=-1;
            tparP = zeros(img_dim,img_dim);
            tparP(i,j)=1;
            parBparx1(:,:,i,j) = tparR.*(Prop(s0(:,:,iteri,iterj),H))+tparP.*Prop(xnr(:,:,2).*Prop(xn(:,:,1).*Prop(s0(:,:,iteri,iterj),H),H),H)+...
                xn(:,:,1).*Prop(xnr(:,:,2).*Prop(tparP.*Prop(s0(:,:,iteri,iterj),H),H),H);
            parBparx2(:,:,i,j) = xn(:,:,1).*Prop(tparR.*Prop(xn(:,:,1).*Prop(s0(:,:,iteri,iterj),H),H),H);
        end
%         disp(i);
    end
    residualn = reshape(residualn,[img_dim*img_dim,1]);
    parBh = zeros(img_dim*img_dim*2,img_dim*img_dim);
    for i = 1:img_dim
        for j = 1:img_dim
            idx1 = (j-1)*img_dim+i;
            tmp = parBparx1(:,:,i,j);
            parBh(idx1,:) = reshape(tmp,[1,img_dim*img_dim]);
            idx2 = (j-1)*img_dim+i+img_dim*img_dim;
            tmp = parBparx2(:,:,i,j);
            parBh(idx2,:) = reshape(tmp,[1,img_dim*img_dim]);
        end
    end
    xnn = mytruncate(xn - reshape(alpha*real(parBh*residualn),img_dim,img_dim,2),0.2,0.9);
    updateratio = norm(xnn(:)-xn(:))/norm(xn(:));
    xn = xnn;
    xnr = 1 - xn;
    [yn,~] = Fwd_BackScattering_MultiSlice( s1(:,:,iteri,iterj), xn, k2, dz );
    residualn = yn-y0(:,:,iteri,iterj);
    mseerror = sum(abs(residualn(:).^2));
    for i = 1:2
    figure(1),subplot(2,2,i),imagesc(abs(xn(:,:,i)),[0,1]),colormap gray, title(['estimation of ',num2str(i),'th slice, after ',num2str(iter),'th iteration']);
    figure(1),subplot(2,2,i+2),imagesc(abs(o_slices_t(:,:,i)),[0,1]),colormap gray, title(['ground truth of ',num2str(i),'th slice']);
    end
    figure(2),subplot(1,2,1),imagesc(abs(y0(:,:,iteri,iterj))),colorbar,colormap gray,title('real measurement');
    figure(2),subplot(1,2,2),imagesc(abs(yn)),colorbar,colormap gray,title('output from current estimation');
    saveas(figure(1),[outdir,'/object_iter_',num2str(iter),'.png']);    
    saveas(figure(2),[outdir,'/hologram_iter_',num2str(iter),'.png']);
    figure(3),plot(iter,mseerror,'r+'),title('mse error'),hold on;
    disp(['iter: ',num2str(iter)]);
    disp(['mse error is: ',num2str(mseerror)]);
    disp(['update ratio is: ',num2str(updateratio)]);
    if updateratio <= delta
        break;
    end
    toc
end








