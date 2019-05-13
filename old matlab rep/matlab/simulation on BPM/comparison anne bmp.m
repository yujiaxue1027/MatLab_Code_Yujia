%% script compare classic bpm and anne's method
%%simulation setup, unit: um
rows = 256;
cols = 256;
nslices = 256;
delta = 0.144;
lambda0 = 0.633;
k0 = 2*pi/lambda0;  %%in free space
n0 = 1.518;
n1 = 1.548;
lambda = lambda0/n0;
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
dz0 = delta;
%% setup the bead, radius 3.5
n_obj = n0*ones(rows,cols,nslices);
dz = repmat(dz0,[1,nslices-1]);
tcenter = [(cols+1)/2*delta,(1+rows)/2*delta,(nslices+1)/2*dz0];
radius = 3.5;
for j = 1:cols
    for i = 1:rows
        for k = 1:nslices
            tpos = [j*delta,i*delta,k*dz0];
            tdist = norm(tpos-tcenter);
            if tdist <= radius
                n_obj(i,j,k) = n1;
            end
        end
    end
end
%% object function in classic bpm (ignore higher order terms of delta_n)
delta_n = n_obj - n0;
object = exp(1i*k0*dz0*delta_n);
%%some setup parameters
FOV = (delta*rows);
du = 1/FOV;
Umax = du * rows/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-rows/2:rows/2-1]*delta);
flim = 1/lambda;
tmpuv = u.^2+v.^2;
pupil = double(1.1*tmpuv<(flim^2));
% k2 = n0*pi*lambda*(u.^2+v.^2);
%%define i0 as the field that is dz0 ahead of the first layer, which can
%%be thought as 0th slice
thetax = 0;
thetay = 0;
i0 = exp(1i*2*pi*(x*thetax+y*thetay));
zinitial = dz0;
Prop = @(f,h,pupil) Ft(F(f).*h.*pupil);
% Prop = @(f,h) Ft(F(f).*h);
k2 = real(-2*pi*(1/lambda-sqrt(1/(lambda*lambda)-(u.^2+v.^2))));
% Hinitial = exp(-1i*2*pi*zinitial*(1/lambda-sqrt(1/(lambda*lambda)-(u.^2+v.^2))));
Hinitial = exp(1i*k2*zinitial);
i1 = Prop(i0,Hinitial,pupil);
[ phi, psi ] = yujiaforward( i1, object, k2, dz,pupil);
transmission = psi(:,:,end);


Hprime = exp(1i*k2*(-13));
farfieldbpm = Prop(transmission,Hprime,pupil);
figure,imagesc(abs(farfieldbpm-1));





%% summation implementation
% L = nslices*delta/2;
% u_in = ones(rows,cols);%% u_in, all one, because of normal incidence
% xx = [rows/2-1:-1:-rows/2]*delta;
% yy = [cols/2-1:-1:-cols/2]*delta;
% zz = L;
% u_scat = zeros(rows,cols);
% F = (k0^2/4/pi)*((n_obj).^2-1);
% HH = F.*psi;
% xcor = xx;
% ycor = yy;
% zcor = [-nslices/2:nslices/2-1]*dz0;
% for i = 1:rows/2
%     for j = 1:cols/2
%         tmp = 0;
%         r = [xx(1,i),yy(1,j),zz];
%         for ii = 128-30:128+30
%             for jj = 128-30:128+30
%                 for kk = 128-30:128+30
%                     tpos = [jj*delta,ii*delta,kk*dz0];
%                     tdist = norm(tpos-tcenter);
%                     if tdist <= radius
%                         rprime = [xcor(1,ii),ycor(1,jj),zcor(1,kk)];
%                         gtmp = exp(1i*k0*norm(r-rprime))/norm(r-rprime);
%                         tmp = tmp + HH(ii,jj,kk)*gtmp;
%                     end
%                 end
%             end
%         end
%         u_scat(i,j) = tmp;
%     end
%     disp(i);
% end

%% FFT implementation (waleed's derivation)
obj3d = n_obj.^2 - n0^2;
pupil3d = double(n_obj == n1);
O = obj3d.*pupil3d;
U = psi;
OU = O.*U;
ztarget = 18;
z =  [-nslices/2:nslices/2-1]*dz0;
tmpsum = zeros(rows,cols);


%% 55555
tmpsum = zeros(rows,cols);
% ttt = zeros(rows,cols,nslices);
for k = 1:nslices
    tmp = F(O(:,:,k).*exp(1i*2*pi/lambda*z(k)));
%     tmp = F(OU(:,:,k));
    factor = exp(1i*2*pi*sqrt(1/(lambda0*lambda0)-(u.^2+v.^2)).*abs(z(k)-ztarget))./sqrt(1/(lambda0*lambda0)-(u.^2+v.^2));
    tmp = tmp.*factor.*dz0.*pupil;
    tmpsum = tmpsum + tmp;
%     ttt(:,:,k) = tmp;
%     disp(k);
end
u_scat = 1i*pi/(lambda0^2)*Ft(tmpsum);
figure,imagesc(abs(u_scat)),axis image off; 






        
        
        
        
        