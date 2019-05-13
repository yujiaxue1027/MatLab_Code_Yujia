%% simulate pinhole setup

%% define functions
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
Prop = @(f,h,pupil) Ft(F(f).*h.*pupil);
H = @(k2,z) exp(1i*k2*z);

%% basic parameters
input_beam_radius = 3751;
efl = 16000;
lambda = 0.632;
pitch = 1;
num_samples = 15000;
fov = pitch * num_samples;
du = 1/fov;
Umax = du * num_samples/2;
[u,v] = meshgrid(-Umax:du:Umax-du);
[x,y] = meshgrid([-num_samples/2:num_samples/2-1]*pitch);
flim = 1/lambda;
tmpuv = u.^2+v.^2;
pupil = double(1.1*tmpuv<(flim^2));
k2 = real(-2*pi*(1/lambda-sqrt(1/(lambda*lambda)-(u.^2+v.^2))));

%% simulation
u0 = zeros(num_samples);
u0(sqrt(x.^2+y.^2)<=input_beam_radius) = 1;
f_infocus = F(u0);
pinhole = zeros(num_samples);
pinhole_diameter = 50;
pinhole(sqrt(x.^2+y.^2)<=(pinhole_diameter/2)) = 1;

%%
screen_pos = 20000;
u1_no_pinhole = Prop(f_infocus, H(k2, screen_pos),pupil);
figure,imagesc(abs(u1_no_pinhole)),title('no pinhole');
u1_pinhole = Prop(f_infocus.*pinhole, H(k2, screen_pos),pupil);
figure,imagesc(abs(u1_pinhole)),title('in focus pinhole');

displace = 500;
tmp = Prop(f_infocus, H(k2, displace),pupil);
tmp = tmp.*pinhole;
u2_pinhole = Prop(tmp, H(k2, screen_pos),pupil);
figure,imagesc(abs(u2_pinhole)),title('displaced');


