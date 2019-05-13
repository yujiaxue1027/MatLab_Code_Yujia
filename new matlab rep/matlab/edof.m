%% unit in um
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

lambda = 0.5;
p = 1;
[x,y] = meshgrid([-1000:999]*p);
fov = 2000*p;
du = 1/fov;
[u,v] = meshgrid([-1000:999]*du);
k2 = pi*lambda*(u.^2+v.^2);

i0 = ones(2000);
f = 10000; %% focal length 10 mm
lens = pi.*(x.^2+y.^2)/lambda/f; %% lens phase
gtmp = i0.*exp(1i*lens); %% field after the lens

Prop = @(f,h) Ft(F(f).*h);
H = exp(1i*k2*f); %% propagate for 10 mm , expect to see a single point on sensor
gout = Prop(gtmp, H);
figure,imagesc(abs(gout));

Prop = @(f,h) Ft(F(f).*h);
H2 = exp(1i*k2*8000); %% propagate for 8 mm , expect to see blurring
gout2 = Prop(gtmp, H2);
figure,imagesc(abs(gout2));


