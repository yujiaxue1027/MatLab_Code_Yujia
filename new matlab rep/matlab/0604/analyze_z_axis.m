%% analyize z-axis decorrelation
VEC = @(x) x(:);
my_mean = @(x) mean(VEC(x));
zero_mean = @(x) x-my_mean(x);
my_cov = @(x,y) my_mean(zero_mean(x).*zero_mean(y));
my_pcc = @(x,y) my_cov(x,y)/(std(VEC(x))*std(VEC(y)));
src_folder = 'psf_zscan/';
num_scans = 61;
rows = 1944;
cols = 2592;
psfs = zeros(rows,cols,num_scans);
for i = 1:num_scans
    tmp = im2double(imread([src_folder,num2str(i,'%.2d'),'.tif']));
    psfs(:,:,i) = tmp;
end

%% plot z-axis magnitude, variance
magnitude = zeros(num_scans,1);
variance = zeros(num_scans,1);
for i = 1:num_scans
    magnitude(i) = norm(VEC(psfs(:,:,i)));
    variance(i) = std(VEC(psfs(:,:,i)));
end
z_co = 3:0.1:9;
figure;
subplot(1,2,1);
plot(z_co, magnitude);
title('magnitude');
grid on;
subplot(1,2,2);
plot(z_co, variance);
title('variance');
grid on;

%% plot correlation C(in_focus_psf, psf(z)), in_focus is idx. 43
pcc_z = zeros(num_scans,1);
ref = psfs(:,:,43);
for i = 1:num_scans
    pcc_z(i) = my_pcc(ref,psfs(:,:,i));
end
pcc_single = zeros(num_scans,1);
ref_single = psfs(801:1100,1801:2100,43);
for i = 1:num_scans
    pcc_single(i) = my_pcc(ref_single,psfs(801:1100,1801:2100,i));
end

figure,
plot(z_co, pcc_z,'b');
hold on;
plot(z_co,pcc_single,'g');
hold on;
plot(z_co, 1/exp(1)*ones(num_scans,1),'r');