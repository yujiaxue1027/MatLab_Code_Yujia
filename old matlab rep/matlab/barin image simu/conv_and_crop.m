%%conv_and_crop
%% generate gt
part1 = 1*sin(2*pi*[0.01:0.01:1]);
part2 = 0.5*sin(2*pi*[0.01:0.01:1]);
gt = [part1';part2'];
% gt = padarray(gt,[0,length(gt)/4]);
figure,plot(gt),title('gt');

%% generate psf
spsf = [0 0  0.1 0.2 0.5 0.2 0.1 0 0];
spsf = spsf./sum(spsf);
spsf = spsf';
spacing = 30;
copy = [3,1];
spsf = padarray(spsf,[round((spacing-length(spsf))/2),...
    0]);
psf = repmat(spsf,copy);
figure,plot(psf),title('psf');

%% conv and crop
output = conv(gt,psf,'same');
output = output(26:175);
fulloutput = conv(gt,psf,'full');
figure, plot(output),title('conv with crop');
figure,plot(fulloutput),title('conv without crop');

%% conv with convmtx
convmat = convmtx(psf,length(gt));
cropmat = eye(150);
cropmat = padarray(cropmat,[0,71]);
cvcp = cropmat*convmat;
output2 = cvcp*gt;
figure,plot(output2),title('cvcp result');

%% deconv
alpha = 0.01;
result = inv(cvcp'*cvcp + alpha.*eye(200,200))*cvcp'*output;
figure,plot(result),title(['tikhonov alpha = ',num2str(alpha)]);
hold on,plot(gt,'r');

