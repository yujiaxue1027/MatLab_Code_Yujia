function cm2_deconv_admm_l1tv_0323(idx,tau_tv,tau_l1)
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));
save_folder = ['save/',num2str(idx),'/'];
load('data_0322.mat');
load('default_parameters.mat');

y = ys(:,:,idx);
psf = psfs(:,:,idx);

% set tau ranges here
para.tau_tv = (1e-5)*(2^(tau_tv-1));
para.tau_l1 = (0.01)*(2^(tau_l1-1));
%

disp(['object idx: ',num2str(idx),', tau tv: ',num2str(tau_tv),', tau l1:',num2str(tau_l1)]);
[xhat,~] = ADMM_LSI_deconv_l1tv(y,psf,para);
imwrite(uint8(255*linear_normalize(xhat)),[save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.png']);
% save([save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.mat'],'xhat');
end