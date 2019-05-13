function cm2_deconv_admm_l1tv(idx,tau_tv,tau_l1)
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));
data_file_name = ['y',num2str(idx),'_s.mat'];
default_para_file_name = 'default_parameters.mat';
save_folder = ['save/',num2str(idx),'/'];

load(data_file_name);
load(default_para_file_name);
load('new_psfs.mat');
psf = npsfs(:,:,idx);
psf(psf<=1e-5) = 0;
psf = psf./sum(psf(:)) .*9;
psf = 9.*modify_psf(psf,[1,1,1;1,1,1;1,1,1]);


% set tau ranges here
para.tau_tv = (1e-5)*(3.6^(tau_tv-1));
para.tau_l1 = (0.001)*(2.3^(tau_l1-1));
%

disp(['object idx: ',num2str(idx),', tau tv: ',num2str(tau_tv),', tau l1:',num2str(tau_l1)]);
[xhat,~] = ADMM_LSI_deconv_l1tv(y,psf,para);
imwrite(uint8(255*linear_normalize(xhat)),[save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.png']);
% save([save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.mat'],'xhat');
end