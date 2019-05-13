function cm2_deconv_admm_3D(idx,tau_tv,tau_l1)
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));


save_folder = ['save/',num2str(idx),'/'];
load('P2_10_25_D0508.mat');

if idx == 1
    y = y;
    psfs = psfs;
elseif idx == 2
    y = y2;
    psfs = psfs2;
end

% set tau ranges here
para.tau_tv = (0.0001)*(10^(tau_tv-1));
para.tau_l1 = (0.00001)*(10^(tau_l1-1));
% set other parameters
para.maxiter = 192;
para.display_flag = 0;
para.img_save_period = 32;
para.img_save_path = ['save/',num2str(idx),'/temp/tv_',num2str(tau_tv),'_l1_',num2str(tau_l1)];

disp(['object idx: ',num2str(idx),', tau tv: ',num2str(tau_tv),', tau l1:',num2str(tau_l1)]);
xhat = ADMM_LSI_deconv_3D(y,psfs,para);
xhat = linear_normalize(xhat);
write_mat_to_tif(xhat,[save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.tif']);
% xhat = uint8(255*linear_normalize(xhat));
% for i = 1:size(psfs,3)
%     imwrite(xhat(:,:,i),...
%         [save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'_layer_',num2str(i),'.png']);
% end
% xhat = single(xhat);
% save([save_folder,'tv_',num2str(tau_tv),'_l1_',num2str(tau_l1),'.mat'],'xhat');
end