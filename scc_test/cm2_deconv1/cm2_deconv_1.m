%% cm2 deconv1 try different parameters for gfp sample
load processed_y_and_selected_psfs.mat
load default_para.mat

y_gfp = ys(:,:,1);
psf_gfp = npsfs(:,:,1);

% tau vary from 1e-6 to 1e0
% clip_max vary from 0.1 to 0.6
idx = 1;
para2.display_flag = false;
while true
    tau_tmp = 10^(rand(1)*6-6);
    clip_max_tmp = rand(1)*0.5 + 0.1;
    para2.tau = tau_tmp;
    para2.clip_max = clip_max_tmp;
    [xhat,~] = ADMM_LSI_deconv_l1(y_gfp,psf_gfp,para2);
    img = uint8(255*linear_normalize(clip(xhat,0,10000)));
    imwrite(img,['save/S1_',num2str(idx),'.png']);
    save(['save/S1_',num2str(idx),'.mat'], 'xhat','para2');
    idx = idx + 1;
end
