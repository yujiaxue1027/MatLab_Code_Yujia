function visualize_tuning(num_tv_taus,num_l1_taus,initial_tv_tau,initial_l1_tau,tv_tau_increase_rate,l1_tau_increase_rate,varargin)
get_tau = @(initial, idx, rate) num2str(round(initial*rate^(idx-1),6));

% num_tv_taus = input('how many tv parameters: ');
% num_l1_taus = input('how many l1 parameters: ');
% initial_tv_tau = input('what is initial tv tau: ');
% initial_l1_tau = input('what is initial l1 tau: ');
% tv_tau_increase_rate = input('what is tv tau increase rate: ');
% l1_tau_increase_rate = input('what is l1 tau increase rate: ');

if length(varargin) == 2
    tv_tau_step = varargin{1};
    l1_tau_step = varargin{2};
else
    tv_tau_step = 1;
    l1_tau_step = 1;
end

disp('loading...');
num_tv = 1 + floor((num_tv_taus-1) / tv_tau_step);
num_l1 = 1 + floor((num_l1_taus-1) / l1_tau_step);

data = zeros(3888, 5184, num_tv, num_l1,'single');
for i = 1:num_tv
    for j = 1:num_l1
        filename = ['tv_',num2str((i-1)*tv_tau_step+1),'_l1_',num2str((j-1)*l1_tau_step+1),'.png'];
        if isfile(filename)
            data(:,:,i,j) = im2single(imread(filename));
        end
    end
end
disp('finish loading.')

f1 = figure('rend','painters','pos',[50 50 1500 900]);
imagesc(data(:,:,1,1)),axis image off; colormap('parula');
disp('select region:')
[x,y] = ginput(2);
imagesc(data(y(1):y(2),x(1):x(2),1,1)),axis image off; colormap(gray(256));
title(['tv tau: ',get_tau(initial_tv_tau,1,tv_tau_increase_rate), ' l1 tau: ',get_tau(initial_l1_tau,1,l1_tau_increase_rate)]);

current_tv_idx=1;
current_l1_idx=1;
while true
    k = waitforbuttonpress;
    if k
        key = double(get(gcf,'CurrentCharacter'));
        switch key
            case 28 %l1 decrease
                if current_l1_idx ~= 1
                    current_l1_idx  = current_l1_idx - 1;
                    imagesc(data(y(1):y(2),x(1):x(2),current_tv_idx,current_l1_idx)),axis image off; colormap(gray(256));
                    title(['tv tau: ',get_tau(initial_tv_tau,(current_tv_idx-1)*tv_tau_step+1,tv_tau_increase_rate), ' l1 tau: ',get_tau(initial_l1_tau,(current_l1_idx-1)*l1_tau_step+1,l1_tau_increase_rate)]);
                end
            case 29 %l1 increase
                if current_l1_idx ~= num_l1
                    current_l1_idx  = current_l1_idx + 1;
                    imagesc(data(y(1):y(2),x(1):x(2),current_tv_idx,current_l1_idx)),axis image off; colormap(gray(256));
                    title(['tv tau: ',get_tau(initial_tv_tau,(current_tv_idx-1)*tv_tau_step+1,tv_tau_increase_rate), ' l1 tau: ',get_tau(initial_l1_tau,(current_l1_idx-1)*l1_tau_step+1,l1_tau_increase_rate)]);
                end
            case 30 %tv decrease
                if current_tv_idx ~= 1
                    current_tv_idx  = current_tv_idx - 1;
                    imagesc(data(y(1):y(2),x(1):x(2),current_tv_idx,current_l1_idx)),axis image off; colormap(gray(256));
                    title(['tv tau: ',get_tau(initial_tv_tau,(current_tv_idx-1)*tv_tau_step+1,tv_tau_increase_rate), ' l1 tau: ',get_tau(initial_l1_tau,(current_l1_idx-1)*l1_tau_step+1,l1_tau_increase_rate)]);
                end
            case 31 %tv increase
                if current_tv_idx ~= num_tv
                    current_tv_idx  = current_tv_idx + 1;
                    imagesc(data(y(1):y(2),x(1):x(2),current_tv_idx,current_l1_idx)),axis image off; colormap(gray(256));
                    title(['tv tau: ',get_tau(initial_tv_tau,(current_tv_idx-1)*tv_tau_step+1,tv_tau_increase_rate), ' l1 tau: ',get_tau(initial_l1_tau,(current_l1_idx-1)*l1_tau_step+1,l1_tau_increase_rate)]);
                end
            otherwise
                close all;
                return
        end
    end
end 
end

