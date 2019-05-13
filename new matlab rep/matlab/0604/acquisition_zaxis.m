%% zaxis scanning and acquisition
%% initialize stage
handle_zaxis;

%% initialize camera
vid = videoinput('winvideo', 1, 'RGB32_2592x1944');
set(vid,'Timeout',50);
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
save_path = '3d_target/';
start(vid);
while get(vid,'FramesAvailable')<1
      unavailable=1;
end
tmp = getdata(vid);
imwrite(uint8(tmp(:,:,2)),[save_path,'initial_test.tif']);

%% scanning
% num_scans = 10;
z_scan = 4.5:0.1:7.5;
handle_zaxis_UPDN.SetAbsMovePos(0,z_scan(1));
handle_zaxis_UPDN.MoveAbsolute(0,1==0);

for i = 1:length(z_scan)
        handle_zaxis_UPDN.SetAbsMovePos(0,z_scan(i));
        handle_zaxis_UPDN.MoveAbsolute(0,1==0);
        pause(3);
        start(vid);
        while get(vid,'FramesAvailable')<1
            unavailable=1;
        end
        tmp = getdata(vid);

        imwrite(uint8(tmp(:,:,2)),...
            [save_path,num2str(i,'%.2d'),'.tif']);
        pause(2);
        disp([num2str(i)]);
end
