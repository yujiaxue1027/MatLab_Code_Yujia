%% initialize stage
handle_horizontal;
handle_vertical;
%% change max and min setting in GUI if needed
horizontal_c = 0.2076;
vertical_c = -0.0610;
handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_c);
handle_horizontal_LPRN.MoveAbsolute(0,1==0);
handle_vertical_UPDN.SetAbsMovePos(0,vertical_c);
handle_vertical_UPDN.MoveAbsolute(0,1==0);

%% initialize camera
vid = videoinput('winvideo', 1, 'RGB32_2592x1944');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
save_path = '0718R10S99/';
start(vid);
tmp = getdata(vid);
imwrite(uint8(tmp(:,:,2)),'initial_img_center.tif');

%% scanning
num_scans = 99;
half_scan_range = 5;
horizontal_scan = linspace(horizontal_c+half_scan_range,horizontal_c-half_scan_range,num_scans);
vertical_scan = linspace(vertical_c-half_scan_range,vertical_c+half_scan_range,num_scans);

order = 1; %% 1: left to right j++; 0: right to left j--
for i = 46:num_scans
    for j = 1:num_scans
        if order == 1
            handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_scan(j));
            handle_horizontal_LPRN.MoveAbsolute(0,1==0);
        else
            handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_scan(num_scans+1-j));
            handle_horizontal_LPRN.MoveAbsolute(0,1==0);
        end
        handle_vertical_UPDN.SetAbsMovePos(0,vertical_scan(i));
        handle_vertical_UPDN.MoveAbsolute(0,1==0);
        pause(2);
        start(vid);
        tmp = getdata(vid);
        if order == 1
            imwrite(uint8(tmp(:,:,2)),...
                [save_path,num2str(i,'%.2d'),...
                '_',num2str(j,'%.2d'),'.tif']);
            pause(2);
        else
            imwrite(uint8(tmp(:,:,2)),...
                [save_path,num2str(i,'%.2d'),...
                '_',num2str(num_scans+1-j,'%.2d'),'.tif']);
            pause(2);
        end
        disp([num2str(i),'__',num2str(j)]);
    end
    if order == 1
        order=0;
    else
        order=1;
    end
end
