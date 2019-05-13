%% initialize stage
handle_horizontal;
handle_vertical;
%% change max and min setting in GUI if needed
horizontal_c = -1.8;
vertical_c = 0.9;
handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_c);
handle_horizontal_LPRN.MoveAbsolute(0,1==0);
handle_vertical_UPDN.SetAbsMovePos(0,vertical_c);
handle_vertical_UPDN.MoveAbsolute(0,1==0);

%% initialize camera
vid = videoinput('winvideo', 1, 'RGB32_2592x1944');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
save_path = '0628R8HC/';
start(vid);
tmp = getdata(vid);
imwrite(uint8(tmp(:,:,2)),'initial_img_center.tif');

%% generate scanning pattern
num_scans = 51; % 51
half_scan_range = 4; %5
% horizontal_scan = linspace(horizontal_c+half_scan_range,horizontal_c-half_scan_range,num_scans);
% vertical_scan = linspace(vertical_c-half_scan_range,vertical_c+half_scan_range,num_scans);
num_scans_v = round(num_scans/sqrt(3));
if mod(num_scans_v, 2) == 0
    num_scans_v = num_scans_v +1;
end
unit_h = half_scan_range*2/(num_scans-1);
unit_v = unit_h*sqrt(3);
[x1,y1] = meshgrid(unit_h*[(num_scans-1)/2:-1:-1*(num_scans-1)/2],...
                    unit_v*[-1*(num_scans_v-1)/2:1:(num_scans_v-1)/2]);
x2 = x1 + 0.5*unit_h;
y2 = y1 + 0.5*unit_v;
x = zeros(size(x1,1)*2, size(x1,2));
y = zeros(size(x1,1)*2, size(x1,2));
x(1:2:end,:) = x1;
x(2:2:end,:) = x2;
y(1:2:end,:) = y1;
y(2:2:end,:) = y2;
horizontal_scan = x + horizontal_c;
vertical_scan = y + vertical_c;



%% acquisition
order = 0; %% 1: left to right j++; 0: right to left j--
for i = 22:size(horizontal_scan, 1)
    for j = 1:size(horizontal_scan ,2)
        if order == 1
            handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_scan(i,j));
            handle_horizontal_LPRN.MoveAbsolute(0,1==0);
        else
            handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_scan(i,size(horizontal_scan,2)+1-j));
            handle_horizontal_LPRN.MoveAbsolute(0,1==0);
        end
        handle_vertical_UPDN.SetAbsMovePos(0,vertical_scan(i,j));
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
                '_',num2str(size(horizontal_scan,2)+1-j,'%.2d'),'.tif']);
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
