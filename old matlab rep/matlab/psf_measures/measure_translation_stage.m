Tian_Lab_Motor;
h_c = -1.1355;
h2_c = -2.7069;
vid = videoinput('winvideo', 1, 'RGB32_2592x1944');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
num_scans = 31;
h_scan = linspace(h_c-5,h_c+5,num_scans);
h2_scan = linspace(h2_c-5,h2_c+5,num_scans);
save_path = '0525/';
for i = 1:num_scans
    for j = 1:num_scans
        h.SetAbsMovePos(0,h_scan(i));
        h.MoveAbsolute(0,1==0);
        h2.SetAbsMovePos(0,h2_scan(j));
        h2.MoveAbsolute(0,1==0);
        pause(2);
        start(vid);
        tmp = getdata(vid);
        imwrite(uint8(tmp(:,:,2)),...
            [save_path,num2str(i,'%.2d'),...
            '_',num2str(j,'%.2d'),'.tif']);
        pause(2);
        disp([num2str(i),'__',num2str(j)]);
    end
end
