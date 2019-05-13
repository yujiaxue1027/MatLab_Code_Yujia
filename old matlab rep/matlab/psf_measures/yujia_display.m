load('pattern1.mat');
pattern = repmat(output./2,[1,1,3]);
pattern_size = size(output,1);
offset = pattern_size;

vid = videoinput('winvideo', 1, 'RGB32_2592x1944');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
num_scans = 31;
[x,y] = meshgrid(linspace(-1,1,num_scans));
displace = 500;
pos_r = round(513 + displace * y);
pos_c = round(641 + displace * x);
% data = zeros(1944,2592,size(pos_r,1),size(pos_r,2));
save_path = '0426/';
for i = 1:size(pos_r,1)
    for j = 1:size(pos_r,2)
        r = pattern_size/2;
        img = zeros(1024+2*offset,1280+2*offset,3);
        img(pos_r(i,j)-r+offset:pos_r(i,j)+r-1+offset,...
            pos_c(i,j)-r+offset:pos_c(i,j)+r-1+offset,:) = pattern; 
        fullscreen(uint8(img(offset+1:offset+1024,offset+1:offset+1280,:)),2);
        pause(2);
        start(vid);
        tmp = getdata(vid);
        imwrite(uint8(tmp(:,:,2)),...
            [save_path,num2str(i,'%.2d'),...
            '_',num2str(j,'%.2d'),'.tif']);
        pause(2);
        disp([num2str(i),'__',num2str(j)]);
    end
    img = zeros(1024,1280,3);
    fullscreen(uint8(img),2);
    pause(2);
    start(vid);
    bg = getdata(vid);
    imwrite(uint8(bg(:,:,2)),[save_path,'bg',num2str(i),'.tif']);
end
img = zeros(1024,1280,3);
fullscreen(uint8(img),2);
pause(2);
start(vid);
bg = getdata(vid);
imwrite(uint8(bg(:,:,2)),[save_path,'bg_final.tif']);
% for i = 1:3
%     for j = 1:3
%         figure,imshow(uint8(data(:,:,i,j)));
%     end
% end
% figure,imshow(uint8(bg));
% save data0330