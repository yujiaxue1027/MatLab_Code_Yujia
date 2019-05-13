function slm_test(img)

%% initialization of the camera
vid = videoinput('ni', 1, 'img1'); %Shuai's camera
%vid = videoinput('ni', 2, 'img0'); % Ayan's camera
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;


%% Load train images
idxx=1;
sz_org=256;
sz_resz=160;
immat=zeros(sz_org,sz_org,idxx,'uint8');

sz_blk=1;

img=imresize(img,[sz_org,sz_org]);
img=uint8(make_block_image(img,sz_blk));

immat(:,:,1)=img;



%% begin capture

Im_stack=[];
strxx=353;
stryy=485; %these two parameters need to be calibrated
threshhold=255;

if threshhold~=255
    immat=uint8(round(double(immat)*threshhold/255));
end


shift=120;
backg=0;
for numimg=1:idxx
    
    img=backg*ones(768,1024,'uint8');
    imgtemp=immat(:,:,numimg);
    imgtemp=imresize(imgtemp,[sz_resz,sz_resz]);
    img(strxx-sz_resz/2+1:strxx+sz_resz/2,stryy-sz_resz/2+1:stryy+sz_resz/2)=imgtemp;
    img=repmat(img,[1,1,3]);
    fullscreen(img,1); %upload image
    pause(0.5);
    start(vid);
    Im= getdata(vid); %capture image
    stop(vid);
    Im_stack=cat(3,Im_stack,Im(:,257-shift:end-shift));%(:,257-32:end-32)
end
closescreen()

%%
delete(vid);



for indexp=1:idxx 
%     close all
    figure;
    modval=mod(indexp,idxx);
    if modval==0
        modval=idxx;
    end
    subplot(1,2,1);
    imagesc(immat(:,:,modval));
    colorbar; axis image
    subplot(1,2,2);
    imagesc(Im_stack(:,:,indexp),2);
    colorbar; axis image;
    
end


