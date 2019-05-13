%% script to display and capture data for SLM deep learning proj

%% training data
NET.addAssembly([pwd, '\Thorlabs.TSI.TLCamera.dll']);
tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;
serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
disp([num2str(serialNumbers.Count), ' camera was discovered.']);
disp('Opening the first camera')
tlCamera = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);

% Set exposure time and gain of the camera.
tlCamera.ExposureTime_us = 1000000;

% Check if the camera supports setting "Gain"
gainRange = tlCamera.GainRange;
if (gainRange.Maximum > 0)
    tlCamera.Gain = 0;
end

% Set the FIFO frame buffer size. Default size is 1.
tlCamera.MaximumNumberOfFramesToQueue = 1;
disp('Starting continuous image acquisition.');
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
tlCamera.Arm;
tlCamera.IssueSoftwareTrigger;
maxPixelIntensity = double(2^tlCamera.BitDepth - 1);


num_imgs = 5;
time1 = 2;
time2 = 0;



img = imread(['D:\Yunzhe\deep learning repository\denseNet\data\MNIST\test_new\black.tif']);
    %for 512 * 512
    %img = padarray(img,[284,704]);
    img = padarray(img,[476,896]);
    img = repmat(img,[1,1,1]);
    fullscreen(img,2);
    pause(time1);
    
    imageFrame = tlCamera.GetPendingFrameOrNull;
 
    capture = uint16(imageFrame.ImageDataUShort1D.ImageData_monoOrBGR)';
    capture = reshape(capture,[1920,1080]);
    capture = capture';
 %   capture = uint8(double(capture)/255);
%     capture = capture(76:1024,538:1463);
%     capture = imresize(capture,[512,512]);
    capture = flip(capture,2);
    capture0 = capture(475:602,895:1022);
    imwrite(capture0,['D:\Yunzhe\deep learning repository\denseNet\data\MNIST\test_new\test\bg.tif'])
    delete(imageFrame);
     pause(time2);


for i = 1:num_imgs
    img = imread(['D:\Yunzhe\deep learning repository\denseNet\data\MNIST\test_new\',num2str(i-1),'.tif']);
    %for 512 * 512
    %img = padarray(img,[284,704]);
    img = padarray(img,[476,896]);
    img = repmat(img,[1,1,1]);
    fullscreen(img,2);
    pause(time1);
    
    imageFrame = tlCamera.GetPendingFrameOrNull;
 
    capture = uint16(imageFrame.ImageDataUShort1D.ImageData_monoOrBGR)';
    capture = reshape(capture,[1920,1080]);
    capture = capture';
%    capture = uint8(double(capture)/255);
%     capture = capture(76:1024,538:1463);
%     capture = imresize(capture,[512,512]);
    capture = flip(capture,2);
    capture = capture(475:602,895:1022);
    capture = capture-capture0;
    imwrite(capture,['D:\Yunzhe\deep learning repository\denseNet\data\MNIST\test_new\test\',num2str(i-1),'.tif'])
    delete(imageFrame);
     pause(time2);
end

% Stop continuous image acquisition
disp('Stopping continuous image acquisition.');
tlCamera.Disarm;
% Release the TLCamera
disp('Releasing the camera');
tlCamera.Dispose;
delete(tlCamera);
% Release the serial numbers
delete(serialNumbers);
% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);
closescreen();

% %% testing data
% 
% 
% NET.addAssembly([pwd, '\Thorlabs.TSI.TLCamera.dll']);
% tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;
% serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
% disp([num2str(serialNumbers.Count), ' camera was discovered.']);
% disp('Opening the first camera')
% tlCamera = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);
% 
% % Set exposure time and gain of the camera.
% tlCamera.ExposureTime_us = 900000;
% 
% % Check if the camera supports setting "Gain"
% gainRange = tlCamera.GainRange;
% if (gainRange.Maximum > 0)
%     tlCamera.Gain = 0;
% end
% 
% % Set the FIFO frame buffer size. Default size is 1.
% tlCamera.MaximumNumberOfFramesToQueue = 1;
% disp('Starting continuous image acquisition.');
% tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
% tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
% tlCamera.Arm;
% tlCamera.IssueSoftwareTrigger;
% maxPixelIntensity = double(2^tlCamera.BitDepth - 1);
% 
% 
% num_imgs = 16;
% time1 = 2;
% time2 = 1;
% for i = 1:num_imgs
%     img = imread(['SISR/test/y/',num2str(i-1),'.tif']);
%     img = padarray(img,[284,704]);
%     img = repmat(img,[1,1,3]);
%     fullscreen(img,2);
%     pause(time1);
%     
%     imageFrame = tlCamera.GetPendingFrameOrNull;
%  
%     capture = uint16(imageFrame.ImageDataUShort1D.ImageData_monoOrBGR)';
%     capture = reshape(capture,[1920,1080]);
%     capture = capture';
%     capture = uint8(double(capture)/255);
% %     capture = capture(76:1024,538:1463);
% %     capture = imresize(capture,[512,512]);
%     capture = flip(capture,1);
%     imwrite(capture,['SISR/test/x/',num2str(i-1),'.tif'])
%     delete(imageFrame);
%     pause(time2);
% end
% % Stop continuous image acquisition
% disp('Stopping continuous image acquisition.');
% tlCamera.Disarm;
% % Release the TLCamera
% disp('Releasing the camera');
% tlCamera.Dispose;
% delete(tlCamera);
% % Release the serial numbers
% delete(serialNumbers);
% % Release the TLCameraSDK.
% tlCameraSDK.Dispose;
% delete(tlCameraSDK);
% closescreen();