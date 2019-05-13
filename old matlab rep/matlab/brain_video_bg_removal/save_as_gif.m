function save_as_gif(video,gif_name,delay)
filename = [gif_name,'.gif'];
if length(size(video))==3
    video = reshape(video,[size(video,1),size(video,2),1,size(video,3)]);
    video = repmat(video,[1,1,3,1]);
end
video = abs(video);
video = uint8(video./max(video(:)).*255);
for idx = 1:size(video,3)
    [A,map] = rgb2ind(video(:,:,:,idx),256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end
end

