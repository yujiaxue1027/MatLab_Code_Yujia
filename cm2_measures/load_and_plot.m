function output = load_and_plot(a,b,c,d,e,f)
clip = @(x, vmin, vmax) max(min(x, vmax), vmin);
filename =  [num2str(a),'_',num2str(b),'_',num2str(c),...
    '_',num2str(d),'_',num2str(e),'_',num2str(f),'.mat'];
load(filename,'est');
figure,imagesc(clip(est,0,1000)),axis image off;colormap gray;truesize;colorbar;
title(filename(1:end-4));
output = est;
end

