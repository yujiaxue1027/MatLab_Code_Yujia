function obr(img, varargin)

newFig = 1;
scale = 1;
caption = [];
dock = 0;
map = [];
range = [];
gamma = 1;
n = size(varargin,2);
i = 1;
while i <= n
    p = lower(varargin{i});
    switch p
        case 'nofig'
            newFig = 0;
        case 'noscale'
            scale = 0;
        case 'caption'
            caption = varargin{i+1};
            i = i + 1;
        case 'dock'
            dock = 1;
        case {'colormap','map'}
            map = varargin{i+1};
            i = i + 1;
        case 'range'
            range = varargin{i+1};
            i = i + 1;
        case 'gamma'
            gamma = varargin{i+1};
            i = i + 1;
        otherwise
            error(['Unknown named parameter (' p ').']);
    end
    i = i + 1;
end

if newFig
    f = figure;
else
    f = gcf;
end
if ~isempty(caption)
    set(gcf,'Name',caption);
end
if dock
    set(gcf,'WindowStyle','docked');
%     drawnow
%     plotedit off
end

if ischar(img)
    img = imread(img);
end
if size(img,3) == 3
    if scale
        img = double(img);
        m = [min(img(:)) max(img(:))];
        if diff(m) == 0
            m(2) = m(1) + 1;
        end
        img = (img - m(1)) / (m(2) - m(1));
    end
    if gamma ~= 1
        img = img.^gamma;
    end
    image(img);
else
    if scale && gamma ~= 1
        img = double(img);
        m = [min(img(:)) max(img(:))];
        if diff(m) == 0
            m(2) = m(1) + 1;
        end
        img = (img - m(1)) / (m(2) - m(1));
        img = img.^gamma;
    end
    imagesc(img);
    if isempty(map)
        map = gray(256);
    end
    if ischar(map)
        map = eval([map '(256)']);
    end
    colormap(map);
    if ~isempty(range)
        caxis(range);
    end
end
axis(gca,'image');
% axis(get(f,'Children'),'image');
