for index = [18]%[1,2,5,6,8,11,12,15,16,17] %1:20
for delta = [10,20,40,80,160,320] %[10,15,20,30,40,60,80,120,160,240,320]
%% input params
br = 10;
crop = [200,400];
names = {'city-1','city-2','city-3','city-4',...
         'face-1','face-2','face-3','forest-1','forest-2','forest-3',...
         'nature-1','nature-2','nature-3','nature-4',...
         'sign-1','sign-2','sign-3','sign-4','sign-5','sign-6'};
name = names{index};
     
%% init input
speed = br/delta;
bs = [1 1]*(2*br+1);
maxR = min(bs)/2;
og = double(imread(['data/' name '.jpg']))/255;
og = 0.299*og(:,:,1) + 0.587*og(:,:,2) + 0.114*og(:,:,3);
og = og(floor((size(og,1)-crop(1))/2)+(1:crop(1)),floor((size(og,2)-crop(2))/2)+(1:crop(2)));
so = [size(og,1),size(og,2)];
col = size(og,3);
x0 = floor((so(2)-delta-maxR)/2);
og = padarray(og,[3*br,3*br]);

B = zeros([prod(bs),delta+1]);
for i = 0:delta
    B(:,i+1) = reshape(cylinderblur(2*i*speed+1,[0,0],bs),[],1);
end
Bc = zeros([prod(bs),so+6*br]);
Bc((prod(bs)+1)/2,:,:) = 1;
Bc(:,:,3*br+x0:3*br+x0+delta) = repmat(reshape(B,[],1,delta+1),[1,size(Bc,2),1]);
Bc(:,:,3*br+x0+delta+1:end) = repmat(B(:,end),[1,size(Bc,2),so(2)-x0-delta+3*br]);
Br = zeros([prod(bs),so+4*br]);
g = zeros([so+4*br,col]);
i = 0;
for x = 1:bs(2)
    for y = 1:bs(1)
        i = i + 1;
        Br(i,:,:) = Bc(i,bs(1)-y+1:end-y+1,bs(2)-x+1:end-x+1);
        g = g + repmat(squeeze(Br(i,:,:)),[1,1,col]).*og(bs(1)-y+1:end-y+1,bs(2)-x+1:end-x+1,:);
    end
end
Nr = squeeze(sum(Br));
Br = Br./repmat(reshape(Nr,[1,size(Nr)]),[prod(bs),1,1]);

data.og = og(br+1:end-br,br+1:end-br,:);
data.g = g;
data.g2 = g./repmat(Nr,[1,1,col]);
data.sPsf = bs;
data.Dc = reshape(Bc(:,br+1:end-br,br+1:end-br),prod(bs),[]);
data.Dr = reshape(Br,prod(bs),[]);
data.hc = genBase(data.Dc,min(100,prod(data.sPsf)));
data.hr = genBase(data.Dr,min(100,prod(data.sPsf)));
data.hc = reshape(data.hc,data.sPsf(1),data.sPsf(2),[]);
data.hr = reshape(data.hr,data.sPsf(1),data.sPsf(2),[]);

%% method params
par.nBase = 10;           % number of base PSFs
par.alpha = 1;
par.mu = 1e5;
par.beta = sqrt(par.mu*par.alpha);
par.gamma = par.beta;
par.verbose = 1;         % verbosity (0 - none, 1 - normal, 2 - thorough)

if par.verbose > 0
    obr(data.g2(3*br+1:end-3*br,3*br+1:end-3*br,:),'Caption',['Input - ' num2str(delta)],'dock','range',[0 1]);
end

%% run all methods
for k = [10]
    par.nBase=k;
    fprintf('idx %i, delta %i, K=%i, ',index,delta,k);
for m = 1:4
    switch m
        case 1
            if isfield(data,'a')
                data = rmfield(data,'a');
            end
            data.D2 = data.Dr;
            data.h2 = data.hr;
            fprintf('RWD:\n');
            result = initRWD(data, par);mtd='RWD';
        case 2
            if isfield(data,'a')
                data = rmfield(data,'a');
            end
            data.D = data.Dc;
            data.h = data.hc;
            fprintf('CWD:\n');
            result = initCWD(data, par);mtd='CWD';
        case 3
            if isfield(data,'a')
                data = rmfield(data,'a');
            end
            data.D2 = flip(data.Dc,1);
            data.h2 = rot90(data.hc,2);
            fprintf('RWD-C:\n');
            result = initRWD(data, par);mtd='RWD-C';
        case 4
            i = round(linspace(1,size(data.g,2),par.nBase));
            data.h2 = reshape(flip(data.Dc(:,i*size(data.g,1)),1),[data.sPsf,par.nBase]);
            data.a = zeros(par.nBase,size(data.Dc,2));
            p = 1;
            for j = 1:length(i)-1
                ii = repmat(linspace(0,1,i(j+1)-i(j)+1),[size(data.g,1),1]);
                ii(:,end) = [];
                data.a(j,p:p+numel(ii)-1) = 1-ii(:);
                data.a(j+1,p:p+numel(ii)-1) = ii(:);
                p = p+numel(ii);
            end
            data.a(end,p:end) = 1;
            fprintf('RWD-BI:\n');
            result = initRWD(data, par);mtd='BI-R';
        case 5
            i = round(linspace(1,size(data.g,2),par.nBase));
            data.h = reshape(data.Dc(:,i*size(data.g,1)),[data.sPsf,par.nBase]);
            data.a = zeros(par.nBase,size(data.Dc,2));
            p = 1;
            for j = 1:length(i)-1
                ii = repmat(linspace(0,1,i(j+1)-i(j)+1),[size(data.g,1),1]);
                ii(:,end) = [];
                data.a(j,p:p+numel(ii)-1) = 1-ii(:);
                data.a(j+1,p:p+numel(ii)-1) = ii(:);
                p = p+numel(ii);
            end
            data.a(end,p:end) = 1;
            fprintf('CWD-BI:\n');
            result = initCWD(data, par);mtd='BI-C';
    end
    %% iterations
    for i = 1:500
        fprintf('Iteration %i: ', length(result.t)+1);
        t = tic;
        switch m
            case {1,3,4}
                result = iterX_RWD(result);
                result = iterU_RWD(result, par);
            case {2,5}
                result = iterX_CWD(result);
                result = iterU_CWD(result, par);
        end
        t = toc(t);
        result.eu = [result.eu; err2orig(result.u,data,3,[x0-br-1,delta+br+2])];
        result.t = [result.t; t];
        fprintf('%.3fdB, %.3fs\n', result.eu(end), t);
    end
    %% view
    tmp = result.u(1:size(data.g,1),1:size(data.g,2),:);
    tmp = tmp(3*br+1:end-3*br,3*br+1:end-3*br,:);
    tmp(tmp < 0) = 0;
    tmp(tmp > 1) = 1;
    if par.verbose > 0
        obr(tmp,'Caption',[mtd '(' num2str(par.nBase) '_' num2str(log10(par.mu)) '_' num2str(length(result.t)) ')'],'dock','range',[0 1]);
        drawnow
        if par.verbose > 1
            imwrite(tmp,sprintf('results/%s-res-%02i_%03i_%s-new.png',name,par.nBase,delta,mtd));
        end
    end
    clear tmp;
    err{par.nBase}(index,:,m,delta/5) = result.eu;
end
end
% f = figure;plot(squeeze(err{par.nBase}(index,:,:,delta/5)));grid on
% % f = figure;plot([5:21,23],squeeze(err(end,:,[5:21,23])),'-o');grid on
% f.Name = ['Graph - delta=' num2str(delta)];
% % f.WindowStyle = 'docked';
end
end

%% final graph
index = [1,2,5,6,8,11,12,15,16,17];
d = [10,20,40,80,160,320]; %[10,15,20,30,40,60,80,120,160,240,320];
f = figure;
p = plot(br./d,squeeze(mean(err{par.nBase}(index,end,2,d/5),1))','-','LineWidth',3);p.DisplayName = 'CWD - SVD'; hold on
p = plot(br./d,squeeze(mean(err{par.nBase}(index,end,1,d/5),1))','--','LineWidth',3);p.DisplayName = 'RWD - SVD';
p = plot(br./d,squeeze(mean(err{par.nBase}(index,end,3,d/5),1))','-.','LineWidth',3);p.DisplayName = 'RWD - SVD-C';
p = plot(br./d,squeeze(mean(err{par.nBase}(index,end,4,d/5),1))',':','LineWidth',3);p.DisplayName = 'RWD - BI';
grid on
f.Name = ['Graph - ' num2str(par.nBase) ' PSFs'];
a = gca;
a.XScale = 'log';
a.XTick = fliplr(br./[10,20,40,80,160,320]);
a.XTickLabel = fliplr({'1','1/2','1/4','1/8','1/16','1/32'});
a.GridLineStyle = '-';
a.FontSize = 14;
a.XLabel.String = 'Speed of PSF change';
a.YLabel.String = 'PSNR (dB)';
a.XLim = [br/d(end),1];
a.YLim = [15,40];
legend('show');
