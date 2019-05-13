function [bg_temporal,bg_spatial,par_C]=par_rank_1_factorization(par_C,indir,step, max_iter, x_offset,y_offset,dx,Nnum)

if nargin==2
    step=1;
end


%%
if exist(indir, 'dir')
    infiles_struct = dir(fullfile(indir, '/*.tif*'));
    [~, order] = sort({infiles_struct(:).name});
    infiles_struct = infiles_struct(order);
else
    disp('indir does not exist');
    return
end

disp('Input files:');
for i=1:size(infiles_struct)
    disp(infiles_struct(i));
end

%%
FileTif=fullfile(indir, infiles_struct(1).name);
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
bg_spatial=ones(nImage,mImage,NumberImages,'double');

if isempty(par_C)
par_C=parpool;
end

infiles_struct = infiles_struct(1:step:end);
N=par_C.NumWorkers;
%%
bg_temporal_par=zeros([par_C.NumWorkers,ceil(length(infiles_struct)/par_C.NumWorkers)]);

bg_spatial=bg_spatial/sqrt(sum(sum(bg_spatial.^2)));
parfor worker=1:par_C.NumWorkers
    bg_spatial_par(worker,:,:)=zeros(size(bg_spatial));
    for i=1:length(worker:N:length(infiles_struct))
        img_rect = double(imread(fullfile(indir, infiles_struct(worker+(i-1)*N).name), 'tiff'));
        parvar=zeros(size(squeeze(bg_temporal_par(worker,:))));
        parvar(i)=sum(sum(bg_spatial.*img_rect));
        bg_temporal_par(worker,:)=bg_temporal_par(worker,:)+parvar;
        bg_spatial_par(worker,:,:)=squeeze(bg_spatial_par(worker,:,:))+img_rect*parvar(i);
    end
    btn(worker)=sum(sum(bg_temporal_par(worker,:).^2));
end
nrm=sqrt(sum(btn));
bg_spatial=squeeze(sum(bg_spatial_par,1))/nrm;
for worker=1:par_C.NumWorkers
    bg_temporal(worker:N:length(infiles_struct))=bg_temporal_par(worker,1:length(worker:N:length(infiles_struct)))';
end
sm_bg_t=smooth(bg_temporal,200);

opt.Algorithm='interior-point';

h=@(x) sub_area_fit(x,sm_bg_t);
xx=fminsearch(h,(min(sm_bg_t)),opt);
r1=fit([1:length(sm_bg_t)]',sm_bg_t-xx,'exp1');
for nn=1:N
    bleach_function_par(nn,1:length(nn:N:length(infiles_struct)))=exp(r1.b*[nn:N:length(infiles_struct)]);
end

disp('anti-bleach');
%%

for iter=1:max_iter
    bg_temporal_par=zeros([par_C.NumWorkers,ceil(length(infiles_struct)/par_C.NumWorkers)]);
    bg_spatial=bg_spatial/sqrt(sum(sum(bg_spatial.^2)));
    parfor worker=1:par_C.NumWorkers
        bg_spatial_par(worker,:,:)=zeros(size(bg_spatial));
        for i=1:length(worker:N:length(infiles_struct))
            img_rect = double(imread(fullfile(indir, infiles_struct(worker+(i-1)*N).name), 'tiff'));
            parvar=zeros(size(squeeze(bg_temporal_par(worker,:))));
            parvar(i)=sum(sum(bg_spatial.*((img_rect-xx)/bleach_function_par(worker,i)+xx)));
            bg_temporal_par(worker,:)=bg_temporal_par(worker,:)+parvar
            bg_spatial_par(worker,:,:)=squeeze(bg_spatial_par(worker,:,:))+((img_rect-xx)/bleach_function_par(worker,i)+xx)*parvar(i);
        end
        btn(worker)=sum(sum(bg_temporal_par(worker,:).^2));               
    end
    nrm=sqrt(sum(btn));
    bg_spatial=squeeze(sum(bg_spatial_par,1))/nrm;
    disp(['finished iteration:' num2str(iter)]);
end

for worker=1:par_C.NumWorkers
    bg_temporal(worker:N:length(infiles_struct))=bg_temporal_par(worker,1:length(worker:N:length(infiles_struct)))';
end

if nargin>4    
    bg_spatial =  ImageRect(bg_spatial, x_offset, y_offset, dx, Nnum,0);    
end

bg_temporal=bg_temporal/nrm;

end
