%%%%%%%%%%%%%%%%%%%%%%%% LFM matrix factorization %%%%%%%%%%%%%%%%%%%%%%%%

%% Input
Input.LFM_folder='lfm-movie-address on HDD';
Input.psf_filename_ballistic='psf-file-address on HDD';
Input.output_path='output-folder-address on HDD';
Input.x_offset=1483.800000;
Input.y_offset=1379.100000;
Input.dx=21.83;
Input.step=1;
Input.temporal_iterations=100000;
Input.spatial_iterations=50000;
Input.bg_iter=2;
Input.noise_level=100;
Input.thres=0.0160;
Input.rectify=1;
Input.Junk_size=100;
Input.bg_sub=1;
Input.sample=1;

psf_ballistic=matfile(Input.psf_filename_ballistic);

%% Compute bg components via rank-1-factorization
tic
if Input.bg_sub==1
    [output.bg_temporal,output.bg_spatial]=par_rank_1_factorization(Input.LFM_folder,Input.step, Input.bg_iter);
else
    output.bg_temporal=[];
    output.bg_spatial=[];
end
toc
%% compute standard-deviation image (std. image)
tic
if Input.rectify==1
    [output.std_image,~]=par_compute_std_image(Input.LFM_folder,Input.step,output.bg_temporal,output.bg_spatial,[], Input.x_offset,Input.y_offset,Input.dx,psf_ballistic.Nnum);
else
    [output.std_image,~]=par_compute_std_image(Input.LFM_folder,Input.step,output.bg_temporal,output.bg_spatial);
end
if (Input.bg_sub==1)&&(Input.rectify==1)
    output.bg_spatial =  ImageRect(output.bg_spatial, Input.x_offset, Input.y_offset, Input.dx, psf_ballistic.Nnum,0);
end
toc
%% filter std. image
Nnum=psf_ballistic.Nnum;

IMG=output.std_image;
for k=1:size(IMG,1)/Nnum
    for j=1:size(IMG,2)/Nnum
        if (k==1)||(j==1)||(k==size(IMG,1)/Nnum)...
                ||(j==size(IMG,2)/Nnum)
            filt= output.bg_spatial((k-1)*Nnum+1:k*Nnum,...
                (j-1)*Nnum+1:j*Nnum)/...
                sqrt(sum(sum(output.bg_spatial((k-1)*Nnum+1:k*Nnum,...
                (j-1)*Nnum+1:j*Nnum).^2)));
            filtered_3((k-1)*Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum)=IMG((k-1)...
                *Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum)-sum(sum(IMG((k-1)...
                *Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum).*filt))*filt;
        else
            filt= output.bg_spatial(((k-1)-1)*Nnum+1:(k+1)...
                *Nnum,((j-1)-1)*Nnum+1:(j+1)*Nnum);
            filt=filt/sqrt(sum(sum(filt.^2)));
            filt1= output.bg_spatial((k-1)*Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum)/sqrt(sum(sum(output.bg_spatial((k-1)...
                *Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum).^2)));
            filtered_3((k-1)*Nnum+1:k*Nnum,(j-1)*...
                Nnum+1:j*Nnum)=IMG((k-1)...
                *Nnum+1:k*Nnum,(j-1)...
                *Nnum+1:j*Nnum)-sum(sum(IMG(((k-1)-1)...
                *Nnum+1:(k+1)*Nnum,((j-1)-1)...
                *Nnum+1:(j+1)*Nnum).*filt))*filt1/9;
        end
    end
    disp(k);
end
filtered_3(filtered_3<0)=0;
filtered_3=filtered_3-Input.noise_level;
filtered_3(filtered_3<0)=0;
clear IMG;
clear filt1;
clear filt;

%% reconstruct std. image
tic
psf_ballistic=load(Input.psf_filename_ballistic);
toc
tic
infile=struct;
infile.LFmovie=filtered_3;
options.gpu_ids=5;
output.reconstructed_std= reconstruction_gpu(infile, psf_ballistic, '/home/', '/tmp', 1, 8, options,0,0,1);
gpuDevice([]);
toc

%% interpolate artifact plane in reconstructed std. image

output.reconstructed_std(:,:,62)=(output.reconstructed_std(:,:,61)+output.reconstructed_std(:,:,63))/2;

%% band-pass-filter reconstructed std. image
psf_ballistic=matfile(Input.psf_filename_ballistic);
Hsize = size(psf_ballistic, 'H');
m=[size(output.std_image,1),size(output.std_image,2),Hsize(5)];
cellSize = 15;

I=zeros(size(output.reconstructed_std)+[0 0 2*cellSize],'single');
I(:,:,cellSize+1:cellSize+Hsize(5))=single(output.reconstructed_std);
for k=0:cellSize-1
    I(:,:,cellSize-k)=I(:,:,cellSize+1-k)*0.96;
    I(:,:,cellSize+Hsize(5)+k)=I(:,:,cellSize+Hsize(5)-1+k)*0.96;
    
end
cellSize = 51;

Ifiltered = I/max(I(:));

Ifiltered=band_pass_filter(Ifiltered, cellSize, 10, 4, 3);
filtered_Image = Ifiltered;
bord=1;

cellSize=15;
segm3=zeros(size(filtered_Image)-[0 0 2*cellSize]);
segm3(bord:size(filtered_Image,1)-bord,bord:size(filtered_Image,2)-bord,:)=filtered_Image(bord:size(filtered_Image,1)-bord,bord:size(filtered_Image,2)-bord,cellSize+1:cellSize+Hsize(5));
%% threshold
segm=segm/max(segm(:));
segm=segm-Input.thres;
segm(segm<0)=0;
%% extract neuron centers
output.centers=[];
B=reshape(segm,[],1);

beads=bwconncomp(segm);
for k=1:beads.NumObjects
    qu=B(beads.PixelIdxList{1,k});
    q=sum(B(beads.PixelIdxList{1,k}));
    [a,b,c]=ind2sub([m(1) m(2) Hsize(5)],beads.PixelIdxList{1,k});
    output.centers(k,:)=([a,b,c]'*qu/q)';
end

%% Initiate forward_model
psf_ballistic=load(Input.psf_filename_ballistic);

tic
forward_model_indices=cell(1,size(output.centers,1));
forward_model_values=forward_model_indices;
N=0;

rr=4;
r=8;

BW=[];
BWW=[];

W=zeros(2*r,2*r,2*rr);
for ii=1:2*r
    for jj=1:2*r
        for kk=1:2*rr
            if  ((ii-((2*r-1)/2+1))^2/r^2+(jj-((2*r-1)/2+1))^2/r^2+(kk-((2*rr-1)/2+1))^2/rr^2)<=1
                W(ii,jj,kk)=1;
            end
        end
    end
end
BW=bwconncomp(W);
[BWW(:,1) BWW(:,2) BWW(:,3)]=ind2sub([2*r,2*r,2*rr],BW.PixelIdxList{1,1});

for k=1:size(output.centers,1)
    B=[];
    for j=1:size(BWW,1)
        bbb=round(BWW(j,:)-[((2*r-1)/2+1)*[1 1] ((2*rr-1)/2+1)]+output.centers(k,:));
        if (bbb(1)<=m(1))&&(bbb(1)>0)&&(bbb(2)<=m(2))&&(bbb(2)>0)&&(bbb(3)<=Hsize(5))&&(bbb(3)>0)
            B=[B' bbb']';
        end
    end
    gt{1,1}=B;
    Q=forwardproject(gt,psf_ballistic,size(output.std_image));
    forward_model_indices{1,k}=find(Q);
    forward_model_values{1,k}=Q(forward_model_indices{1,k});
    N=N+length(forward_model_values{1,k});
    disp(k);
end
I=zeros(N,1);
J=I;
S=I;
jj=0;
for k=1:size(forward_model_indices,2)
    J(jj+1:jj+size(forward_model_values{1,k},2))= forward_model_indices{1,k};
    I(jj+1:jj+size(forward_model_values{1,k},2))=k*ones(size(forward_model_values{1,k}));
    S(jj+1:jj+size(forward_model_values{1,k},2))=forward_model_values{1,k};
    jj=jj+size(forward_model_values{1,k},2);
    disp(k)
end
output.forward_model=sparse(I,J,S,size(output.centers,1),size(output.std_image,1)*size(output.std_image,2));
toc;
%% Non negative matrix factorization

II=[];JJ=[];
tic
for neuron=1:size(output.forward_model,1)
    img=reshape(output.forward_model(neuron,:),size(output.std_image));
    img_=zeros(size(output.std_image)/Nnum);
    for k=1:size(output.std_image,1)/Nnum
        for j=1:size(output.std_image,2)/Nnum
            img_(k,j)=mean(mean(img((k-1)*Nnum+1:k*Nnum,(j-1)*Nnum+1:j*Nnum)));
        end
    end
    img_=img_/max(img_(:));
    img_(img_<0.0001)=0;
    [I_,J_,~]=find(img_);
    I=[];J=[];
    for k=1:length(I_)
        s=(I_(k)-1)*Nnum+1:I_(k)*Nnum;
        for l=1:Nnum
            I=[I' (ones(Nnum,1)*s(l))']';
            J=[J' ((J_(k)-1)*Nnum+1:J_(k)*Nnum)]';
        end
    end
    
    II=[II' (ones(size(I))*neuron)']';
    JJ=[JJ' sub2ind(size(img),I,J)']';
    disp(neuron);
end
toc
template=sparse(II,JJ,ones(size(II)),size(output.forward_model,1),size(output.std_image,1)*size(output.std_image,2));
%%
tic
sensor_movie=read_sensor_movie(Input.LFM_folder,Input.x_offset,Input.y_offset,Input.dx,Nnum,Input.sample,Input.rectify);
toc

tic
output.timeseries=output.forward_model*sensor_movie;
o=output.timeseries;
B=output.forward_model*output.forward_model';
for k=1:Input.temporal_iterations
output.timeseries=output.timeseries.*o./(B*output.timeseries);
end
toc
output.timeseries_=output.timeseries;

for iter=1:3
    tic
    output.forward_model_=update_spatial_component(output.timeseries_, sensor_movie, template, optz);
    toc
    if Input.bg_sub==1
        bg_spatial_=(squeeze(max(template,[],1))==0);
        bg_spatial_=bg_spatial_.*(sensor_movie*output.timeseries_(end,:)')'/(output.timeseries_(end,:)*output.timeseries_(end,:)');
        bg_spatial_=bg_spatial_+squeeze(output.forward_model_(end,:));
        [I,J,S]=find(output.forward_model_);
        f=find(I~=size(output.timeseries_,1));
        I=I(f);
        J=J(f);
        S=S(f);
        output.forward_model_=sparse(I,J,S,size(output.centers,1),size(output.forward_model_,2));
    end
    tic
    output.timeseries_=output.forward_model_*sensor_movie;
    o=output.timeseries_;
    B=output.forward_model_*output.forward_model_';
    for k=1:Input.temporal_iterations
        output.timeseries_=output.timeseries_.*o./(B*output.timeseries_);
    end
    toc
    disp(iter);
end

%% save output

save(Input.output_path,'-struct','output','-v7.3');

