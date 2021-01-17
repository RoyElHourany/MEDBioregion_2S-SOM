%% Apprentissage
close all; clear all; clc
load('CHLSSTMLD_Data.mat');

U(U<0)=NaN;
U(U>30)=NaN;
U(U==0)=NaN;
U(:,1:12)=log10(U(:,1:12));
 U=U(~any(isnan(U(:,1)),2),:);
Shid_1c=som_data_struct(U);
Shid_1cc=som_normalize(Shid_1c,'var');
 n=50;

DimData=[12 12 12];
for di=1:length(DimData)
    DimBloc(di).Dim=DimData(di);

end

lambda=1;
eta=2000;

radius=[100,1,0.5];
trainlen=[100,2000];
radius2s=[0.5,0.01];
trainlen2s=3000;

[sMap_hid1, sMap_denorm, Result] = learn_2s_som(Shid_1cc.data,n,'lattice','hexa'...
   ,'radius',radius,'trainlen',trainlen,'radius-2s-som',radius2s,'trainlen-2s-som',trainlen2s...
   ,'s2-som','dimdata',DimData,'lambda', lambda,'eta',eta);

% Alpha weights
figure;
Poid=Result.Alpha;
for i=1:3
    subplot(1,3,i), som_cplane(sMap_hid1,Poid(:,i));
    caxis([0 1])
end
% Variable plotting on the 2S-Som map
figure;
for i=1:36
    subplot(3,12,i), som_cplane(sMap_hid1,sMap_denorm.codebook(:,i));
%     caxis([0 1])
end

%% Dendrogram
Z = som_cllinkage(sMap_hid1.codebook,'ward','topol','neighbors');
D = pdist(sMap_hid1.codebook);
leafOrder = optimalleaforder(Z.tree,D);
figure()
dendrogram(Z.tree,10);
%Clustering results on the 2S-Som map
nc=10;
c = cluster(Z.tree,'maxclust',nc);
figure;
h=som_cplane(sMap_hid1,c);
colormap(jet(nc));

%% Projection phase
load sstchl_cycleDATA-MED.mat
U(U<0)=NaN;
U(U>30)=NaN;
U(U==0)=NaN;
U(:,1:12)=log10(U(:,1:12));


U=[U, nan(length(U),12)];
U=som_normalize(U,Shid_1cc);
bmus2=som_bmus(sMap_hid1,U);
bmc = nan*ones(size(bmus2));
for i=1:size(sMap_hid1.codebook,1)
    bmc(bmus2==i)=c(i);
end


load lon_med.mat
load lat_med.mat

figure;
bms=reshape(bmc,1086,392);
bms(800:end,1:126)=nan; bms(1:160,1:75)=nan;

pcolor (lon,lat,bms)
shading flat
axis('equal','tight')
colormap(jet(nc))

