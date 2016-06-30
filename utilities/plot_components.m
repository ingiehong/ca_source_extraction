%% Plot spatial components
figure;
unit=zeros(d1,d2);
for i=1:size(A_or,2)
    unit=unit + reshape(A(:,i), d1, d2);
    imagesc(unit)
    pause(0.1)
end

%% Compute correlation matrix between spatial components
corr_matrix = CreateCrossCorrMatrix(C_or');
figure
%imshow(corr_matrix,[0.9 1]);
imagesc(corr_matrix, [0.9 1]);
colormap(jet);
axis on;

%% Cluster spatial components

%# remove diagonal elements
corr_matrix(eye(size(corr_matrix))~=0)=0;
%# and convert to a vector (as pdist)
dissimilarity = 1 - corr_matrix(find(corr_matrix))';

%# decide on a cutoff
%# remember that 0.4 corresponds to corr of 0.6!
cutoff = 0.1; 

%# perform complete linkage clustering
Z = linkage(dissimilarity, 'single', 'correlation'); %,'complete');

%# group the data into clusters
%# (cutoff is at a correlation of 0.5)
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');

figure;
dendrogram(Z,0,'colorthreshold',cutoff)

%% Plot spatial components in same cluster

figure;
scrsz = get(groot,'ScreenSize');
set(gcf, 'Position',[50 50 scrsz(3)-100 scrsz(4)*2/3-150])
pan = panel();
pan.pack('h', {1/3 []} )

for i=1:length(unique(groups))
    groupi = find(groups==i);
    pan(1).select()
    unit=zeros(d1,d2);
    for j=1:size(groupi,1)
        unit=unit + reshape(A(:,groupi(j)), d1, d2);
        imagesc(unit)
        axis image; 
    end
    title(['Spatial components: ' num2str(groupi')]);
    pan(2).select()
    plot(C_df(groupi,:)');
    pause
end


%% Regroup ROIs based on clustering
[A_gr,C_gr,S_gr,P_gr,srt] = group_ROIs(A,C,S,P,groups);

%% Compute correlation matrix between spatial components
corr_matrix_grouped = CreateCrossCorrMatrix(C_gr');
figure
imshow(corr_matrix_grouped,[0.9 1]);
%colormap(jet);
axis on;


%% Cluster spatial components

%# remove diagonal elements
corr_matrix_grouped(eye(size(corr_matrix_grouped))~=0)=0;
%# and convert to a vector (as pdist)
dissimilarity = 1 - corr_matrix_grouped(find(corr_matrix_grouped))';

%# decide on a cutoff
%# remember that 0.4 corresponds to corr of 0.6!
cutoff = 0.1; 

%# perform complete linkage clustering
Z = linkage(dissimilarity,'complete');

%# group the data into clusters
%# (cutoff is at a correlation of 0.5)
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');

figure;
dendrogram(Z,0,'colorthreshold',cutoff)
