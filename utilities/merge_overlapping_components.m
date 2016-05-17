function [A,C,nr,merged_ROIs,P,S] = merge_overlapping_components(Y,A,b,C,f,P,S,Cn,options, display_numbers, max_number)
% Merge components that overlap spatially in a semi-manual fashion.
% This is desirable when multiple highly overlapping spatial components are found
% by the CNMF algorithm.
%
% Code by Eftychios A. Pnevmatikakis from plot_contours()and merge_comonents() reused heavily.
%
% Adapted by:
% Ingie Hong, Johns Hopkins Medical Institute, 2016

if ~isfield(options,'spatial_merge_thr') || isempty(options.merge_thr); spatial_merge_thr = 0.7; else spatial_merge_thr = options.spatial_merge_thr; end     % spatial overlap-based merging threshold
if ~isfield(options,'spatial_merge_bManual') || isempty(options.spatial_merge_bManual); spatial_merge_bManual = false; else spatial_merge_bManual = options.spatial_merge_bManual; end     % manual merging or not

if nargin < 11 || isempty(max_number)
    max_number = size(A,2);
else
    max_number = min(max_number,size(A,2));
end
if nargin < 10 || isempty(display_numbers)
    display_numbers = 1;
end

units = 'centimeters';
fontname = 'helvetica';
thr = 0.995;
d1=options.d1;
d2=options.d2;
cmap = hot(3*size(A,2));
CC = cell(size(A,2),1);
CR = cell(size(A,2),2);

% Get coordinates and plot spatial components
figure; imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
title(['Before merging, ' num2str(size(A,2)) ' components'])
axis tight; axis equal; 
hold on;

for i = 1:size(A,2)
    A_temp = full(reshape(A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend'); 
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        [CC{i},h(i)] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(i+size(A,2),:));
        %CC{i} = contourc(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)));
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        %CR{i,2} = A_temp(fp)';
    end
    hold on;
end

% Plot labels
cm = com(A(:,1:end),d1,d2);
if display_numbers
    lbl = strtrim(cellstr(num2str((1:size(A,2))')));
    %text(round(cm(1:max_number,2)),round(cm(1:max_number,1)),lbl(1:max_number),'color',[0,0,0],'fontsize',16,'fontname',fontname,'fontweight','bold');
    text(round(cm(1:max_number,2)),round(cm(1:max_number,1)),lbl(1:max_number),'color',[0,0,0],'fontsize',8,'fontname',fontname);
end

% Generate jsf
for i = 1:size(A,2);
    if ~isempty(CR{i,1})
        jsf(i) = struct('id',i,...
                    'coordinates',CR{i,1}'); %,...
%                     'values',CR{i,2},...
%                     'bbox',[min(CR{i,1}(1,:)),max(CR{i,1}(1,:)),min(CR{i,1}(2,:)),max(CR{i,1}(2,:))],...
%                     'centroid',cm(i,:));
    end
    if i == 1
        jsf = repmat(jsf,size(A,2),1);
    end
end


%% Find spatial overlap of spatial components and merge based on threshold

% Spatial overlap calculation (normalized to smaller spatial component)
%spatial_overlap_pixels=zeros(size(jsf,1),size(jsf,1)); 
spatial_overlap=zeros(size(jsf,1),size(jsf,1));

for i=1:size(jsf,1)
    for j=i+1:size(jsf,1)
        common_pixels=intersect(jsf(i).coordinates,jsf(j).coordinates,'rows');

        %spatial_overlap_pixels(i,j)= size(common_pixels,1); % For absolute pixel number quantification
        spatial_overlap(i,j)= size(common_pixels,1) / min([size(jsf(i).coordinates,1) size(jsf(j).coordinates,1)]); % normalized to smaller spatial component
    end
end

%% Cluster spatial components

dissimilarity = 1 - (spatial_overlap+spatial_overlap');
dissimilarity(eye(size(dissimilarity))~=0)=0;

% Spatial overlap cutoff
% remember that 0.3 corresponds to overlap of 0.7!
cutoff = 1-spatial_merge_thr; 

% perform complete linkage clustering
Z = linkage(dissimilarity,'single', 'correlation');

% group the data into clusters
% (cutoff is at a overlap of 0.7 as default)
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');

% plot clustering based on spatial overlap
figure;dendrogram(Z,0,'colorthreshold',cutoff)
title('Clustered groups based on spatial overlap')
ylabel('1 - spatial overlap')
xlabel('Spatial component')

%% Plot spatial components in same cluster

figure;
scrsz = get(groot,'ScreenSize');
set(gcf, 'Position',[50 50 scrsz(3)-100 scrsz(4)*2/3-150])
bUsePanel = exist('panel', 'file') ==2;

if bUsePanel 
    pan = panel();
    pan.pack('h', {1/3 []} )
else
    subplot(1,2,1);
end

spatial_overlap=spatial_overlap.*(spatial_overlap>spatial_merge_thr); % Ignore all overlap below threshold

for i=1:length(unique(groups))
    groupi = find(groups==i);
    if size(groupi,1)>1
        
        if bUsePanel
            pan(1).select()
        else
            subplot(1,2,1);
        end
        
        hold off;
        imagesc(Cn);
        hold on;
        cmap = autumn(size(groupi,1)+1);
        
        for j = 1:size(groupi,1)
            cont = medfilt1(CC{groupi(j)}')';
            if size(cont,2) > 1
                plot(cont(1,2:end),cont(2,2:end),'Color',cmap(j,:));  
            end
            
            if display_numbers
                lbl = strtrim(cellstr(num2str((1:size(A,2))')));
                text(round(cm(groupi(j),2)),round(cm(groupi(j),1)),lbl(groupi(j)),'Color',cmap(j,:),'fontsize',10,'fontname',fontname);
                
            end

        end
        
        title(['Spatial components: ' num2str(groupi')]);
        axis image;

        % Plot spatial components using A (unnormalized)       
        %         for j=1:size(groupi,1)
        %             unit=unit + reshape(A_or(:,groupi(j)), d1, d2);
        %             imagesc(unit)
        %             axis image; 
        %         end
        
        if bUsePanel
            pan(2).select()
        else
            subplot(1,2,2);
        end
        
        ax = gca;
        set(ax,'ColorOrder',cmap);
        ax.ColorOrderIndex = 1;
        plot(C(groupi,:)');
        legend([num2str(groupi)]);
        %pause
        drawnow
        
        if spatial_merge_bManual==true
        
            % Use a questdlg for user input on merging options
            choice = questdlg('Merge the following spatial componenets? (quit panel to abort merging)', ...
                'Merge_overlapping_components()', ...
                'Merge', 'Merge subset','No', 'No');
            % Handle response
            switch choice
                case 'Merge'
                    disp(['Merging: ' num2str(groupi') ])

                case 'Merge subset'
                    %disp([choice ' of:' num2str(groupi')])

                    [chosen_index,status] = listdlg('PromptString','Select spatial components to merge:',...
                    'InitialValue',[1:length(groupi)],...
                    'ListString',cellstr(num2str(groupi(:))));

                    if status==0||length(chosen_index)<2
                        disp('Not appropriate input, skipping..')
                        spatial_overlap(groupi,:)=0;   % Remove entries from spatial_overlap to avoid any merging
                    else
                        disp(['Merging subset: ' num2str(groupi(chosen_index)')]) 
                        spatial_overlap(setdiff([1:length(groupi)],chosen_index),:)=0;  % To merge only selected spatial componenets
                        spatial_overlap(:,setdiff([1:length(groupi)],chosen_index))=0;  % To merge only selected spatial componenets
                    end

                case 'No'
                    disp('Not merged.')
                    spatial_overlap(groupi,:)=0;  % Remove entries from spatial_overlap to avoid any merging
                case ''
                    error('Aborting')
            end
        end
        
    end
end


% Merge components and finish!
%size(A) % Size before merging
disp(['Removing ' num2str(nnz( max((spatial_overlap+spatial_overlap')>0,[],2))) ' components']); % how many units should be merged
[A,C,nr,merged_ROIs,P,S] = merge_components(Y,A,b,C,f,P,S,options, spatial_overlap);
%size(A)  % Size after merging

%% Replot spatial components after merge
% Get coordinates and plot spatial components
figure; imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
title(['After merging, ' num2str(size(A,2)) ' components'])
axis tight; axis equal; 
hold on;
CC = cell(size(A,2),1);
CR = cell(size(A,2),2);
cmap = hot(3*size(A,2));

for i = 1:size(A,2)
    A_temp = full(reshape(A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend'); 
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        [CC{i},h(i)] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(i+size(A,2),:));
        %CC{i} = contourc(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)));
%         fp = find(A_temp >= A_temp(ind(ff)));
%         [ii,jj] = ind2sub([d1,d2],fp);
%         CR{i,1} = [ii,jj]';
%         CR{i,2} = A_temp(fp)';
    end
    hold on;
end

if nargin < 11 || isempty(max_number)
    max_number = size(A,2);
else
    max_number = min(max_number,size(A,2));
end

% Plot labels
cm = com(A(:,1:end),d1,d2);
if display_numbers
    lbl = strtrim(cellstr(num2str((1:size(A,2))')));
    text(round(cm(1:max_number,2)),round(cm(1:max_number,1)),lbl(1:max_number),'color',[0,0,0],'fontsize',8,'fontname',fontname);
end

figure(gcf)



%% Distance measure based on center position 
% 
% center_dist=zeros(size(center,1),size(center,1));
% 
% for i=1:size(center,1)
%     for j=i+1:size(center,1)       
%         center_dist(i,j)= sum((center(i,:)-center(j,:)).^2)^0.5; 
%     end
% end



function [A,C,nr,merged_ROIs,P,S] = merge_components(Y,A,b,C,f,P,S,options, spatial_overlap)
% merging of spatially overlapping components (Adapted from
% merge_components.m)
%
% Inputs:
% Y:            raw data
% A:            matrix of spatial components
% b:            spatial background
% C:            matrix of temporal components
% f:            temporal background
% P:            struct for neuron parameters
% S:            deconvolved activity/spikes (optional)
% options:      struct for algorithm parameters

% Outputs:
% A:            matrix of new spatial components
% C:            matrix of new temporal components
% nr:           new number of components
% merged_ROIs:  list of old components that were merged
% P:            new parameter struct
% S:            matrix of new deconvolved/activity spikes

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% Adapted by:
% Ingie Hong, Johns Hopkins Medical Institute, 2016

defoptions = CNMFSetParms;
if nargin < 8; options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end       % # of columns
%if ~isfield(options,'merge_thr') || isempty(options.merge_thr); thr = defoptions.merge_thr; else thr = options.merge_thr; end     % merging threshold
%if ~isfield(options,'max_merg'); mx = 50; else mx = options.max_merg; end           % maximum merging operations
if ~isfield(options,'deconv_method') || isempty(options.deconv_method); options.deconv_method = defoptions.deconv_method; end
if ~isfield(options,'fast_merge') || isempty(options.fast_merge); options.fast_merge = defoptions.fast_merge; end  % flag for using fast merging

nr = size(A,2);
%[d,T] = size(Y);
d = size(A,1);
T = size(C,2);
%C_corr = corr(full(C(1:nr,:)'));                    % calculate correlation between temporal components
%FF1 = triu(C_corr)>= thr;                           % find graph of strongly correlated temporal components

%A_corr = triu(A(:,1:nr)'*A(:,1:nr));                % get spatial overlap between spatial components
%A_corr(1:nr+1:nr^2) = 0;                            % set diagonal to zero
%FF2 = A_corr > 0;                                   % find graph of overlapping spatial components

%FF3 = and(FF1,FF2);                                 % intersect the two graphs
FF3=spatial_overlap>0;
[l,c] = graph_connected_comp(sparse(FF3+FF3'));     % extract connected components
MC = [];
for i = 1:c
    if length(find(l==i))>1
        MC = [MC,(l==i)'];
    end
end

cor = zeros(size(MC,2),1);
for i = 1:length(cor)
    fm = find(MC(:,i));
    for j1 = 1:length(fm)
        for j2 = j1+1:length(fm)
            cor(i) = cor(i) + spatial_overlap(fm(j1),fm(j2));
        end
    end
end

[~,ind] = sort(cor,'descend');
%nm = min(length(ind),mx);   % number of merging operations
nm = length(ind);
merged_ROIs = cell(nm,1);
A_merged = zeros(d,nm);
C_merged = zeros(nm,T);
S_merged = zeros(nm,T);
if strcmpi(options.deconv_method,'constrained_foopsi')
    P_merged.gn = cell(nm,1);
    P_merged.b = cell(nm,1);
    P_merged.c1 = cell(nm,1);
    P_merged.neuron_sn = cell(nm,1);
end
if ~options.fast_merge
    Y_res = Y - A*C;
end

for i = 1:nm
    merged_ROIs{i} = find(MC(:,ind(i)));
    nC = sqrt(sum(C(merged_ROIs{i},:).^2,2));
    if options.fast_merge
        aa = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);
        for iter = 1:10
            cc = (aa'*A(:,merged_ROIs{i}))*C(merged_ROIs{i},:)/sum(aa.^2);
            aa = A(:,merged_ROIs{i})*(C(merged_ROIs{i},:)*cc')/norm(cc)^2;
        end
        na = sqrt(sum(aa.^2));
        aa = aa/na;
        %[cc,b_temp,c1_temp,g_temp,sn_temp,ss] = constrained_foopsi(cc);
        cc = na*cc';
        ss = cc;
    else
        A_merged(:,i) = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);    
        Y_res = Y_res + A(:,merged_ROIs{i})*C(merged_ROIs{i},:);
        ff = find(A_merged(:,i));
        Pmr = P;
        if isfield(Pmr,'unsaturatedPix');
            px = intersect(Pmr.unsaturatedPix,ff);
            Pmr.unsaturatedPix = zeros(length(px),1);
            for pxi = 1:length(px)
                Pmr.unsaturatedPix(pxi) = find(ff == px(pxi));
            end
        end
        cc = update_temporal_components(Y_res(ff,:),A_merged(ff,i),b(ff,:),median(spdiags(nC,0,length(nC),length(nC))\C(merged_ROIs{i},:)),f,Pmr,options);
        [aa,bb] = update_spatial_components(Y_res,cc,f,A_merged(:,i),P,options);    
        [cc,~,Ptemp,ss] = update_temporal_components(Y_res(ff,:),aa(ff),bb(ff,:),cc,f,Pmr,options);
    end
    A_merged(:,i) = aa;    
    C_merged(i,:) = cc;
    S_merged(i,:) = ss;
    if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
        if options.fast_merge
            P_merged.gn{i} = 0; %g_temp;   % do not perform deconvolution during merging
            P_merged.b{i} = 0;  %b_temp;
            P_merged.c1{i} = 0; %c1_temp;
            P_merged.neuron_sn{i} = 0; %sn_temp;
        else
            P_merged.gn{i} = Ptemp.gn{1};
            P_merged.b{i} = Ptemp.b{1};
            P_merged.c1{i} = Ptemp.c1{1};
            P_merged.neuron_sn{i} = Ptemp.neuron_sn{1};
            if i < nm
                Y_res(ff,:) = Y_res(ff,:) - aa(ff)*cc;
            end
        end
    end
end

neur_id = unique(cell2mat(merged_ROIs));

A = [A(:,1:nr),A_merged,A(:,nr+1:end)];
C = [C(1:nr,:);C_merged;C(nr+1:end,:)];
A(:,neur_id) = [];
C(neur_id,:) = [];

if nargin < 7
    S = [];
    if nargout == 6
        warning('Merged spikes matrix is returned as empty because the original matrix was not provided.');
    end
else
    S = [S(1:nr,:);S_merged];
    S(neur_id,:) = [];
end

if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
    P.b(neur_id) = [];
    P.b(nr - length(neur_id) + (1:nm)) = P_merged.b;
    P.gn(neur_id) = [];
    P.gn(nr - length(neur_id) + (1:nm)) = P_merged.gn;
    P.c1(neur_id) = [];
    P.c1(nr - length(neur_id) + (1:nm)) = P_merged.c1;
    P.neuron_sn(neur_id) = [];
    P.neuron_sn(nr - length(neur_id) + (1:nm)) = P_merged.neuron_sn;
end
nr = nr - length(neur_id) + nm;