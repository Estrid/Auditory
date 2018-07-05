function[] = AuditoryParcellationInter(subjectID, session, datadir, mmpdir, outdir)
%
subjectID = num2str(subjectID);

mkdir /home/estrid/Data/HCP/MMP_indiv

% % uncomment this for testing on a single subject
% datadir = ('/home/estrid/Data/HCP/HCP_S1200_glyphsets/');
% mmpdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');
% outdir = ('/home/estrid/Data/HCP/MMP_indiv/');
% subjectID = '100307';
% session = ('REST1');

Hemi = ('LR');

%% get surface for individual subject

surf = SurfStatReadSurf( {...
[datadir subjectID '/lh.inflated'],...
[datadir subjectID '/rh.inflated']} );


%% get connectivity matrix for individual subject
disp('loading connectivity matrix...')
tic
fid = fopen([datadir subjectID '/rfMRI_' session '_corr_avg.gii.data'], 'r');
M = fread(fid,[59412 59412], 'float32');
toc

%% get MMP labels (Glasser parcellation) and extract auditory ROIs

% to check which parcel is which, run the following command in terminal:
%/home/estrid/Software/workbench/bin_linux64/wb_command -file-information /home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii% >181 is left hemi
%181 and up is left hemi

mmp = ciftiopen('/home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii', '/home/estrid/Software/workbench/bin_linux64/wb_command');
mmp = mmp.cdata;

lmask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/L.atlasroi.32k_fs_LR.shape.gii');
lmask = lmask.cdata;
rmask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/R.atlasroi.32k_fs_LR.shape.gii');
rmask = rmask.cdata;
hemimask = vertcat(lmask, rmask);
hemisize = zeros(size(rmask));
lmask = vertcat(lmask, hemisize);
rmask = vertcat(hemisize, rmask);

%% use 9 areas

PAC_L = sum(double(mmp == [204 353 354]),2);
PAC_R = sum(double(mmp == [24 173 174]),2);

SAC_L = sum(double(mmp == [304 287 355]),2);
SAC_R = sum(double(mmp == [124 107 175]),2);

areas_val_L = [284 305 303 308 309 356 310];
areas_val_R = [104 125 123 128 129 176 130];

rois_L = zeros(length(mmp), length(areas_val_L));
for i = 1:length(areas_val_L)
    n = areas_val_L(i);
    roi = mmp == n;
    rois_L(:,i) = roi;
end

rois_R = zeros(length(mmp), length(areas_val_R));
for i = 1:length(areas_val_R)
    n = areas_val_R(i);
    roi = mmp == n;
    rois_R(:,i) = roi;
end

areas_L = {'PAC_L', 'SAC_L', 'RI_L', 'A5_L', 'STGa_L', 'STSda_L', 'STSdp_L', 'STSva_L', 'STSvp_L'};
areas_R = {'PAC_R', 'SAC_R', 'RI_R', 'A5_R', 'STGa_R', 'STSda_R', 'STSdp_R', 'STSva_R', 'STSvp_R'};

areas = horzcat(areas_L, areas_R);
rois = horzcat(horzcat(PAC_L, SAC_L, rois_L), horzcat(PAC_R, SAC_R, rois_R));

%% or combine similar regions according to spatial correlation
% 
% PAC_L = sum(double(mmp == [204 284 304 353 354]),2);
% PAC_R = sum(double(mmp == [24 104 124 173 174]),2);
% 
% SAC1_L = sum(double(mmp == [355 287]),2);
% SAC1_R = sum(double(mmp == [175 107]),2);
% 
% SAC2_L = sum(double(mmp == [303 305]),2);
% SAC2_R = sum(double(mmp == [123 125]),2);
% 
% SAC3_L = sum(double(mmp == [308 309]),2);
% SAC3_R = sum(double(mmp == [128 129]),2);
% 
% SAC4_L = sum(double(mmp == [310]),2);
% SAC4_R = sum(double(mmp == [130]),2);
% 
% SAC5_L = sum(double(mmp == [356]),2);
% SAC5_R = sum(double(mmp == [176]),2);
% 
% areas_L = {'PAC_L', 'SAC1_L', 'SAC2_L', 'SAC3_L', 'SAC4_L', 'SAC5_L'};
% areas_R = {'PAC_R', 'SAC1_R', 'SAC2_R', 'SAC3_R', 'SAC4_R', 'SAC5_R'};
% rois_L = horzcat(PAC_L, SAC1_L, SAC2_L, SAC3_L, SAC4_L, SAC5_L);
% rois_R = horzcat(PAC_R, SAC1_R, SAC2_R, SAC3_R, SAC4_R, SAC5_R);
% 
% areas = horzcat(areas_L, areas_R);
% rois = horzcat(rois_L, rois_R);

%% include nonareas
% 
% nonareas_L = {'PFcm_L', 'Ig_L', '52_L', 'PSL_L', 'STV_L', 'TPOJ1_L', ...
%     'PHT_L', 'TE1p_L', 'TE1m_L', 'TE1a_L', 'TGd_L', 'PI_L'};
% nonareas_R = {'PFcm_R','Ig_R', '52_R', 'PSL_R', 'STV_R', 'TPOJ1_R', ...
%     'PHT_R', 'TE1p_R', 'TE1m_R', 'TE1a_R', 'TGd_R', 'PI_R'};
% nonareas_val_L = [285 348 283 205 208 319 317 313 357 312 311 358];
% nonareas_val_R = [105 168 103 25 28 139 137 133 177 132 131 178];
% 
% nonarearois_L = zeros(length(mmp), length(nonareas_val_L));
% for i = 1:length(nonareas_val_L)
%     n = nonareas_val_L(i);
%     roi = mmp == n;
%     nonarearois_L(:,i) = roi;
% end
% nonarearois_R = zeros(length(mmp), length(nonareas_val_R));
% for i = 1:length(nonareas_val_R)
%     n = nonareas_val_R(i);
%     roi = mmp == n;
%     nonarearois_R(:,i) = roi;
% end
% 
% 
% % non1_L = sum(double(mmp == [285 283 358 348]),2);
% % non1_R = sum(double(mmp == [105 103 178 168]),2);
% % 
% % non2_L = sum(double(mmp == [208 319 205]),2);
% % non2_R = sum(double(mmp == [28 139 25]),2);
% % 
% % non3_L = sum(double(mmp == [311 312]),2);
% % non3_R = sum(double(mmp == [131 132]),2);
% % 
% % non4_L = sum(double(mmp == [313]),2);
% % non4_R = sum(double(mmp == [133]),2);
% % 
% % non5_L = sum(double(mmp == [357]),2);
% % non5_R = sum(double(mmp == [177]),2);
% % 
% % non6_L = sum(double(mmp == [317]),2);
% % non6_R = sum(double(mmp == [137]),2);
% % 
% % nonareas_L = {'non1_L', 'non2_L', 'non3_L', 'non4_L', 'non5_L', 'non6_L'};
% % nonareas_R = {'non1_R', 'non2_R', 'non3_R', 'non4_R', 'non5_R', 'non6_R'};
% % nonarearois_L = horzcat(non1_L, non2_L, non3_L, non4_L, non5_L, non6_L);
% % nonarearois_R = horzcat(non1_R, non2_R, non3_R, non4_R, non5_R, non6_R);
% 
% areas = horzcat(areas_L, nonareas_L, areas_R, nonareas_R);
% rois = horzcat(rois_L, nonarearois_L, rois_R, nonarearois_R);

%% 
map = zeros(length(rois),1);
for i = 1:length(areas)
    map(rois(:,i)== 1) = i;
end
labels = zeros(length(hemimask),1);
labels(hemimask==1) = map';
ROI = map>0;

%%  dilate the MMP labels by computing geodesic distance from each label on midthickness surface

surf2 = SurfStatReadSurf( {...
[datadir subjectID '/lh.midthickness'],...
[datadir subjectID '/rh.midthickness']} );

surf2 = surfGetNeighbors(surf2);
surf2.nbr(surf2.nbr==0)=1; %necessary for correct indexing

labels_dilated = zeros(size(labels,1), max(labels));
for i = 1:max(labels)
    label = labels;
    label(label~=i) = 0;
    label(label==i)=1;
    % calculate geodesic distance from label
    dist = surfGeoDist_parcellation(surf2, label');
    label_dilated = zeros(length(label), 1);
    label_dilated(dist<=5) = 1; %% change this to adjust dilation distance
    if i <= length(areas)/2
        label_dilated = label_dilated.*lmask;
    else
        label_dilated = label_dilated.*rmask;
    end
    labels_dilated(:,i) = label_dilated;
end

%figure; SurfStatViewData_LRml(label_dilated, surf, [0 1], 'label_dilated');
%figure; SurfStatViewData_LRml(dist, surf, [min(dist) max(dist)], 'geodist');


%% mask connectivity matrix by ROI

M_masked = M(:,ROI);
% M_masked(isnan(M_masked))=0;

%% get template connectivity maps

% load in each area template connectivity map and save it into a matrix
templates = zeros(size(M,1), size(areas,2));
for i = 1:length(areas)
    area = areas(i);
    area = char(area);
    conn = importdata([mmpdir 'mmp_meanconn_' area '.1D']);
    conn = conn(logical(hemimask));
    templates(:,i) = conn;
end

%% run vertex-wise partial correlation between vertices within each dilated label and group template maps
disp('running partial correlations for group templates...')
tic
partcorr_areas = zeros(size(M,1), size(templates,2));
for i = 1:size(templates, 2)
    rm = templates;
    rm(:,i) = [];
    roi = logical(labels_dilated(:,i));
    roi = roi(logical(hemimask));
    M_area = M(:,roi);
    partcorr = partialcorr(M_area, templates(:,i), rm);
    results = zeros(length(M),1);
    results(roi) = partcorr;
    partcorr_areas(:,i) = results;
end
toc

%% extract individual-level template maps using maximum partial correlation

templates_indiv = zeros(size(M,1), size(areas,2));
for i = 1:size(partcorr_areas, 2)
    [~,maxnode] = max(partcorr_areas(:,i));
    indivconn = M(:,maxnode);
    templates_indiv(:,i) = indivconn;
end

%% run vertex-wise partial correlation with indiv template maps for entire ROI
disp('running partial correlations for individual templates...')
tic
partcorr_areas_indiv = zeros(size(M,1), size(templates,2));
for i = 1:size(templates_indiv, 2)
    rm = templates_indiv;
    rm(:,i) = [];
    partcorr = partialcorr(M_masked, templates_indiv(:,i), rm);
    results = zeros(length(M),1);
    results(ROI) = partcorr;
    partcorr_areas_indiv(:,i) = results;
end
partcorr_areas_indiv(isnan(partcorr_areas_indiv)) = 0;
toc

%% apply winner-take-all partition to individual-level template partial correlation maps

[val,part] = max(partcorr_areas_indiv,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

results = zeros(length(hemimask),1);
results(hemimask==1) = part;

% figure;%('visible','off');
% SurfStatViewData_LR(results, surf, [min(results) max(results)], 'WTA indiv');

% filename = sprintf('%s%s_mmp_partcorr_WTA_indiv_%s.', outdir, subjectID, Hemi);
% saveas(gcf, [filename 'png']); close all;
% fid = fopen([filename, '1D'],'w');
% fprintf(fid, '%u\n', results);
% fclose(fid);

%% identify the best-overlapping clusters for each label and calculate geodesic distance from them
% then apply a spatial weighting to the partial correlation maps using inverse geodesic distance

partcorr_weighted = zeros(size(results,1), max(results));

for h = 1:max(results)
    label = results;
    label(label~=h) = 0;
    partcorr_weighted(logical(hemimask),h) = partcorr_areas_indiv(:,h);
    if any(label) == false % if label contains no ones i.e. doesn't exist
        partcorr_weighted(:,h) = 0; % set to zero
    elseif any(label) == true % if label exists, continue
        % create edge list for label
        edgelist = [];
        for i = 1:length(label)
            nbrs = surf2.nbr(:,i);
            nbrslabel = label(nbrs);
            for j = 1:length(nbrs)
                if label(i)~=0 && label(i) == nbrslabel(j)
                    edge = [i nbrs(j)];
                    edgelist = vertcat(edgelist, edge);
                end
            end
        end
        
        if isempty(edgelist) == true % e.g. this will happen if label consists of only single independent points
            partcorr_weighted(:,h) = 0; % set to zero
        else
            % build adjacency matrix from edge list
            A = sparse([edgelist(:,1),edgelist(:,2)],[edgelist(:,2),edgelist(:,1)],1);

            % get network components for label
            [~,~,members] = networkComponents(A);
%             % find the largest component of label
%             largest = members{1};
%             maxclust = zeros(size(results));
%             maxclust(largest) = 1;

%             % calculate geodesic distance from label
%             dist = surfGeoDist_parcellation(surf2, maxclust');
            
            % calculate spatial overlap of dilated label and each cluster to find the best fit
            DC = zeros(10,1);
            for k = 1:10
                clust = zeros(size(results));
                clust(members{k}) = 1;
                area = labels_dilated(:,h);
                area1 = single(clust==1);
                area2 = single(area==1);
                areasum = area1 + area2;
                areaoverlap = sum(single(areasum==2));
                area1 = sum(area1);
                area2 = sum(area2);
                dc = 2*areaoverlap/(area1+area2);
                DC(k,1)=dc;
            end
            
            % calculate geodesic distance from label with max overlap
            [~,maxoverlap] = max(DC);
            maxclust = zeros(size(results));
            maxclust(members{maxoverlap}) = 1;
            dist = surfGeoDist_parcellation(surf2, maxclust');

            % invert, normalize, and threshold distances
            dist = dist(:,logical(labels));
            dist = dist.*-1; %invert
            dist = 2*(dist - min(dist)) ./ (max(dist) - min(dist))-1; %normalize -1 to 1
            dist(dist < 0) = 0; %remove negative values
        %     dist = (dist - min(dist)) ./ (max(dist) - min(dist));
            geodist_label = zeros(length(labels),1);
            geodist_label(logical(labels)) = dist;
            if h <= length(areas)/2
                geodist_label = geodist_label.*lmask;
            else
                geodist_label = geodist_label.*rmask;
            end
            % SurfStatViewPAC_LR(geodist_label, surf, [min(geodist_label) max(geodist_label)], 'geodist');
            % apply spatial weighting to partial correlations using normalized geodesic distances
            partcorr_weighted(:,h) = partcorr_weighted(:,h).*geodist_label;
        end
    end
end

% test = partcorr_weighted(:,9);
% test = geodist_label;
% SurfStatViewPAC_LR(test, surf, [min(test) max(test)], 'partcorr');

%% rerun winner-take-all partition with spatial weighting

[val,part] = max(partcorr_weighted,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

%% Match left and right areas

partviz = part;
for i=1:length(partviz)
    if partviz(i)>length(areas)/2
        partviz(i)= partviz(i)-length(areas)/2;
    end
end

%% remove non-areas

% % if using all areas
% for i= [14:25 39:50]
%     partviz(partviz==i)=0;
% end

% % if using combined areas
% if length(areas)>6
%     partviz(partviz>6)=0;
% end
        
%% Viz and save

figure('visible','off');
SurfStatViewData_LR(partviz, surf, [min(partviz) max(partviz)], 'WTA indiv');

filename = sprintf('%s%s_mmp_partcorr_%s_%s.', outdir, subjectID, Hemi, session);
saveas(gcf, [filename 'png']); close all;
fid = fopen([filename, '1D'],'w');
fprintf(fid, '%u\n', part);
fclose(fid);

%% Calculate mean connectivity for each parcellated area and save

for i = 1:max(part)
    roi = sum(double(part == i),2);
    
    % remove medial wall (59k)
    roi = roi(logical(hemimask));
    
    % mask connectivity matrix
    Z = M(logical(roi),:);
    r = tanh(Z); %z to r transform
    Avgconn = mean(r, 1);
    
    % replace medial wall (64k)
    results = zeros(length(hemimask),1);
    results(hemimask==1) = Avgconn';
    
    % save as png and 1D
    figure('visible','off');
    SurfStatViewData_LRml(results, surf, [0 1], ['Mean Connectivity' areas{i}]);

    filename = sprintf('%s%s_mmp_meanconn_%s_%s_%s.', outdir, subjectID, Hemi, session, areas{i});
    saveas(gcf, [filename 'png']); close all;
    fid = fopen([filename, '1D'],'w');
    fprintf(fid, '%u\n', results);
    fclose(fid);
end

