function[] = AuditoryParcellationInter(subjectID, session, datadir, mmpdir, outdir)

% script for individual-level parcellation of auditory subregions
% outputs the parcellation itself and average connectivity across each parcel
% NOTE: requires gifti toolbox and modified version of surfstat containing additional functions surfGetNeighbors, surfGeoDist, and networkComponents
% REQUIRED INPUTS
% 'subjectID' - subject identifier (e.g. 100307 for HCP data)
% 'session' - name of RS session (for HCP data, 'REST1' or 'REST2')
% 'datadir' - path to individual subject data directory (should contain one subfolder per subject e.g. 100307)
    % should include:
        % individual subject directories named as subjectIDs containing:
            % correlation matrix file - rfMRI_REST_left_corr_avg.gii.data
            % left and right hemisphere surfaces - e.g. 'lh.inflated' and 'rh.inflated'
% 'mmpdir' - path to group-level data and ROIs
    % should include:
        % group connectivity maps
        % MMP results from Glasser et al. 2016: 'Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'
        % masks of left and right hemispheres: 'L.atlasroi.32k_fs_LR.shape.gii' and 'R.atlasroi.32k_fs_LR.shape.gii'
% 'outdir' - path to the desired output directory

subjectID = num2str(subjectID);
wb_command = '/home/estrid/Software/workbench/bin_linux64/wb_command'; % should point to directory where wb_command is installed

% % uncomment this for testing on a single subject
% datadir = ('/home1/estrid/Data/HCP/HCP_S1200_glyphsets/');
% mmpdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');
% outdir = ('/home/estrid/Data/HCP/AuditoryParcellation/');
% subjectID = '100307';
% session = ('REST1');

%% get surface for individual subject

disp('importing surfaces...')
surf = SurfStatReadSurf( {...
[datadir subjectID '/lh.inflated'],...
[datadir subjectID '/rh.inflated']} );

% get masks of left and right hemispheres excluding medial wall
lmask = gifti([mmpdir 'L.atlasroi.32k_fs_LR.shape.gii']);
lmask = lmask.cdata;
rmask = gifti([mmpdir 'R.atlasroi.32k_fs_LR.shape.gii']);
rmask = rmask.cdata;
hemimask = vertcat(lmask, rmask);
hemisize = zeros(size(rmask));
lmask = vertcat(lmask, hemisize);
rmask = vertcat(hemisize, rmask);

%% get connectivity matrix for individual subject

disp('loading connectivity matrix...')
tic
fid = fopen([datadir subjectID '/rfMRI_' session '_corr_avg.gii.data'], 'r');
M = fread(fid,[59412 59412], 'float32');
toc

%% get MMP labels (Glasser parcellation) and extract auditory ROIs
% group spatially contiguous areas that have almost identical connectivity based on group level analysis (results in 9 areas per hemi to parcellate)
% to check which parcel is which, run wb_command -file-information
% indices 181 and up is left hemisphere

% mmp = ciftiopen([mmpdir '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'], wb_command);
% mmp = mmp.cdata;
% filename = sprintf('%sMMP.', mmpdir);
% fid = fopen([filename, '1D'],'w');
% fprintf(fid, '%u\n', mmp);
% fclose(fid);

mmp = importdata([mmpdir 'MMP.1D']);

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

%% combine the 9 areas per hemi into a single vector and exclude the medial wall

map = zeros(length(rois),1);
for i = 1:length(areas)
    map(rois(:,i)== 1) = i;
end
labels = zeros(length(hemimask),1);
labels(hemimask==1) = map'; % the 18 areas/labels
ROI = map>0; % the temporal ROI/area to be parcellated

%%  dilate the labels by computing geodesic distance from each label on midthickness surface

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
    dist = surfGeoDist(surf2, label');
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
% results in 18 partial correlation maps representing the similarity of each vertex to each of the connectivity templates
disp('running partial correlations for group templates...')
tic
partcorr_areas = zeros(size(M,1), size(templates,2));
for i = 1:size(templates, 2)
    rm = templates;
    rm(:,i) = [];
    roi = logical(labels_dilated(:,i));
    roi = roi(logical(hemimask));
    M_area = M(:, roi);
    partcorr = partialcorr(M_area, templates(:,i), rm);
    results = zeros(length(M),1);
    results(roi) = partcorr;
    partcorr_areas(:,i) = results;
end
toc

%% extract individual-level template maps using maximum partial correlation
% for each of the 18 partial correlation maps, the maximum is identified and the connectivity of that vertex in the individual is extracted to be
% used as the individual-level connectivity template

templates_indiv = zeros(size(M,1), size(areas,2));
for i = 1:size(partcorr_areas, 2)
    [~,maxnode] = max(partcorr_areas(:,i));
    indivconn = M(:,maxnode);
    templates_indiv(:,i) = indivconn;
end

% % or do the same but exclude local connections
% templates_indiv = zeros(size(M,1), size(areas,2));
% for i = 1:size(partcorr_areas, 2)
%     [~,maxnode] = max(partcorr_areas(:,i));
%     indivconn = M(:,maxnode);
%     roi = logical(rois(:,i));
%     indivconn(roi) = 0;
%     templates_indiv(:,i) = indivconn;
% end

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

            % get network components (clusters) for label
            [~,~,members] = networkComponents(A);
            
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
            
            % calculate geodesic distance from the label with highest overlap
            [~,maxoverlap] = max(DC);
            maxclust = zeros(size(results));
            maxclust(members{maxoverlap}) = 1;
            dist = surfGeoDist_parcellation(surf2, maxclust');

            % invert, normalize, and threshold distances
            dist = dist(:,logical(labels));
            dist = dist.*-1; %invert
            dist = 2*(dist - min(dist)) ./ (max(dist) - min(dist))-1; %normalize -1 to 1
            dist(dist < 0) = 0; %remove negative values
            geodist_label = zeros(length(labels),1);
            geodist_label(logical(labels)) = dist;
            if h <= length(areas)/2
                geodist_label = geodist_label.*lmask;
            else
                geodist_label = geodist_label.*rmask;
            end
            % apply spatial weighting to partial correlations using normalized geodesic distances
            partcorr_weighted(:,h) = partcorr_weighted(:,h).*geodist_label;
        end
    end
end

%% rerun winner-take-all partition with spatially weighted partial correlation maps

[val,part] = max(partcorr_weighted,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

%% Match left and right areas

partviz = part;
for i=1:length(partviz)
    if partviz(i)>length(areas)/2
        partviz(i)= partviz(i)-length(areas)/2;
    end
end
        
%% Viz and save

figure('visible','off');
SurfStatViewData_LR(partviz, surf, [min(partviz) max(partviz)], 'WTA indiv');

filename = sprintf('%s%s_AuditoryParcellation_%s.', outdir, subjectID, session);
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
    
    % viz and save
    figure('visible','off');
    SurfStatViewData_LRml(results, surf, [0 1], ['Mean Connectivity' areas{i}]);

    filename = sprintf('%s%s_AuditoryParcellation_%s_meanconn_%s.', outdir, subjectID, session, areas{i});
    saveas(gcf, [filename 'png']); close all;
    fid = fopen([filename, '1D'],'w');
    fprintf(fid, '%u\n', results);
    fclose(fid);
end

