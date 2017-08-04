function[] = AuditoryParcellation(subjectID, hemi, datadir, mmpdir, outdir)
subjectID = num2str(subjectID);

% mkdir /home/estrid/Data/HCP/MMP_indiv
% outdir = ('/home/estrid/Data/HCP/MMP_indiv/');
% datadir = ('/home/estrid/pamina_mount/HCP/HCP_S1200_glyphsets/');
% mmpdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');
% subjectID = '124422';
% hemi = 2;

if hemi == 1
    Hemi = ('L');
elseif hemi == 2
    Hemi = ('R');
end

%% get surface for individual subject

if hemi == 1
    surf = SurfStatReadSurf1([datadir subjectID '/lh.inflated']);
elseif hemi == 2
    surf = SurfStatReadSurf1([datadir subjectID '/rh.inflated']);
end

%% get connectivity matrix for individual subject
disp('loading connectivity matrix...')
tic
if hemi == 1
    fid = fopen([datadir subjectID '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
    M = fread(fid,[32492 32492], 'float32');
elseif hemi == 2
    fid = fopen([datadir subjectID '/rfMRI_REST_right_corr_avg.gii.data'], 'r');
    M = fread(fid,[32492 32492], 'float32');
end
toc

%% get MMP labels (Glasser parcellation) and extract auditory ROIs

% to check which parcel is which, run the following command in terminal:
%/home/estrid/Software/workbench/bin_linux64/wb_command -file-information /home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii% >181 is left hemi
%181 and up is left hemi

mmp = ciftiopen('/home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii', '/home/estrid/Software/workbench/bin_linux64/wb_command');
mmp = mmp.cdata;

if hemi == 1
    areas = {'A1_L', 'RI_L', 'PBelt_L', 'MBelt_L', 'LBelt_L', ...
        'TA2_L', 'STGa_L', 'A5_L', 'STSda_L', 'STSdp_L', 'STSvp_L', 'A4_L', 'STSva_L'};
    areas_val = [204 284 304 353 354 287 303 305 308 309 310 355 356];
    mmp = mmp(mmp>=181);
    hemimask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/L.atlasroi.32k_fs_LR.shape.gii');
    hemimask = hemimask.cdata;
elseif hemi == 2        
    areas = {'A1_R', 'RI_R', 'PBelt_R', 'MBelt_R', 'LBelt_R', ...
        'TA2_R', 'STGa_R', 'A5_R', 'STSda_R', 'STSdp_R', 'STSvp_R', 'A4_R', 'STSva_R'};
    areas_val = [24 104 124 173 174 107 123 125 128 129 130 175 176];
    mmp = mmp(mmp<=180);
    hemimask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/R.atlasroi.32k_fs_LR.shape.gii');
    hemimask = hemimask.cdata;
end

rois = double(mmp == areas_val);
map = zeros(length(rois),1);
for i = 1:length(areas_val)
    map(rois(:,i)== 1) = i;
end
labels = zeros(length(hemimask),1);
labels(hemimask==1) = map';

ROI = labels>0;

%%  dilate the MMP labels by computing geodesic distance from each label on midthickness surface

if hemi == 1
    surf2 = SurfStatReadSurf1([datadir subjectID '/lh.midthickness']);
elseif hemi == 2
    surf2 = SurfStatReadSurf1([datadir subjectID '/rh.midthickness']);
end

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
    labels_dilated(:,i) = label_dilated;
end

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

%% run vertex-wise partial correlation with indiv template maps
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

% SurfStatViewPAC_L(partcorr_areas_indiv(:,1), surf, [min(partcorr_areas_indiv(:,1)) max(partcorr_areas_indiv(:,1))], 'Partial Correlation');

%% apply winner-take-all partition to individual-level template partial correlation maps

[val,part] = max(partcorr_areas_indiv,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% figure;%('visible','off');
% if hemi == 1
%     SurfStatViewPAC_L(part, surf, [min(part) max(part)], 'WTA indiv L');
% elseif hemi == 2
%     SurfStatViewPAC_R(part, surf, [min(part) max(part)], 'WTA indiv R');
% end

% filename = sprintf('%s%s_mmp_partcorr_WTA_indiv_%s.', outdir, subjectID, Hemi);
% saveas(gcf, [filename 'png']); close all;
% fid = fopen([filename, '1D'],'w');
% fprintf(fid, '%u\n', part);
% fclose(fid);

%% identify the largest clusters for each label and calculate geodesic distance from them
% then apply a spatial weighting to the partial correlation maps using inverse geodesic distance

geodist_part = zeros(size(part,1), max(part));
partcorr_weighted = partcorr_areas_indiv;

for h = 1:max(part)
    label = part;
    label(label~=h) = 0; 
    if any(label) == false % if label contains no ones
        partcorr_weighted(:,h) = partcorr_weighted(:,h); % do nothing
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

        % build adjacency matrix from edge list
        A = sparse([edgelist(:,1),edgelist(:,2)],[edgelist(:,2),edgelist(:,1)],1);

        % get network components for label
        [~,~,members] = networkComponents(A);
        % find the largest component of label
        largest = members{1};
        maxclust = zeros(size(part));
        maxclust(largest) = 1;

        % calculate geodesic distance from label
        dist = surfGeoDist_parcellation(surf2, maxclust');

        % invert, normalize, and threshold distances
        dist = dist(:,ROI);
        dist = dist.*-1;
        dist = 2*(dist - min(dist)) ./ (max(dist) - min(dist))-1;
        dist(dist < 0) = 0;
    %     dist = (dist - min(dist)) ./ (max(dist) - min(dist));
        geodist_label = zeros(length(M),1);
        geodist_label(ROI) = dist;
    %     SurfStatViewPAC_L(geodist_label, surf, [min(geodist_label) max(geodist_label)], 'geodist');
        geodist_part(:,h) = geodist_label;
        % apply spatial weighting to partial correlations using normalized geodesic distances
        partcorr_weighted(:,h) = partcorr_weighted(:,h).*geodist_label;
    end
end

%% rerun winner-take-all partition with spatial weighting

[val,part] = max(partcorr_weighted,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

figure('visible','off');
if hemi == 1
    SurfStatViewPAC_L(part, surf, [min(part) max(part)], 'WTA indiv L');
elseif hemi == 2
    SurfStatViewPAC_R(part, surf, [min(part) max(part)], 'WTA indiv R');
end

filename = sprintf('%s%s_mmp_partcorr_WTA_indiv_%s.', outdir, subjectID, Hemi);
saveas(gcf, [filename 'png']); close all;
fid = fopen([filename, '1D'],'w');
fprintf(fid, '%u\n', part);
fclose(fid);

%% get sulc map, normalize, and save as 1D for viz

if hemi == 1
    sulc = gifti(['/home/estrid/pamina_mount/HCP/HCP_S1200_subjects/' subjectID '/MNINonLinear/fsaverage_LR32k/' subjectID '.L.sulc.32k_fs_LR.shape.gii']);
elseif hemi == 2
    sulc = gifti(['/home/estrid/pamina_mount/HCP/HCP_S1200_subjects/' subjectID '/MNINonLinear/fsaverage_LR32k/' subjectID '.R.sulc.32k_fs_LR.shape.gii']);
end
sulc = sulc.cdata;

datamin = min(sulc);
datamax = max(sulc);
datarange = datamax - datamin;
sulcnorm = (sulc - datamin) / datarange;

filename = sprintf('%s/%s/sulcnorm_%s.', datadir, subjectID, Hemi);
fid = fopen([filename, '1D'],'w');
fprintf(fid, '%u\n', sulcnorm);
fclose(fid);