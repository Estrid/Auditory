% script to comput average connectivity across MMP subregions in group average connectome
% takes into account inter-hemispheric connectivity, as opposed to only intra

%% define variables

mkdir /home/estrid/Data/HCP/GroupAvg820_MMP
outdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');

%% load in data

% load in surfaces
surf = SurfStatReadSurf( {...
    '/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated',...
    '/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated'} );

% or load in left and right hemi masks
lmask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/L.atlasroi.32k_fs_LR.shape.gii');
lmask = lmask.cdata;
rmask = gifti('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/R.atlasroi.32k_fs_LR.shape.gii');
rmask = rmask.cdata;
mask = vertcat(lmask, rmask);
%SurfStatViewPAC_L(lmask, surf_l, [0 1], 'mask')

% load in connectivity
conn = ciftiopen(['/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/HCP_S900_820_rfMRI_MSMAll_groupPCA_d4500ROW_zcorr_CORTEX.dconn.nii'], '/home/estrid/Software/workbench/bin_linux64/wb_command');
conn = conn.cdata;

%% get MMP labels (Glasser parcellation) and extract auditory ROIs

% to check which parcel is which, run the following command in terminal:
%/home/estrid/Software/workbench/bin_linux64/wb_command -file-information /home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii% >181 is left hemi
%181 and up is left hemi

mmp = ciftiopen('/home/estrid/Data/HCP/Glasser_et_al_2016_HCP_MMP1.0/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii', '/home/estrid/Software/workbench/bin_linux64/wb_command');
mmp = mmp.cdata;

areas = [204 284 304 353 354 287 303 305 308 309 310 355 356 24 104 124 173 ...
    174 107 123 125 128 129 130 175 176];
key = {'A1_L', 'RI_L', 'PBelt_L', 'MBelt_L', 'LBelt_L', 'TA2_L', 'STGa_L',...
    'A5_L', 'STSda_L', 'STSdp_L', 'STSvp_L', 'A4_L', 'STSva_L', 'A1_R',...
    'RI_R', 'PBelt_R', 'MBelt_R', 'LBelt_R', 'TA2_R', 'STGa_R', 'A5_R',...
    'STSda_R', 'STSdp_R', 'STSvp_R', 'A4_R', 'STSva_R'};
areas_mask = ismember(mmp, areas);

%% calculate mean connectivity across each mmp area
 
for i = 1:length(areas)
    n = areas(i);
    roi = mmp == n;
    Z = conn(roi,:);
    r = tanh(Z); %z to r transform
    Avgconn = mean(r, 1);
    
    % make results double surface-size (64k)
    results = zeros(length(mask),1);
    results(mask==1) = Avgconn';
    
    figure;%('visible','off');
    SurfStatViewData_LRml(results, surf, [0 1], ['Mean Connectivity' key{i}]);

    filename = sprintf('%smmp_meanconn_%s.', outdir, key{i});
    saveas(gcf, [filename 'png']); close all;
    fid = fopen([filename, '1D'],'w');
    fprintf(fid, '%u\n', results);
    fclose(fid);
end



