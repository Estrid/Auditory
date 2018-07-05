function[] = GroupProbParcellationInter(subject_list, datadir, outdir, session, outext)

% subject_list = importdata(subject_list);
% datadir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
% outdir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
% session = 'REST1';

ext = '_mmp_partcorr_LR_'; %extension of files to import

surf = SurfStatReadSurf( {...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated'),...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated')} );

areas_L = {'PAC_L', 'SAC_L', 'RI_L', 'A5_L', 'STGa_L', 'STSda_L', 'STSdp_L', 'STSva_L', 'STSvp_L'};
areas_R = {'PAC_R', 'SAC_R', 'RI_R', 'A5_R', 'STGa_R', 'STSda_R', 'STSdp_R', 'STSva_R', 'STSvp_R'};

areas = horzcat(areas_L, areas_R);

%% create group probability maps for each area

part_prob = zeros(size(surf.coord,2), length(areas));
errors = 0;
fid = fopen('errors.txt','w');
for h = 1:length(subject_list)
    subjectID = num2str(subject_list(h));
    try
    part = importdata([datadir subjectID ext session '.1D']);
    catch
        fprintf(fid, 'could not open subject # %u\n', subject_list(h));
        errors = errors+1;
        continue
    end
    for i = 1:max(part)
        part_area = double(part==i);
        if max(part_area)<1
            fprintf(fid, 'area %s does not exist in subject # %u\n', char(areas(i)), subject_list(h));
            errors = errors+1;
        else
        end
        part_prob(:,i) = part_prob(:,i) + part_area;
    end
%     fprintf('subject # %u\n',subject_list(h))
end
fclose(fid);

% combine left and right prob maps for visualization
part_prob_LR = part_prob(:,1:size(part_prob,2)/2) + part_prob(:,(size(part_prob,2)/2+1):size(part_prob,2));
areas = {'PAC', 'SAC', 'RI', 'A5', 'STGa', 'STSda', 'STSdp', 'STSva', 'STSvp'};

% create figure and save 1D
for i = 1:size(part_prob_LR,2)
    area = char(areas(i));
    maxL = max(part_prob_LR(1:length(part_prob)/2,i));
    maxR = max(part_prob_LR(length(part_prob)/2+1:length(part_prob),i));
    part_prob_LR(:,i) = part_prob_LR(:,i)./length(subject_list); % calculate probability from overlap of labeled subjects
    figure('visible','off');
    SurfStatViewData_LR(part_prob_LR(:,i), surf, [0 1], [area ' LR --maxL' num2str(maxL) ' maxR' num2str(maxR)]);

    filename = sprintf('%sGroupProb_%s_LR_%s', outdir, area, session);
    if exist('outext', 'var')
        saveas(gcf, [filename outext '.png']); close all;
        fid = fopen([filename outext '.1D'],'w');
    else
        saveas(gcf, [filename '.png']); close all;
        fid = fopen([filename '.1D'],'w');
    end
    fprintf(fid, '%u\n', part_prob_LR(:,i));
    fclose(fid);
end

% %threshold at 0.75 and visualize to check how spatially distinct areas are
% part_prob50 = part_prob_LR;
% part_prob50(part_prob50<50)=0;
% part_prob50 = sum(part_prob50, 2);
% figure('visible','off');
% SurfStatViewData_LR(part_prob50, surf, [min(part_prob50) max(part_prob50)], '50% thresholded probability');
% filename = sprintf('%sGroupProb_allareas50_LR_%s', outdir, session);
% if exist('outext', 'var')
%     saveas(gcf, [filename outext '.png']); close all;
%     fid = fopen([filename outext '.1D'],'w');
% else
%     saveas(gcf, [filename '.png']); close all;
%     fid = fopen([filename '.1D'],'w');
% end
% fprintf(fid, '%u\n', part_prob50);
% fclose(fid);


    
