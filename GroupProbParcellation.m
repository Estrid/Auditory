% function[] = GroupProbParcellation(subject_list, hemi, datadir)

datadir = ('/home/estrid/Data/HCP/MMP_indiv/');
outdir = ('/home/estrid/Data/HCP/MMP_indiv/');
subject_list = importdata('/home/estrid/scripts/subjectlist.txt');
hemi = 2;

if hemi == 1
    Hemi = ('L');
    surf = SurfStatReadSurf1('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated');
elseif hemi == 2
    Hemi = ('R');
    surf = SurfStatReadSurf1('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated');
end

areas = {'A1', 'RI', 'PBelt', 'MBelt', 'LBelt', ...
'TA2', 'STGa', 'A5', 'STSda', 'STSdp', 'STSvp', 'A4', 'STSva'};
%% create group probability maps for each area

part_prob = zeros(size(surf.coord,2), length(areas));
errors = 0;
for h = 1:length(subject_list)
    subjectID = num2str(subject_list(h));
    try
    part = importdata([datadir subjectID '_mmp_partcorr_WTA_indiv_' Hemi '.1D']);
    for i = 1:max(part)
        part_area = double(part==i);
        part_prob(:,i) = part_prob(:,i) + part_area;
    end
    catch
        fprintf('error in subject # %u\n', subject_list(h))
        errors = errors+1;
        continue
    end
    fprintf('subject # %u\n',subject_list(h))
end


% create figure and save 1D
for i = 1:size(part_prob,2)
    area = char(areas(i));
    num = max(part_prob(:,i));
    part_prob(:,i) = part_prob(:,i)./num; % calculate probability from overlap of labeled subjects
    figure('visible','off');
    if hemi == 1
        SurfStatViewPAC_L(part_prob(:,i), surf, [0 1], [area ' ' Hemi ' --' num2str(num) ' subjects']);
    elseif hemi == 2
        SurfStatViewPAC_R(part_prob(:,i), surf, [0 1], [area ' ' Hemi ' --' num2str(num) ' subjects']);
    end

    filename = sprintf('%sGroupProb_%s_%s.', outdir, area, Hemi);
    saveas(gcf, [filename 'png']); close all;
    fid = fopen([filename, '1D'],'w');
    fprintf(fid, '%u\n', part_prob(:,i));
    fclose(fid);
end





    
