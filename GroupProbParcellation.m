function[] = GroupProbParcellation(subject_list, session, datadir, mmpdir, outdir)

surf = SurfStatReadSurf( {...
[mmpdir 'lh.very_inflated'],...
[mmpdir 'rh.very_inflated']} );

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
    part = importdata([datadir subjectID '_AuditoryParcellation_' session '.1D']);
    catch
        fprintf(fid, 'could not open subject # %u\n', subject_list(h));
        errors = errors+1;
        continue
    end
    for i = 1:max(part)
        part_area = double(part==i);
        if max(part_area)<1
            fprintf(fid, 'error in subject # %u - %s does not exist\n', subject_list(h), char(areas(i)));
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

    filename = sprintf('%sGroupProb_%s_%s', outdir, area, session);
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
