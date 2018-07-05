%% do a winner-take-all across all area probability maps to get a group parcellation

% import group prob maps for each area and write into a matrix

datadir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
outdir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
session = ('REST2');

surf = SurfStatReadSurf( {...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated'),...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated')} );

areas = {'PAC', 'SAC', 'RI', 'A5', 'STGa', 'STSda', 'STSdp', 'STSva', 'STSvp'};

probmat = zeros(size(surf.coord,2), length(areas));

for i = 1:length(areas)
    area = areas(i);
    area = char(area);
    prob = importdata([datadir 'GroupProb_' area '_LR_' session '.1D']);
    probmat(:,i)=prob;
end

[val,part] = max(probmat,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

figure('visible','off');
SurfStatViewData_LR(part, surf, [min(part) max(part)], 'Group parcellation');

filename = sprintf('%sGroupProbWTA_LR_%s', outdir, session);
saveas(gcf, [filename '.png']); close all;
fid = fopen([filename '.1D'],'w');
fprintf(fid, '%u\n', part);
fclose(fid);