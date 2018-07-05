function[fitscore, tmax] = GroupConnParcellationInter(subject_list, datadir, outdir, session, outext)

% subject_list = importdata('/home/estrid/scripts/bash/subjectlist.txt');
% datadir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
% outdir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
% session = ('REST1');

% subject_list = importdata(subject_list);
ext = '_mmp_meanconn_LR_'; %extension of files to import

surf = SurfStatReadSurf( {...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated'),...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated')} );

areas_L = {'PAC_L', 'SAC_L', 'RI_L', 'A5_L', 'STGa_L', 'STSda_L', 'STSdp_L', 'STSva_L', 'STSvp_L'};
areas_R = {'PAC_R', 'SAC_R', 'RI_R', 'A5_R', 'STGa_R', 'STSda_R', 'STSdp_R', 'STSva_R', 'STSvp_R'};

areas = horzcat(areas_L, areas_R);
fitscore = zeros(length(subject_list), length(areas));
tmax = zeros(size(areas));

for i = 1:length(areas)
    area = areas(i);
    area = char(area);
    groupconn = zeros(size(surf.coord,2), length(subject_list));
    errors = 0;
    for h = 1:length(subject_list)
        subjectID = num2str(subject_list(h));
        try
        conn = importdata([datadir subjectID ext session '_' area '.1D']);
        groupconn(:,h) = conn;
        catch
            fprintf('error in subject # %u\n', subject_list(h))
            errors = errors+1;
            continue
        end
        fprintf('subject # %u\n',subject_list(h))
    end

    groupconn = groupconn';
    slm = SurfStatLinMod(groupconn, 1, surf );
    slm.t = mean(groupconn);
    slm.tri = surf.tri;
    mask = ones(1,length(slm.t));
    contrast = ones([length(subject_list) 1]); %set length to number of subjects
    slm = SurfStatT(slm, contrast);
    resels = SurfStatResels(slm, mask);
    statthresh = stat_threshold(resels, length(slm.t), 1, slm.df);

    t_thresh = slm.t;
    t_thresh(t_thresh<statthresh)=0;
    figure('visible','off');
%     SurfStatView(t_thresh, surf);
%     SurfStatViewData_LRml(t_thresh, surf, [min(t_thresh) max(t_thresh)], ['Group Connectivity ' area]);
    SurfStatViewData_LRml(t_thresh, surf, [0 23], ['Group Connectivity ' area]);
%     cm=colormap(jet);
%     hold on; SurfStatColormap(cm);
    colorbar('visible', 'off');

    filename = sprintf('%sGroupConn_LR_%s_%s', outdir, session, area);
    if exist('outext', 'var')
        saveas(gcf, [filename outext '.png']); close all;
        fid = fopen([filename outext '.1D'],'w');
    else
        saveas(gcf, [filename '.png']); close all;
        fid = fopen([filename '.1D'],'w');
    end
    fprintf(fid, '%u\n', t_thresh);
    fclose(fid);
    
    %normalize the group t-map
    t_thresh_norm = (t_thresh - min(t_thresh)) / (max(t_thresh)-min(t_thresh));
        
    for h=1:size(groupconn,1)
        %normalize indiv connectivity map
        indivconn_norm = (groupconn(h,:) - min(groupconn(h,:))) / (max(groupconn(h,:))-min(groupconn(h,:)));
        %correlate each individual map with the group map to get a score of individual deviance from the group average
        fitscore_indiv = corr(indivconn_norm', t_thresh_norm');
        %write into fitscore
        fitscore(h, i) = fitscore_indiv;
    end
    tmax(i) = max(t_thresh); %use this to set max for visualization in brainGL
end



