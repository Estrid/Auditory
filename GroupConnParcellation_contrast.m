function[tmax] = GroupConnParcellation_contrast(subject_list, datadir, outdir, session, outext)

subject_list = importdata('/home/estrid/scripts/bash/subjectlist.txt');
datadir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
outdir = ('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/');
session = ('REST1');

% subject_list = importdata(subject_list);
ext = '_mmp_meanconn_LR_'; %extension of files to import

surf = SurfStatReadSurf( {...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated'),...
('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated')} );

areas = {'PAC', 'SAC', 'RI', 'A5', 'STGa', 'STSda', 'STSdp', 'STSva', 'STSvp'};

tmax = zeros(size(areas));

for i = 1:length(areas)
    area = areas(i);
    area = char(area);
    groupconn_L = zeros(size(surf.coord,2), length(subject_list));
    groupconn_R = zeros(size(surf.coord,2), length(subject_list));
    errors = 0;
    for h = 1:length(subject_list)
        subjectID = num2str(subject_list(h));
        try
        conn_L = importdata([datadir subjectID ext session '_' area '_L.1D']);
        groupconn_L(:,h) = conn_L;
        conn_R = importdata([datadir subjectID ext session '_' area '_R.1D']);
        groupconn_R(:,h) = conn_R;
        catch
            fprintf('error in subject # %u\n', subject_list(h))
            errors = errors+1;
            continue
        end
        fprintf('subject # %u\n',subject_list(h))
    end
    
    groupconn = horzcat(groupconn_L, groupconn_R);
    groupconn = groupconn';
    left = vertcat(ones([length(subject_list) 1]), zeros([length(subject_list) 1]));
    right = vertcat(zeros([length(subject_list) 1]), ones([length(subject_list) 1]));
    hemi = horzcat(left, right);
    hemi = term(hemi);
    M = 1 + hemi;
    slm = SurfStatLinMod(groupconn, M, surf );
    
    slm.t = mean(groupconn);
    slm.tri = surf.tri;
    mask = ones(1,length(slm.t));
    
    %left - right
    contrast = hemi.hemi1 - hemi.hemi2;
    slm = SurfStatT(slm, contrast);
    resels = SurfStatResels(slm, mask);
    statthresh = stat_threshold(resels, length(slm.t), 1, slm.df);
    t_thresh = slm.t;
    t_thresh(t_thresh<statthresh)=0;
    
    LR = t_thresh;
    LR(LR>0)=-1;
    
%     figure('visible','off');
%     SurfStatViewData_LRml(slm.t, surf, [min(slm.t) max(slm.t)], ['Group Connectivity Left-Right ' area]);
%     saveas(gcf, [outdir 'GroupConn_LR_' session '_' area '_contrastL-R.png']); close all;
%     fid = fopen([outdir 'GroupConn_LR_' session '_' area '_contrastL-R.1D'],'w');
%     fprintf(fid, '%u\n', t_thresh);
%     fclose(fid);
    
%     figure('visible','off');
%     SurfStatViewData_LRml(t_thresh, surf, [min(t_thresh) max(t_thresh)], ['Group Connectivity Left-Right ' area]);
%     saveas(gcf, [outdir 'GroupConn_LR_' session '_' area '_contrastL-R_tthresh.png']); close all;
%     fid = fopen([outdir 'GroupConn_LR_' session '_' area '_contrastL-R_tthresh.1D'],'w');
%     fprintf(fid, '%u\n', t_thresh);
%     fclose(fid);
    
    %right - left
    contrast = hemi.hemi2 - hemi.hemi1;
    slm = SurfStatT(slm, contrast);
    resels = SurfStatResels(slm, mask);
    statthresh = stat_threshold(resels, length(slm.t), 1, slm.df);
    t_thresh = slm.t;
    t_thresh(t_thresh<statthresh)=0;
    
    RL = t_thresh;
    RL(RL>0)=1;
    
%     figure('visible','off');
%     SurfStatViewData_LRml(slm.t, surf, [min(slm.t) max(slm.t)], ['Group Connectivity Right-Left ' area]);
%     saveas(gcf, [outdir 'GroupConn_LR_' session '_' area '_contrastR-L.png']); close all;
%     fid = fopen([outdir 'GroupConn_LR_' session '_' area '_contrastR-L.1D'],'w');
%     fprintf(fid, '%u\n', t_thresh);
%     fclose(fid);
    
%     figure('visible','off');
%     SurfStatViewData_LRml(t_thresh, surf, [min(t_thresh) max(t_thresh)], ['Group Connectivity Right-Left ' area]);
%     saveas(gcf, [outdir 'GroupConn_LR_' session '_' area '_contrastR-L_tthresh.png']); close all;
%     fid = fopen([outdir 'GroupConn_LR_' session '_' area '_contrastR-L_tthresh.1D'],'w');
%     fprintf(fid, '%u\n', t_thresh);
%     fclose(fid);

    vizcontrast = LR+RL;
    figure;%('visible','off');
    SurfStatViewData_LRml(vizcontrast, surf, [-1.1 1.1], ['Contrast ' area]);
    cm=[zeros(1,3);
        zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.8;
        ones(127,1)    (0:126)'/127   zeros(127,1)];
    hold on; SurfStatColormap(cm);
    colorbar('visible', 'off');
    saveas(gcf, [outdir 'GroupConn_LR_' session '_' area '_contrastviz.png']); close all;
    fid = fopen([outdir 'GroupConn_LR_' session '_' area '_contrastviz.1D'],'w');
    fprintf(fid, '%u\n', t_thresh);
    fclose(fid);
    
    tmax(i) = max(t_thresh); %use this to set max for visualization in brainGL
end



%% Or import thresholded t-maps and subtract

% session = 'REST1_';
% 
% surf = SurfStatReadSurf( {...
% ('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/lh.very_inflated'),...
% ('/home/estrid/Data/HCP/HCP_S900_GroupAvg_v1/rh.very_inflated')} );
% 
% areas = {'PAC', 'SAC', 'RI', 'A5', 'STGa', 'STSda', 'STSdp', 'STSva', 'STSvp'};
% 
% % import left and right group connectivity maps
% for i = 1:length(areas)
%     area = areas(i);
%     area = char(area);
%     left = importdata([filename session area '_L.1D']);
%     right = importdata([filename session area '_R.1D']);
%     t_contrastLR = left - right;
%     t_contrastRL = right - left;
%     
% %     norm = 2*(t_contrastLR - min(t_contrastLR))/(max(t_contrastLR) - min(t_contrastLR)) - 1;
% %     figure;%('visible','off');
% %     SurfStatViewData_LRml(norm, surf, [min(norm) max(norm)], ['Group Connectivity Contrast L-R ' area]);
% %     saveas(gcf, [filename session area '_contrastL-Rnorm.png']); close all;
% %     fid = fopen([filename session area '_contrastL-Rnorm.1D'],'w');
% %     fprintf(fid, '%u\n', norm);
% %     fclose(fid);
%     
%     figure;%('visible','off');
%     SurfStatViewData_LRml(t_contrastLR, surf, [min(t_contrastLR) max(t_contrastLR)], ['Group Connectivity Contrast L-R ' area]);
%     saveas(gcf, [filename session area '_contrastL-R.png']); close all;
%     fid = fopen([filename session area '_contrastL-R.1D'],'w');
%     fprintf(fid, '%u\n', t_contrastLR);
%     fclose(fid);
%     
%     figure;%('visible','off');
%     SurfStatViewData_LRml(t_contrastRL, surf, [min(t_contrastRL) max(t_contrastRL)], ['Group Connectivity Contrast R-L ' area]);
%     saveas(gcf, [filename session area '_contrastR-L.png']); close all;
%     fid = fopen([filename session area '_contrastR-L.1D'],'w');
%     fprintf(fid, '%u\n', t_contrastRL);
%     fclose(fid);
% end