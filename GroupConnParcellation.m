function[fitscore] = GroupConnParcellation(subject_list, session, datadir, mmpdir, outdir, correction)

surf = SurfStatReadSurf( {...
[mmpdir 'lh.very_inflated'],...
[mmpdir 'rh.very_inflated']} );

mask = importdata([mmpdir 'medialwallmask.1D']);
mask = mask';
mask = logical(mask);

areas = {'PAC', 'SAC', 'RI', 'A5', 'STGa', 'STSda', 'STSdp', 'STSva', 'STSvp'};

fitscore = zeros(length(subject_list), length(areas));

for i = 1:length(areas)
    area = areas(i);
    area = char(area);
    groupconn_L = zeros(size(surf.coord,2), length(subject_list));
    groupconn_R = zeros(size(surf.coord,2), length(subject_list));
    errors = 0;
    for h = 1:length(subject_list)
        subjectID = num2str(subject_list(h));
        try
        conn_L = importdata([datadir subjectID '_AuditoryParcellation_' session '_meanconn_' area '_L.1D']);
        groupconn_L(:,h) = conn_L;
        catch
            fprintf('error in subject # %u - left %s does not exist\n', subject_list(h), area)
            errors = errors+1;
            continue
        end
        try
        conn_R = importdata([datadir subjectID '_AuditoryParcellation_' session '_meanconn_' area '_R.1D']);
        groupconn_R(:,h) = conn_R;
        catch
            fprintf('error in subject # %u - right %s does not exist\n', subject_list(h), area)
            errors = errors+1;
            continue
        end
        fprintf('subject # %u\n',subject_list(h))
    end
    
    for j= 1:2
        if j==1
            groupconn = groupconn_L;
        elseif j==2
            groupconn = groupconn_R;
        end
        groupconn = groupconn';
        slm = SurfStatLinMod(groupconn, 1, surf );
        contrast = ones([length(subject_list) 1]);
        slm = SurfStatT(slm, contrast);
        [~, ~] = SurfStatP(slm, mask, 0.05);
        
        if isequal(correction, 'cluster') % applies cluster correction
            t_thresh = slm.t;
            resels = SurfStatResels(slm, mask);
            statthresh = stat_threshold(resels, length(slm.t), 1, slm.df);
            t_thresh(t_thresh<statthresh)=0;
        elseif isequal(correction, 'FDR') % applies FDR correction
            qval = SurfStatQ(slm, mask);
            t_thresh = slm.t;
            t_thresh(qval.Q>0.05) = 0;
        end
        
        % viz and save
        figure('visible','off');
        SurfStatViewData_LRml(t_thresh, surf, [0 23], ['Group Connectivity ' correction area]);
        colorbar('visible', 'off');

        filename = sprintf('%sGroupConn_%s_%s_%s', outdir, session, area, correction);
        saveas(gcf, [filename '.png']); close all;
        fid = fopen([filename '.1D'],'w');
        fprintf(fid, '%u\n', t_thresh);
        fclose(fid);
        
        %normalize group t-map
        t_thresh_norm = (t_thresh - min(t_thresh)) / (max(t_thresh)-min(t_thresh));
        
        %calculate fit score
        for h=1:size(groupconn,1)
            %normalize indiv connectivity map
            indivconn_norm = (groupconn(h,:) - min(groupconn(h,:))) / (max(groupconn(h,:))-min(groupconn(h,:)));
            %correlate each individual map with the group map to get a score of individual deviance from the group average
            fitscore_indiv = corr(indivconn_norm', t_thresh_norm');
            %write into fitscore
            fitscore(h, i) = fitscore_indiv;
        end
    end
    
    %% compute contrast between left and right
    
    groupconn = horzcat(groupconn_L, groupconn_R);
    groupconn = groupconn';
    left = vertcat(ones([length(subject_list) 1]), zeros([length(subject_list) 1]));
    right = vertcat(zeros([length(subject_list) 1]), ones([length(subject_list) 1]));
    hemi = horzcat(left, right);
    hemi = term(hemi);
    M = 1 + hemi;
    slm = SurfStatLinMod(groupconn, M, surf);

    contrast1 = hemi.hemi1 - hemi.hemi2;
    contrast2 = hemi.hemi2 - hemi.hemi1;
    
    for j= 1:2
        if j==1
            contrast = contrast1; %left - right
        elseif j==2
            contrast = contrast2; %right - left
        end

        slm = SurfStatT(slm, contrast);
        if isequal(correction, 'cluster')
            resels = SurfStatResels(slm, mask);
            statthresh = stat_threshold(resels, length(slm.t), 1, slm.df);
            t_thresh = slm.t;
            t_thresh(t_thresh<statthresh)=0;
        elseif isequal(correction, 'FDR')
            qval = SurfStatQ(slm, mask);
            t_thresh = slm.t;
            t_thresh(qval.Q>0.05) = 0;
        end

        if j==1
            LR = t_thresh;
            LR(LR>0)=-1; %left - right
        elseif j==2
            RL = t_thresh;
            RL(RL>0)=1; %right - left
        end
    end

    vizcontrast = LR+RL;
    figure;%('visible','off');
    SurfStatViewData_LRml(vizcontrast, surf, [-1.1 1.1], ['Contrast ' correction area]);
    cm=[zeros(1,3);
        zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.8;
        ones(127,1)    (0:126)'/127   zeros(127,1)];
    hold on; SurfStatColormap(cm);
    colorbar('visible', 'off');
    filename = sprintf('%sGroupConn_contrast_%s_%s_%s', outdir, session, area, correction);
    saveas(gcf, [filename '.png']); close all;
    fid = fopen([filename '.1D'],'w');
    fprintf(fid, '%u\n', t_thresh);
    fclose(fid);

end


