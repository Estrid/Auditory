subject_list = importdata('/home/estrid/scripts/bash/subjectlist.txt');

datadir = ('/home/estrid/Data/HCP/HCP_S1200_glyphsets/');
mmpdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');

%% run Broca parcellation

for i = 1:length(subject_list)
    try
    AuditoryParcellationInter(subject_list(i), 'REST1', datadir, mmpdir, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/')
    AuditoryParcellationInter(subject_list(i), 'REST2', datadir, mmpdir, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/')
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

%% create group probability maps for each subregion
subject_list = importdata('/home/estrid/scripts/bash/subjectlist.txt');

GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST1')
GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST2')

% %% do the same as above but separate right and non-right-handed subjects
% subject_list = importdata('/home/estrid/scripts/bash/righthanded.txt');
% GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST1', '_righthanded')
% GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST2', '_righthanded')
% 
% subject_list = importdata('/home/estrid/scripts/bash/nonrighthanded.txt');
% GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST1', '_nonrighthanded')
% GroupProbParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST2', '_nonrighthanded')


%% create group connectivity maps for each subregion

[fitscore1, tmax1] = GroupConnParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST1')
csvwrite('fitscore_REST1.csv', fitscore1)
csvwrite('tmax_REST1.csv', tmax1)
[fitscore2, tmax2] = GroupConnParcellationInter(subject_list, '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', '/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', 'REST2')
csvwrite('fitscore_REST2.csv', fitscore2)
csvwrite('tmax_REST2.csv', tmax2)

%% calculate spatial overlap between RS sessions
% NOTE: ignores areas that don't exist in an individual
ext = '_mmp_partcorr_LR_';

AvgDC100 = zeros(length(subject_list),1);
DC100 = [];
NumNodes100 = [];

for i = 1:length(subject_list)
    try
    [AvgDC, DC, NumNodes] = ParcellationDC('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', subject_list(i), ext);
    AvgDC100(i,1) = AvgDC;
    DC100 = horzcat(DC100,DC);
    NumNodes100 = horzcat(NumNodes100,NumNodes);
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

csvwrite('AvgDC.csv',AvgDC100)
csvwrite('DC100.csv',DC100)
csvwrite('NumNodes100.csv', NumNodes100)

% figure; hist(AvgDC100);
% 
% meanDC = [mean(DC100(1,:),'omitnan') mean(DC100(7,:),'omitnan');...
%     mean(DC100(2,:),'omitnan') mean(DC100(8,:),'omitnan');...
%     mean(DC100(3,:),'omitnan') mean(DC100(9,:),'omitnan');...
%     mean(DC100(4,:),'omitnan') mean(DC100(10,:),'omitnan');...
%     mean(DC100(5,:),'omitnan') mean(DC100(11,:),'omitnan');...
%     mean(DC100(6,:),'omitnan') mean(DC100(12,:),'omitnan');]
% stdDC = [std(DC100(1,:),'omitnan') std(DC100(7,:),'omitnan');...
%     std(DC100(2,:),'omitnan') std(DC100(8,:),'omitnan');...
%     std(DC100(3,:),'omitnan') std(DC100(9,:),'omitnan');...
%     std(DC100(4,:),'omitnan') std(DC100(10,:),'omitnan');...
%     std(DC100(5,:),'omitnan') std(DC100(11,:),'omitnan');...
%     std(DC100(6,:),'omitnan') std(DC100(12,:),'omitnan');]
% 
% figure
% hold on
% bar(meanDC)
% errorbar(meanDC,stdDC,'.')

%% calculate spatial overlap between hemispheres

ext = '_mmp_partcorr_LR_REST1.1D';
AvgDC100 = zeros(length(subject_list),1);
DC100 = [];

for i = 1:length(subject_list)
    try
    [AvgDC, DC] = ParcellationDC_hemis('/home/estrid/Data/HCP/MMP_indiv/9areas_nonon/', subject_list(i), ext);
    AvgDC100(i,1) = AvgDC;
    DC100 = horzcat(DC100,DC);
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

csvwrite('AvgDC_hemis.csv',AvgDC100)
csvwrite('DC100_hemis.csv',DC100)
