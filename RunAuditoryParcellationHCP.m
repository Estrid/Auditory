subject_list = importdata(fullfile('example_data', 'subjectlist.txt')); % this should point to a text file containing subject IDs
datadir = fullfile('example_data','HCP','HCP_S1200_glyphsets'); % should point to where individual subject data is
mmpdir = fullfile('example_data','HCP','GroupAvg820_MMP'); % should point to where group average data is
outdir = fullfile('example_data','HCP','AuditoryParcellation'); % should point to where you want results to be saved
session = ('REST1'); % identifies the rs session you wan to use - 'REST1' or 'REST2'
correction = 'FDR'; %indicates the type of correction to apply for the group connectivity maps - can be 'FDR' or 'cluster'

%% run auditory parcellation for individual subjects
for i = 1:length(subject_list)
    try
    AuditoryParcellationInter(subject_list(i), session, datadir, mmpdir, outdir)
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

%% create group probability maps for each subregion
disp('Creating group probability maps...')
datadir = outdir;
outdir = fullfile(outdir, 'GroupResults');
mkdir(outdir);
GroupProbParcellation(subject_list, session, datadir, mmpdir, outdir)

%% compute fit scores for each subject and create group connectivity maps for each subregion
disp('Creating group connectivity maps...')
fitscores = GroupConnParcellation(subject_list, session, datadir, mmpdir, outdir, correction);
csvwrite(fullfile(outdir, 'fitscores.csv'), fitscores);
