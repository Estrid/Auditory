subject_list = importdata('/home/estrid/scripts/subjectlist.txt');

outdir = ('/home/estrid/Data/HCP/MMP_indiv/');
datadir = ('/home/estrid/pamina_mount/HCP/HCP_S1200_glyphsets/');
mmpdir = ('/home/estrid/Data/HCP/GroupAvg820_MMP/');

%% run Broca parcellation

for i = 1:length(subject_list)
    try
    AuditoryParcellation(subject_list(i), 1, datadir, mmpdir, outdir)
    AuditoryParcellation(subject_list(i), 2, datadir, mmpdir, outdir)
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

