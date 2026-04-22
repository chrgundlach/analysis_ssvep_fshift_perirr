% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'C:\Dropboxdata\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftPerIrr\spm\SPM_06e_Define_Estimate_FlexibleFactorialDesign_precueonly_26Hzonly_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
