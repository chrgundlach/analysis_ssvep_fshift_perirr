%% script to run conversion and preparation of eeg files for spm
% run.pathin      = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\input1\';
% run.pathout     = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata\';

run.pathin      = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\input2\';
run.pathout     = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata_attavereraged\';
run.files       = dir([run.pathin '*.set']);
% run.index       = 295:364;
run.index       = 1:numel(run.files);
run.indfiles    = run.files(run.index);




% List of open inputs
cd(run.pathout) % change to output directory
nrun = numel(run.indfiles); % enter the number of runs here
jobfile = {'C:\Dropboxdata\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftPerIrr\spm\SPM_02_PreprocessConvert_forSPM_FShiftPerIrr_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
    inputs{1, crun} = {[run.pathin run.indfiles(crun).name]}; % Conversion: File Name ? input file
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
