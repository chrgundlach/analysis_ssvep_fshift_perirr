% script to run source reconstruction and create output image for specified
% time and frequency

% run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata\';
run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata_attavereraged\';
run.freqs_all   = [23 23;26 26; 29 29]; % upper/lower limit of frequencies to extract
run.freq_name   = {'23Hz';'26Hz';'29Hz'};
% run.con         = {'_attended';'_unattended'};
run.con         = {'att_collapsed'};
run.D_index     = 1; % use specified sourceconstruction file
run.timewindows = [0 1000; 1500 2500];

run.freq2use    = [1:3]; % [1:2]
run.con2use     = [1]; % [1:2]
run.time2use    = [1]; % [1:2]

for i_freq = 1:numel(run.freq2use)
    for i_con = 1:numel(run.con2use)
        for i_time = 1:numel(run.time2use)
            run.files{i_freq,i_con,i_time}      = dir(sprintf('%sspmeeg*%s_%s*.mat',run.path,run.con{run.con2use(i_con)},run.freq_name{run.freq2use(i_freq)}));
            run.freqs{i_freq,i_con,i_time}      = run.freqs_all(i_freq,:);
            run.times{i_freq,i_con,i_time}      = run.timewindows(i_time,:);
        end
    end
end

%%
% List of open inputs
cd(run.path) % change to output directory
nrun = numel(run.files); % enter the number of runs here
jobfile = {'C:\Dropboxdata\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftPerIrr\spm\SPM_04_source_images_FShiftPerIrr_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);
for crun = 1:nrun
    inputs{1, crun} = {run.files{crun}.name}'; % files
    inputs{2, crun} = run.D_index; % D
    inputs{3, crun} = run.times{crun}; % times
    inputs{4, crun} = run.freqs{crun}; % files
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
