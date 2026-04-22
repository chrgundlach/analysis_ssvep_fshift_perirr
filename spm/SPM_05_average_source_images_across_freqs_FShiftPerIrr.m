% script to average the source images for the two RDK frequencies (frequencies shouldn't matter)
clear run
% run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata\';
% run.pathout     = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\EEG\SPM\spmconvdata_freqaveraged';
run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata_attavereraged\';
run.pathout     = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata_attavereraged_freqaveraged';
t.files         = dir(sprintf('%s*.nii',run.path));
run.subs        = unique(regexp({t.files.name}, 'VP\d+', 'match', 'once'));
% run.con         = {'_attended';'_unattended'};
run.con         = {'att_collapsed'};
run.times       = {'_t0_1000';'_t1500_2500'};

run.con2        = {'_isolum';'_offset'};



run.subs2use    = [1:50]; % [1:50]
run.con2use     = [1]; % [1:2]
run.time2use    = [1]; % [1:2]

for i_sub = 1:numel(run.subs2use)
    for i_con = 1:numel(run.con2use)
        for i_time = 1:numel(run.time2use)
            run.files{i_sub,i_con,i_time}       = dir( ...
                sprintf('%s*%s*%s*%s*.nii', ...
                run.path, run.subs{run.subs2use(i_sub)}, run.con{run.con2use(i_con)},  run.times{run.time2use(i_time)}) ...
                );
            % check for offset or isoluminance
            t.con2name = run.con2{cellfun(@(x) contains(run.files{i_sub,i_con,i_time}(1).name,x),run.con2)};
            run.nameout{i_sub,i_con,i_time} = sprintf('spmeeg_%s%s_allSSVEPs%s%s.nii', ...
                run.subs{run.subs2use(i_sub)}, run.con{run.con2use(i_con)},t.con2name,run.times{run.time2use(i_time)});
        end
    end
end

%%
% List of open inputs
cd(run.path) % change to output directory
nrun = numel(run.files); % enter the number of runs here
jobfile = {'C:\Dropboxdata\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftPerIrr\spm\SPM_05_average_source_images_across_freqs_FShiftPerIrr_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = {run.files{crun}.name}'; % files
    inputs{2, crun} = run.nameout{crun}; % name out
    inputs{3, crun} = {run.pathout}; % path out
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
