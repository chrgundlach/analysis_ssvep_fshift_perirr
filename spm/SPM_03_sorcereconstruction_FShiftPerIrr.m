run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata\';
run.path        = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spmconvdata_attavereraged\';
run.freq_name   = {'23';'26';'29'};
run.con1         = {'_attended_';'_unattended_'};
run.con2         = {'_offset';'_isolum'};

run.freq2use    = [1 2 3]; %[1 2]
run.con12use     = [1 2]; %[1:7]
run.con22use     = [1 2]; %[1:7]

run.freq_name   = {'23';'26';'29'};
run.con1         = {'_att_collapsed_'};
run.con2         = {'_offset';'_isolum'};

run.freq2use    = [1 2 3]; %[1 2]
run.con12use     = [1]; %[1:7]
run.con22use     = [1 2]; %[1:7]

for i_con1 = 1:numel(run.con12use)
    for i_con2= 1:numel(run.con22use)
        for i_freq = 1:numel(run.freq2use)
            run.files       = dir(sprintf('%sspmeeg*%s%sHz_*%s.mat', ...
                run.path, run.con1{run.con12use(i_con1)}, run.freq_name{run.freq2use(i_freq)}, run.con2{run.con22use(i_con2)}));
            if ~isempty(run.files)

                run.indfiles    = run.files;


                % List of open inputs
                cd(run.path) % change to output directory
                nrun = 1; % enter the number of runs here
                jobfile = {'C:\Dropboxdata\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftPerIrr\spm\SPM_03_sorcereconstruction_FShiftPerIrr_job.m'};
                jobs = repmat(jobfile, 1, nrun);
                inputs = cell(1, nrun);
                for crun = 1:nrun
                    inputs{1, crun} = {run.indfiles.name}'; % files
                end
                spm('defaults', 'EEG');
                spm_jobman('run', jobs, inputs{:});
            end
        end
    end
end


