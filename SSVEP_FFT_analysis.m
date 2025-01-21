%% analysis script for FShiftPerirr
% SSVEP analysis based on FFT 


clearvars
%% parameter
F.PathInEEG = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\epoch\';
F.PathInBehavior = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\behavior\';
F.PathOut = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\fft\';

F.trigger               = {[10 ]; ... %RDK1 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                            [20 ]; ... %RDK2 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                            [30 ]; ... %RDK1 attended; RDK1 and RDK3 colors in periphery peri attended + irrelevant
                            [40 ]; ... %RDK2 attended; RDK2 and RDK3 colors in periphery peri attended + irrelevant
                            [50 ]; ... %RDK1 attended; RDK2 and RDK3 colors in periphery peri unattended + irrelevant
                            [60 ]};  %RDK2 attended; RDK1 and RDK3 colors in periphery peri unattended + irrelevant

F.FFT_timewins = {[-1 0]; [0.5 1.5]}; % relevant time windows in s
F.FFT_freqres           = 2^12; % frequency resolution
F.SSVEPfreqrange = [-0.1 +0.1]; % frequency range of SSVEPs analyzed


F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% F.sub2use               = [6:13 15:18];%:53;
F.sub2use               = [22:52];%
F.sub2use               = [33:52];%


for i_sub = 1:numel(F.sub2use)
    %% load files
    % load eeg
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);
    % pop_eegplot(EEG,1,1,1,1)
    % eeglab redraw % load eeg gui + current data

    % load behavior(loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));

    %% do csd transform
    if  i_sub == 1 % calculate CSD matrix
        try
            CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
        catch
            CSD.chanmat=ExtractMontage('C:\Users\EEG\Documents\MATLAB\christopher\general_functions\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
        end
        [CSD.G,CSD.H] = GetGH(CSD.chanmat);
    end
    fprintf(1,'\n###\ncalculating CSD transform\n###\n')
    for i_tr = 1:EEG.trials
        % csd of raw data
        EEG.data(:,:,i_tr)= CSDTransform(EEG.data(:,:,i_tr), CSD.G, CSD.H);
    end

    %% do FFT transform
    % bookkeeping
    Results.electrodes = EEG.chanlocs; % save electrode positions
    Results.RDK = behavior.RDK; % save
    Results.all_trialnum = EEG.trials;
    Results.con_trigger = F.trigger;
    Results.con_trialnum = cellfun(@(x) sum(ismember([EEG.event.type],x)), F.trigger); % count trial number for each condition

    Results.FFT.res = F.FFT_freqres;
    Results.FFT.timewin = F.FFT_timewins;
    Results.FFT.data_evo= nan(Results.FFT.res, EEG.nbchan,numel(F.trigger),numel(Results.FFT.timewin));

    % index trials
    t.trialindex = cellfun(@(x) [EEG.event([ismember([EEG.event.type],x)]).epoch], F.trigger,'UniformOutput',false);

    % do FFT transform for evoked data
    for i_win = 1:numel(Results.FFT.timewin)
        % select data for specific time window
        EEGt = pop_select(EEG,'time',Results.FFT.timewin{i_win});
        % detrend
        EEGt = eegF_Detrend(EEGt,[]);

        % loop across conditions
        for i_con = 1:numel(Results.con_trigger)
            t.data = detrend(mean(EEGt.data(:,:,t.trialindex{i_con}),3)')'; % average across trials
            Results.FFT.data_evo(:,:,i_con,i_win) = squeeze(abs(fft(t.data,Results.FFT.res,2))*2/size(t.data,2))'; % fft transform + normalization to get amplitude spectrum with unit muV/cm²
        end
    end
    Results.FFT.freqs = ((0:size(Results.FFT.data_evo,1)-1)/size(Results.FFT.data_evo,1)) * EEGt.srate;

    % figure; plot(Results.FFT.freqs, Results.FFT.data_evo(:,29,1,1)) % plot all frequencies at Oz for con 1 and time window 1
    % figure; plot(Results.FFT.freqs, mean(Results.FFT.data_evo(:,29,:,1),3))

    %% extract SSVEP amplitudes
    Results.FFT.SSVEPamps = nan(numel(Results.RDK.RDK),numel(Results.electrodes),numel(F.trigger),numel(Results.FFT.timewin));

    % loop across RDK
    for i_rdk = 1:numel(Results.RDK.RDK)
        % what frequency?
        t.freq = Results.RDK.RDK(i_rdk).freq + F.SSVEPfreqrange;
        % index of frequencies in spectra data
        t.freq_i = dsearchn(Results.FFT.freqs',t.freq');
        % extract SSVEP amplitude in frequency range
        Results.FFT.SSVEPamps(i_rdk,:,:,:) = mean(Results.FFT.data_evo(t.freq_i(1):t.freq_i(2),:,:,:),1);
    end

    %% save all of it
    Results.savetime=datestr(now);
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    fprintf(1,'|| saving file ||  %s\\VP%s_tfa.mat ||\n', ...
        F.PathOut,F.subjects{F.sub2use(i_sub)})
    save(fullfile(F.PathOut, sprintf('VP%s_tfa.mat',F.subjects{F.sub2use(i_sub)})), 'Results')
end






%% additional functions
% plot electrode positions by label and number
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% plot the actual data
pl.RDKi = 1;
pl.data = squeeze(mean(Results.FFT.SSVEPamps(pl.RDKi,:,:,1),3));
pl.clim = [0 max(pl.data)]; % range of data to be plotted

figure;
topoplot(pl.data, Results.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',viridis,...
    'whitebk','on');
title(sprintf('RDK %1.0f | SSVEP freq %1.0f Hz',pl.RDKi,Results.RDK.RDK(pl.RDKi).freq))
colorbar
