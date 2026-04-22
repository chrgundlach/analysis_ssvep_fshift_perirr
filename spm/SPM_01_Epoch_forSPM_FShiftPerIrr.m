function [] = SPM_01_Epoch_forSPM_FShiftPerIrr()
%SPM_EPOCH_FORSPM prepare data for spm source reconstruction
%   loads different subjects
%   epochs data for relevant events and averages those
%   filters data
%   prepares EEG files for SPM source reconstruction
%   creates set file for each

%% parameters
clearvars
F.Pathlocal             = 'E:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\';


F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch\');
F.PathInBehavior        = fullfile(F.Pathlocal, 'behavior\');
F.PathOut               = fullfile(F.Pathlocal, 'eeg\SPM\input1\'); % with FWHM 0.5
F.PathOut               = fullfile(F.Pathlocal, 'eeg\SPM\input2\'); % with FWHM 0.5 % ignore ATTENTION


F.Subjects          = arrayfun(@(x) sprintf('%02.0f',x),1:60,'UniformOutput',false)';
F.Sub2Use           = [1:13 15:41 43:52];
% F.Sub2Use           = [1:3];

% keep track of peripheral conditions
F.ConTrigger        = {[10 ]; ... %RDK1 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                        [20 ]; ... %RDK2 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                        [30 ]; ... %RDK1 attended; RDK1 and RDK3 colors in periphery peri attended + irrelevant
                        [40 ]; ... %RDK2 attended; RDK2 and RDK3 colors in periphery peri attended + irrelevant
                        [50 ]; ... %RDK1 attended; RDK2 and RDK3 colors in periphery peri unattended + irrelevant
                        [60 ]};  %RDK2 attended; RDK1 and RDK3 colors in periphery peri unattended + irrelevant
F.CondNames         = {'RDK1_att_periRDK1_2'; 'RDK2_att_periRDK1_2';
                        'RDK1_att_periRDK1_3';'RDK2_att_periRDK1_3';
                        'RDK1_att_periRDK2_3';'RDK2_att_periRDK2_3'};
F.CondNames2         = [{'attended'; 'unattended';'attended'; 'unattended';'attended'; 'unattended'}, ...
                        {'unattended';'attended';'unattended';'attended';'unattended';'attended'}];

% ignore peripheral conditions
F.ConTrigger        = {[10 30 50]; ... %RDK1 attended; 
                        [20 40 60]};  %RDK2 attended;
F.CondNames         = {'RDK1_att_all_peri';'RDK2_att_all_peri'};
F.CondNames2         = [{'attended'; 'unattended'}, ...
                        {'unattended';'attended'}];
F.SubCons           = {[1:13 15:21]; ... % offset
                        [22:41 43:52]}; % isolum
F.SubConNames       = {'offset';'isolum'};

% ignore peripheral conditions and ignore attention
F.ConTrigger        = {[10 30 50 20 40 60]; ... %RDK1 attended; 
                        };  %RDK2 attended;
F.CondNames         = {'peri_collapsed';'peri_collapsed'};
F.CondNames2         = [{'att_collapsed', 'att_collapsed'}];
F.SubCons           = {[1:13 15:21]; ... % offset
                        [22:41 43:52]}; % isolum
F.SubConNames       = {'offset';'isolum'};


F.FilterFreqs       = [2 6];

F.TimeOfInterest    = [-1 2];


%% actual preparation of data for SPM
step=1;
h.wb = waitbar(0,'Preparing data for SPM...please wait...');
for i_sub = 1:numel(F.Sub2Use)
    fprintf('...loading data VP%s\n',F.Subjects{F.Sub2Use(i_sub),:})
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.Subjects{F.Sub2Use(i_sub)}),'filepath',F.PathInEEG);
    %pop_eegplot(EEG,1,1,1)

    % behavior (loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.Subjects{F.Sub2Use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));

    % pull out RDK frequencies
    t.rdkfreqs = [behavior.RDK.RDK(1:2).freq];
    
    % loop for triggers
    for i_tr = 1:numel(F.ConTrigger)
        EEG_ep = pop_selectevent( EEG, 'type',F.ConTrigger{i_tr} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        %pop_eegplot(EEG_ep,1,1,1)
        
        %detrend
        EEG_epd = EEG_ep;
        for i_trial = 1:EEG_ep.trials
            EEG_epd.data(:,:,i_trial)=detrend(EEG_ep.data(:,:,i_trial)')';
        end
        %pop_eegplot(EEG_epd,1,1,1)
        
        for i_RDK = 1:numel(t.rdkfreqs)
            waitbar(step/(numel(t.rdkfreqs)*numel(F.ConTrigger)*numel(F.Sub2Use)))
            EEG_epdm = EEG_epd;
            EEG_epdm = pop_select( EEG_epdm,'trial',[1 2] );
            EEG_epdm.epochdescription=F.CondNames{i_tr};
            EEG_epdm=eeg_checkset(EEG_epdm);
            
            EEG_epdm.data(:,:,1) = mean((EEG_epd.data),3);
            % pop_eegplot(EEG_epdm,1,1,1)
            % filter
            EEG_epdmf = pop_eegfiltnew(EEG_epdm, t.rdkfreqs(i_RDK)-.5, t.rdkfreqs(i_RDK)+.5, 8*EEG_epdm.srate, 0, [], 0);
            % pop_eegplot(EEG_epdmf,1,1,1)
            % select time
            EEG_epdmft = pop_select(EEG_epdmf, 'time', F.TimeOfInterest);
            EEG_epdmft = pop_select( EEG_epdmft,'trial',[1] );
            EEG_epdmft.xmin = min(EEG_epdmft.times)/1000;  EEG_epdmft.xmax = max(EEG_epdmft.times)/1000;
            % save
            [FreqInt,FreqFract]=deal(fix(t.rdkfreqs(i_RDK)), t.rdkfreqs(i_RDK)-fix(t.rdkfreqs(i_RDK)));
            tt=sprintf('%1.1f',FreqFract);
            % FileName_New=sprintf('%s_RDK%1.0f_%1.0fp%sHz_%s',EEG.filename(1:end-4),i_RDK,FreqInt,tt(end),F.CondNames{i_tr});
            FileName_New=sprintf('%s_%s_%1.0fHz_RDK%1.0f_%s_%s', ...
                EEG.filename(1:end-4), ...
                F.CondNames2{i_tr,i_RDK}, ...
                FreqInt, ...
                i_RDK, ...
                F.CondNames{i_tr}, ...
                F.SubConNames{cellfun(@(x) any(ismember(x,F.Sub2Use(i_sub))), F.SubCons)});
            if ~exist(F.PathOut); mkdir(F.PathOut);end
            EEG_epdmft = pop_saveset(EEG_epdmft, 'filename',FileName_New,'filepath',F.PathOut);
            
            step=step+1;
        end
    end
end
close(h.wb)

