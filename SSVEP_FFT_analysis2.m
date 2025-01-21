%% script to plot and extract data
%% parameters
clearvars
F.PathInEEG             = 'C:\Users\psy05cvd\Desktop\fft'; % with FWHM 0.5
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\fft'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% participant 42 has lower trial number
F.Subs2use              = [22:41 43:52];
% F.Subs2use              = [22:32];
                        

F.SSVEP_Freqs           = [14 17 20 23 26]; % sub 22 onwars
F.RDK_pos               = [0 -255 255];
F.RDK_pos_label         = {'center';'left';'right'};

F.conlabel_att = {'att RDK1';'att RDK2'; 'att RDK1';'att RDK2'; 'att RDK1';'att RDK2'};
F.conlabel_pres_periRDK = {'RDK3 RDK4';'RDK3 RDK4'; 'RDK3 RDK5'; 'RDK4 RDK5';'RDK4 RDK5'; 'RDK3 RDK5'};
F.conlabel_pres_periRDK1 = {'RDK3';'RDK3'; 'RDK3'; 'RDK4';'RDK4'; 'RDK3'};
F.conlabel_pres_periRDK2 = {'RDK4';'RDK4'; 'RDK5'; 'RDK5';'RDK5'; 'RDK5'};
F.conlabel_att_periRDK1 = {'att';'unatt'; 'att'; 'att';'unatt'; 'unatt'};
F.conlabel_att_periRDK2 = {'unatt';'att'; 'irrel'; 'irrel';'irrel'; 'irrel'};
F.conRDKpresent = logical([1 1 1 1 0; 1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 1 0 1 1; 1 1 1 0 1]);
F.conRDKcolattended = [1 0 1 0 2; 0 1 0 1 2; 1 0 1 0 2; 0 1 0 1 2; 1 0 1 0 2; 0 1 0 1 2];
F.conRDKpresent_label = repmat({'not presented'},size(F.conRDKpresent));
F.conRDKpresent_label(F.conRDKpresent)={'presented'};
F.conRDKcolattended_label = repmat({'not attended'},size(F.conRDKcolattended));
F.conRDKcolattended_label(F.conRDKcolattended==1)={'attended'};
F.conRDKcolattended_label(F.conRDKcolattended==2)={'irrelevant'};

% an additional category [manually]
F.conRDKcolattended_label2 = F.conRDKcolattended_label;
[F.conRDKcolattended_label2{[13 20]}] = deal('attended + not attended');
[F.conRDKcolattended_label2{[14 19]}] = deal('not attended + attended');
[F.conRDKcolattended_label2{[15 22]}] = deal('attended + irrelevant');
[F.conRDKcolattended_label2{[18 23]}] = deal('not attended + irrelevant');
[F.conRDKcolattended_label2{[27 28]}] = deal('irrelevant + attended');
[F.conRDKcolattended_label2{[29 30]}] = deal('irrelevant + not attended');

%% load data for all participants
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_results.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})

    temp.results = open(fullfile(F.PathInEEG,sprintf('VP%s_tfa.mat',F.Subs{F.Subs2use(i_sub)})));

    % preallocate memory
    if i_sub == 1
        RESULTS.electrodes = temp.results.Results.electrodes;
        RESULTS.con_trialnum = nan([numel(temp.results.Results.con_trialnum) numel(F.Subs2use)]);
        RESULTS.fftdata_evo = nan([size(temp.results.Results.FFT.data_evo),numel(F.Subs2use)]);
        RESULTS.ffttimewin = temp.results.Results.FFT.timewin;
        RESULTS.fftfreqs = temp.results.Results.FFT.freqs;
        RESULTS.SSVEPamps = nan([size(temp.results.Results.FFT.SSVEPamps),numel(F.Subs2use)]);
    end

    % fill data
    RESULTS.fftdata_evo(:,:,:,:,i_sub) = temp.results.Results.FFT.data_evo;
    RESULTS.SSVEPamps(:,:,:,:,i_sub) = temp.results.Results.FFT.SSVEPamps;
    RESULTS.RDK(i_sub) = temp.results.Results.RDK;

    % add position of RDK to data
    t.label = cellfun(@(x) F.RDK_pos_label(x(1)==F.RDK_pos),{RESULTS.RDK(i_sub).RDK.centershift});
    [RESULTS.RDK(i_sub).RDK.poslabel] = deal(t.label{:});
    
    clear temp

end

%% plot spectra for central stimuli
% large center as in tango | periphery: central and lateral 
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({RESULTS.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


pl.time2plot = [1]; % which window of baseline [1]
pl.sub2plot = 1:numel(F.Subs2use);
% pl.sub2plot = 1:11;

pl.freq2plot=F.SSVEP_Freqs(4);

% extract data
pl.data_evo = nan(size(RESULTS.fftdata_evo,1), numel(pl.sub2plot)); 
t.elidx = [];
for i_sub = 1:numel(pl.sub2plot)
    % index position of frequency
    t.rdkidx = [RESULTS.RDK(pl.sub2plot(i_sub)).RDK.freq]==pl.freq2plot;
    t.posidx = RESULTS.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).poslabel;
    t.elidx(:,i_sub) = strcmp(pl.elec2plot(:,2), t.posidx);
    % index conditions for which the RDK is shown
    t.conidx = any(t.rdkidx & F.conRDKpresent,2);

    % extract data
    pl.data_evo(:,i_sub) = mean(RESULTS.fftdata_evo(:,pl.elec2plot_i{find(t.elidx(:,i_sub))},t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[2,3,4]);
end


% plotting
figure;
set(gcf,'Position',[100 100 600 400],'PaperPositionMode','auto')
plot(RESULTS.fftfreqs,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(RESULTS.fftfreqs,mean(pl.data_evo,2),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f | FOI = %1.1f Hz', ...
    numel(pl.sub2plot), pl.freq2plot),'Interpreter','none')
vline(pl.freq2plot,'k:')
% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(1:64,RESULTS(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i(any(t.elidx'))),1)),'o','r',4,1});

%% plot Grand Mean FFT data | topoplot for positions (as frequencies are random)
pl.time2plot = [1];
pl.pos2plot='center';
% pl.pos2plot='right';
% pl.pos2plot='left';
pl.sub2sel = 1:numel(F.Subs2use);
% pl.sub2sel = find(F.Subs2use<22); % luminance offset
pl.sub2sel = find(F.Subs2use>21); % isoluminant to background

pl.sub2plot = pl.sub2sel( ...
    cellfun(@(x) any(strcmp({x.poslabel},pl.pos2plot)), ... % index for which subject RDKs were presented at pl.pos2plot
    {RESULTS.RDK(pl.sub2sel).RDK}));

for i_sub=1:numel(pl.sub2plot)
    % which RDK is relevant
    t.rdkidx = strcmp({RESULTS.RDK(pl.sub2plot(i_sub)).RDK.poslabel},pl.pos2plot);
    t.rdkidx2 = find(t.rdkidx );
    % preallocate data
    if i_sub == 1
        t.data = nan(numel(RESULTS.electrodes),numel(t.rdkidx2) ,numel(pl.sub2plot));
    end

    for i_rdk = 1:numel(t.rdkidx2) % loop across RDKs
        % which condition is relevant?
        t.conidx = F.conRDKpresent(:,t.rdkidx2(i_rdk));
        t.data(:,i_rdk,i_sub) = mean(RESULTS.SSVEPamps(t.rdkidx2(i_rdk),:,t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[3,4]);
    end
end

pl.data_evo = mean(t.data,[2,3]);



figure;
set(gcf,'Position',[100 100 300 300],'PaperPositionMode','auto')

pl.clim = [0 max(pl.data_evo)];
topoplot(pl.data_evo, RESULTS.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',viridis,...
    'whitebk','on');
title(sprintf('evoked SSVEP amps for stim pos %s\nN = %1.0f | [%1.0f %1.0f]ms', ...
    pl.pos2plot, numel(pl.sub2plot), min([RESULTS.ffttimewin{pl.time2plot}])*1000, max([RESULTS.ffttimewin{pl.time2plot}])*1000))
colorbar

%% extract data for R
% large center as in tango | periphery: central and lateral 
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({RESULTS.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.sub2sel = 1:numel(F.Subs2use);
pl.baseline = 1;

pl.RDKlabel = {'RDK1';'RDK2';'RDK3';'RDK4';'RDK5'};
pl.RDKposlabel2 = {'center';'center';'peri';'peri';'peri'};
pl.timelabel = cellfun(@(x) vararg2str(x),RESULTS.ffttimewin(:),'UniformOutput',false);

% preallocate data
dataout.SSVEPamps = nan(5,numel(F.conlabel_att),numel(RESULTS.ffttimewin),numel(pl.sub2sel));

% get all descriptives
% the ones that are fixed
pl.RDK.con = permute(repmat((1:numel(F.conlabel_att))',[1 size(dataout.SSVEPamps,[1 3 4])]),[2 1 3 4]);
pl.RDK.RDK_id = repmat(pl.RDKlabel,[1 size(dataout.SSVEPamps, [2 3 4 ])]);
pl.RDK.RDK_ispresent = permute(repmat(F.conRDKpresent_label,[1 1 size(dataout.SSVEPamps,[3 4])]),[2 1 3 4]);
pl.RDK.RDK_isattended = permute(repmat(F.conRDKcolattended_label,[1 1 size(dataout.SSVEPamps,[3 4])]),[2 1 3 4]);
pl.RDK.RDK_isattended2 = permute(repmat(F.conRDKcolattended_label2,[1 1 size(dataout.SSVEPamps,[3 4])]),[2 1 3 4]);
pl.RDK.timewin = permute(repmat(pl.timelabel,[1 size(dataout.SSVEPamps,[1 2 4])]), [2 3 1 4]);
pl.RDK.sub = permute(repmat(F.Subs2use(pl.sub2sel)',[1 size(dataout.SSVEPamps,[1 2 3])]),[2 3 4 1]);
pl.RDK.RDK_pos2 = repmat(pl.RDKposlabel2,[1 size(dataout.SSVEPamps, [2 3 4 ])]);

% the onses that change
pl.RDK.RDK_freq = dataout.SSVEPamps;
pl.RDK.RDK_color = repmat({''},size(dataout.SSVEPamps));
pl.RDK.RDK_pos1 = repmat({''},size(dataout.SSVEPamps));
pl.RDK.RDK_electrodes = repmat({''},size(dataout.SSVEPamps));

% loop across participants
for i_sub = 1:numel(pl.sub2sel)
    for i_rdk = 1:5
        % extract data for respective RDK at respective electrode cluster
        % which position
        t.posidx = strcmp(pl.elec2plot(:,2),RESULTS.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel);

        % extract data by averaging across electrodes
        dataout.SSVEPamps(i_rdk,:,:,i_sub) = mean(RESULTS.SSVEPamps(i_rdk,pl.elec2plot_i{t.posidx},:,:,pl.sub2sel(i_sub)),2);

        % write some data
        pl.RDK.RDK_freq(i_rdk,:,:,i_sub) = RESULTS.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).freq;
        pl.RDK.RDK_color(i_rdk,:,:,i_sub) = {RESULTS.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).colnames};
        pl.RDK.RDK_pos1(i_rdk,:,:,i_sub) = {RESULTS.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel};
        pl.RDK.RDK_electrodes(i_rdk,:,:,i_sub) = {vararg2str(pl.elec2plot(t.posidx,1))};
    end
end

% baseline corrected data in percent
dataout.SSVEPmods = 100*(bsxfun(@rdivide, dataout.SSVEPamps, dataout.SSVEPamps(:,:,pl.baseline,:))-1);

R_Mat.all = [{'amplitude_evoked','modulation_evoked', ...
    'subjects', 'condition', 'time', ...
    'RDK_id', 'RDK_position1','RDK_position2', 'RDK_freq', 'RDK_color', 'RDK_ispresented', 'RDK_isattended', 'RDK_isattended2', 'RDK_electrodes'}; ...
    num2cell([dataout.SSVEPamps(:) dataout.SSVEPmods(:) pl.RDK.sub(:) pl.RDK.con(:)]) ...
    pl.RDK.timewin(:) pl.RDK.RDK_id(:) pl.RDK.RDK_pos1(:) pl.RDK.RDK_pos2(:) num2cell(pl.RDK.RDK_freq(:)) pl.RDK.RDK_color(:) ...
    pl.RDK.RDK_ispresent(:) pl.RDK.RDK_isattended(:) pl.RDK.RDK_isattended2(:) pl.RDK.RDK_electrodes(:)
    ];

R_Mat.all_table = cell2table(R_Mat.all(2:end,:), "VariableNames",R_Mat.all(1,:));
t.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclust_%s.csv',t.datestr)),'Delimiter',';')

