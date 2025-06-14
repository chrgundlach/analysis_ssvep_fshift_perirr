%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\EEG\TFA'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% F.Subs2use              = [1:13 15:21];
% changed experiment from participant 22 onwards (stimuli isoluminant to
% background and used other frequencies
% participant 42 has lower trial number
F.Subs2use              = [1:13 15:41 43:52];
                        
F.TFA.baseline          = [-500 -250];

F.SSVEP_Freqs{1}        = [17 20 23 26 29]; % sub 1 to 21
F.SSVEP_Freqs{2}        = [14 17 20 23 26]; % sub 22 onwars
F.RDK_pos               = [0 -255 255];
F.RDK_pos_label         = {'center';'left';'right'};

F.SSVEP_Freqsforsubs    = ones(size(F.Subs2use)); F.SSVEP_Freqsforsubs(F.Subs2use>21)=2;


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
%[F.conRDKcolattended_label2{strcmp(F.conRDKpresent_label,'not presented')}] = deal('not presented');


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_tfa.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.tfa = open(fullfile(F.PathInEEG,sprintf('VP%s_tfa.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % convert to single
    temp.tfa.TFA.data_evo = single(temp.tfa.TFA.data_evo);
    temp.tfa.TFA.data_ind = single(temp.tfa.TFA.data_ind);
    temp.tfa.TFA.FFT.data_evo = single(temp.tfa.TFA.FFT.data_evo);
    temp.tfa.TFA.FFT.data_ind = single(temp.tfa.TFA.FFT.data_ind);
    
    
    % preallocate memory
    if i_sub == 1
        TFA.data_evo = single(nan([size(temp.tfa.TFA.data_evo),numel(F.Subs2use)]));
        TFA.data_ind = single(nan([size(temp.tfa.TFA.data_ind),numel(F.Subs2use)]));
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.con_trialnum = temp.tfa.TFA.con_trialnum;
        TFA.srate = temp.tfa.TFA.params.srate/2;
        TFA.fftdata_ind = nan([size(temp.tfa.TFA.FFT.data_ind),numel(F.Subs2use)]);
        TFA.fftdata_evo = nan([size(temp.tfa.TFA.FFT.data_evo),numel(F.Subs2use)]);
        TFA.ffttimewin = temp.tfa.TFA.FFT.timewin;
        TFA.fftfreqs = temp.tfa.TFA.FFT.freqs;
    end
    
    % assign data
    TFA.data_evo(:,:,:,:,i_sub) = temp.tfa.TFA.data_evo; % evoked data
    TFA.data_ind(:,:,:,:,i_sub) = temp.tfa.TFA.data_ind; % induced data
    %     TFA(i_exp).data_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
    %         mean(temp.tfa.TFA.data(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:,:),2));
%     TFA(i_exp).data_bc(:,:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data, ...
%         mean(temp.tfa.TFA.data(:,eeg_time2points(F.TFA.baseline(1),TFA(i_exp).time):eeg_time2points(F.TFA.baseline(2),TFA(i_exp).time),:,:,:),2)))-1);
    TFA.fftdata_evo(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_evo;
    TFA.fftdata_ind(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_ind;
    TFA.RDK(i_sub) = temp.tfa.TFA.RDK;

    % add position of RDK to data
    t.label = cellfun(@(x) F.RDK_pos_label(x(1)==F.RDK_pos),{TFA.RDK(i_sub).RDK.centershift});
    [TFA.RDK(i_sub).RDK.poslabel] = deal(t.label{:});
    
    clear temp
    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];



%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'};
% pl.elec2plot = {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'};
% pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'};
% pl.elec2plot = {'POz';'O1';'Oz';'I2';'Iz'};
% pl.elec2plot = {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'};
% pl.elec2plot = {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});

%% plot grand mean FFT data | spectra | for distinct frequencies (lookup of respective electrode cluster)
% plotting parameters
pl.elec2plot = {{'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
    {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

pl.elec2plot = {{'POz';'O1';'Oz';'I2';'Iz'}, 'center';...
    {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
    {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% small adapted
pl.elec2plot = {{'Oz';'Iz'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% small adapted 2
pl.elec2plot = {{'Oz';'Iz';'O1';'O2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% % very large center
pl.elec2plot = {{'P3';'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P4';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% large center as in tango | periphery: central and lateral 
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO4';'POz';'Oz';'O2'}, 'right'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


pl.time2plot = [1];
% pl.time2plot = [1 2 3];
pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot = find(F.Subs2use<22); % luminance offset
pl.freq2plot=F.SSVEP_Freqs{1}(1);

pl.sub2plot = find(F.Subs2use>21); % isoluminant to background
pl.freq2plot=F.SSVEP_Freqs{2}(4);

% extract data
pl.data_ind = nan(size(TFA.fftdata_ind,1), numel(pl.sub2plot)); pl.data_evo = pl.data_ind;
t.elidx = false(size(pl.elec2plot(:,2),1),numel(pl.sub2plot));
for i_sub = 1:numel(pl.sub2plot)
    % index position of frequency
    t.rdkidx = [TFA.RDK(pl.sub2plot(i_sub)).RDK.freq]==pl.freq2plot;
    t.posidx = TFA.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).poslabel;
    t.elidx(:,i_sub) = strcmp(pl.elec2plot(:,2), t.posidx);
    % index conditions for which the RDK is shown
    t.conidx = any(t.rdkidx & F.conRDKpresent,2);

    % extract data
    pl.data_ind(:,i_sub) = mean(TFA.fftdata_ind(:,pl.elec2plot_i{t.elidx(:,i_sub)},t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[2,3,4]);
    pl.data_evo(:,i_sub) = mean(TFA.fftdata_evo(:,pl.elec2plot_i{t.elidx(:,i_sub)},t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[2,3,4]);
end


% plotting
figure;
set(gcf,'Position',[100 100 600 600],'PaperPositionMode','auto')
subplot(2,1,1)
plot(TFA.fftfreqs,pl.data_ind,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,mean(pl.data_ind,2),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('induced GrandMean FFT spectra | N = %1.0f | FOI = %1.1f Hz', ...
    numel(pl.sub2plot), pl.freq2plot),'Interpreter','none')
vline(pl.freq2plot,'k:')
% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i(any(t.elidx,2))),1)),'o','r',4,1});

subplot(2,1,2)
hold on;
plot(TFA.fftfreqs,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(TFA.fftfreqs,mean(pl.data_evo,2),'Color','k','LineWidth',2)
xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f | FOI = %1.1f Hz', ...
    numel(pl.sub2plot), pl.freq2plot),'Interpreter','none')
vline(pl.freq2plot,'k:')
box on

% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.28 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i(any(t.elidx,2))),1)),'o','r',4,1});


any(t.elidx,2)

%% plot Grand Mean FFT data | topoplot for positions (as frequencies are random)
pl.time2plot = [1:3];
pl.time2plot = [1];
% pl.pos2plot='center';
pl.pos2plot='right';
% pl.pos2plot='left';
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);
% pl.sub2sel = find(F.Subs2use<22); % luminance offset
pl.sub2sel = find(F.Subs2use>21); % isoluminant to background


pl.sub2plot = pl.sub2sel( ...
    cellfun(@(x) any(strcmp({x.poslabel},pl.pos2plot)), ... % index for which subject RDKs were presented at pl.pos2plot
    {TFA.RDK(pl.sub2sel).RDK}));

% extract data differently [is it converging?]
% first, extract data for all of them
pl.data_ind = []; pl.data_evo = [];
t.data_ind = nan(numel(TFA.electrodes),numel(TFA.RDK(1).RDK),numel(F.conlabel_att),numel(pl.sub2sel));
t.data_evo = t.data_ind;
t.data_freq = t.data_ind; 
t.dat_position = repmat({''},size(t.data_ind));
t.data_ispresent = permute(repmat(F.conRDKpresent',[1 1 numel(TFA.electrodes) numel(pl.sub2sel)]), [3 1 2 4]); 
for i_sub = 1:numel(pl.sub2sel)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        t.freq = TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).freq; % which SSVEP frequency?
        t.idx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
        % extract data
        t.data_ind(:,i_rdk,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.idx(1):t.idx(2),:,:,pl.time2plot,pl.sub2sel(i_sub)),[1 4]));
        t.data_evo(:,i_rdk,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.idx(1):t.idx(2),:,:,pl.time2plot,pl.sub2sel(i_sub)),[1 4]));
        
        % write some data
        t.data_freq(:,i_rdk,:,i_sub) = t.freq;
        [t.dat_position{:,i_rdk,:,i_sub}] = deal(TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel);
    end
end

% now extract the data
pl.idx = strcmp(t.dat_position,pl.pos2plot);
pl.data_ind = t.data_ind; pl.data_evo = t.data_evo;
pl.data_ind(~pl.idx)=nan;
pl.data_evo(~pl.idx)=nan;
pl.data_ind = squeeze(mean(pl.data_ind,[2,3],'omitnan'));
pl.data_evo = squeeze(mean(pl.data_evo,[2,3],'omitnan'));


figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')

% induced
h.s(1) = subplot(1,2,1);
pl.mdata = mean(pl.data_ind,2,'omitnan');
pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, TFA.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('induced SSVEP amps at %s\nN = %1.0f | [%1.0f %1.0f]ms', ...
    pl.pos2plot, sum(~isnan(pl.data_ind(1,:))), min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
colorbar

h.s(2) = subplot(1,2,2);
pl.mdata = mean(pl.data_evo,2,'omitnan');
pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, TFA.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('evoked SSVEP amps at %s\nN = %1.0f | [%1.0f %1.0f]ms', ...
    pl.pos2plot, sum(~isnan(pl.data_ind(1,:))), min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
colorbar


%% plot FFT data modulation | topoplot effects on central stimuli
pl.time_base = 1;
pl.time2plot = [2];
pl.pos2plot='center';
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);
% pl.sub2sel = find(F.Subs2use<22); % luminance offset
pl.sub2sel = find(F.Subs2use>21); % isoluminant to background
pl.con2plot = F.conRDKcolattended_label(1:2);


pl.sub2plot = pl.sub2sel( ...
    cellfun(@(x) any(strcmp({x.poslabel},pl.pos2plot)), ... % index for which subject RDKs were presented at pl.pos2plot
    {TFA.RDK(pl.sub2sel).RDK}));

% extract data differently [is it converging?]
% first, extract data for all of them
pl.data_ind = []; pl.data_evo = [];
t.data_ind = nan(numel(TFA.electrodes),numel(TFA.RDK(1).RDK),numel(F.conlabel_att),numel(TFA.ffttimewin),numel(pl.sub2sel));
t.data_evo = t.data_ind;
t.data_freq = t.data_ind; 
t.dat_position = repmat({''},size(t.data_ind));
t.data_ispresent = permute(repmat(F.conRDKpresent',[1 1 numel(TFA.electrodes) numel(pl.sub2sel)]), [3 1 2 4]); 
for i_sub = 1:numel(pl.sub2sel)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        t.freq = TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).freq; % which SSVEP frequency?
        t.idx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
        % extract data
        t.data_ind(:,i_rdk,:,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.idx(1):t.idx(2),:,:,:,pl.sub2sel(i_sub)),[1]));
        t.data_evo(:,i_rdk,:,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.idx(1):t.idx(2),:,:,:,pl.sub2sel(i_sub)),[1]));
        
        % write some data
        t.data_freq(:,i_rdk,:,i_sub) = t.freq;
        [t.dat_position{:,i_rdk,:,:,i_sub}] = deal(TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel);
        
    end
end

% now extract the data and collapse across RDKs
t.rdkidx = find(sum(strcmp(t.dat_position,pl.pos2plot),[1,3 4 5 6]));
t.data_ind_coll = nan(numel(TFA.electrodes),numel(t.rdkidx),numel(pl.con2plot),numel(TFA.ffttimewin),numel(pl.sub2sel));
t.data_evo_coll = t.data_ind_coll;
for i_rdk = 1:numel(t.rdkidx)
    for i_con = 1:numel(pl.con2plot)
        t.conidx = strcmp(F.conRDKcolattended_label(:,t.rdkidx(i_rdk)),pl.con2plot{i_con});
        t.data_ind_coll(:,i_rdk,i_con,:,:) = mean(t.data_ind(:,t.rdkidx(i_rdk),t.conidx,:,:),3);
        t.data_evo_coll(:,i_rdk,i_con,:,:) = mean(t.data_evo(:,t.rdkidx(i_rdk),t.conidx,:,:),3);
    end
end

t.data_ind_coll_bc = ((mean(t.data_ind_coll(:,:,:,pl.time2plot,:),4) ./ t.data_ind_coll(:,:,:,pl.time_base,:))-1)*100;
t.data_evo_coll_bc = ((mean(t.data_evo_coll(:,:,:,pl.time2plot,:),4) ./ t.data_evo_coll(:,:,:,pl.time_base,:))-1)*100;
% t.data_ind_coll_bc = mean(t.data_ind_coll(:,:,:,pl.time2plot,:),4) - t.data_ind_coll(:,:,:,pl.time_base,:);
% t.data_evo_coll_bc = mean(t.data_evo_coll(:,:,:,pl.time2plot,:),4) - t.data_evo_coll(:,:,:,pl.time_base,:);

% collapse across RDKs
t.data_ind_coll_bc = squeeze(mean(t.data_ind_coll_bc,2));
t.data_evo_coll_bc = squeeze(mean(t.data_evo_coll_bc,2));

% difference between conditions
t.data_ind_coll_bc(:,3,:) = diff(fliplr(t.data_ind_coll_bc),1,2);
t.data_evo_coll_bc(:,3,:) = diff(fliplr(t.data_evo_coll_bc),1,2);

t.data_ind_coll_bc_m = mean(t.data_ind_coll_bc,3);
t.data_evo_coll_bc_m = mean(t.data_evo_coll_bc,3);

% do ttests
[t.data_ind_coll_ttp t.data_evo_coll_ttp] = deal(nan(size(t.data_ind_coll_bc_m)));

for i_con = 1:size(t.data_ind_coll_ttp,2)
    [tt.h t.data_ind_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_ind_coll_bc(:,i_con,:))');
    [tt.h t.data_evo_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_evo_coll_bc(:,i_con,:))');
end

% evoked
figure;
set(gcf,'Position',[100 100 800 500],'PaperPositionMode','auto')
pl.conlabel = [F.conRDKcolattended_label(1:2) 'diff'];

% modulations
s

% uncorrected t-values
for i_con = 1:size(t.data_evo_coll_ttp,2)
    t.data = abs(log10(t.data_evo_coll_ttp(:,i_con)));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end

    h.s2(i_con)=subplot(2,size(t.data_evo_coll_ttp,2),i_con+size(t.data_ind_coll_ttp,2));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap,'whitebk','on');
    title(sprintf('respective t-values'))
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))

end


% induced
figure;
set(gcf,'Position',[100 100 800 500],'PaperPositionMode','auto')
pl.conlabel = [F.conRDKcolattended_label(1:2) 'diff'];

% modulations
for i_con = 1:size(t.data_ind_coll_bc_m,2)
    h.s(i_con) = subplot(2,size(t.data_ind_coll_bc_m,2),i_con);
    pl.clim = [-1 1] *max(abs(t.data_ind_coll_bc_m),[],'all');
    topoplot(t.data_ind_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',flipud(cbrewer2('RdBu')),...
        'whitebk','on');
    title(sprintf('induced SSVEP mod | %s\n[%1.0f %1.0f]ms', ...
        pl.conlabel{i_con}, TFA.ffttimewin{pl.time2plot}*1000))
    colorbar
end

% uncorrected t-values
for i_con = 1:size(t.data_evo_coll_ttp,2)
    t.data = abs(log10(t.data_ind_coll_ttp(:,i_con)));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end

    h.s2(i_con)=subplot(2,size(t.data_ind_coll_ttp,2),i_con+size(t.data_ind_coll_ttp,2));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap, 'whitebk','on');
    title(sprintf('respective t-values'))
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))

end

%% plot FFT data modulation | topoplot effects for peripheral stimuli [right hemisphere is contralateral to stimulus]
pl.time_base = 1;
pl.time2plot = [2];
pl.pos2plot= {'left'; 'right'};
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);
% pl.sub2sel = find(F.Subs2use<22); % luminance offset
pl.sub2sel = find(F.Subs2use>21); % isoluminant to background
pl.con2plot = unique(F.conRDKcolattended_label(:,3:end));


% extract data differently [is it converging?]
% first, extract data for all of them
pl.data_ind = []; pl.data_evo = [];
t.data_ind = nan(numel(TFA.electrodes),numel(TFA.RDK(1).RDK),numel(F.conlabel_att),numel(TFA.ffttimewin),numel(pl.sub2sel));
t.data_evo = t.data_ind;
t.data_freq = t.data_ind; 
t.dat_position = repmat({''},size(t.data_ind));
t.data_ispresent = permute(repmat(F.conRDKpresent',[1 1 numel(TFA.electrodes) numel(pl.sub2sel)]), [3 1 2 4]); 
for i_sub = 1:numel(pl.sub2sel)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        t.freq = TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).freq; % which SSVEP frequency?
        t.idx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
        
        % check for position of peripheral stimulus
        if strcmp(TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel,'right')
            t.elec_idx = eeg_elec_hemisphere_collapse({TFA.electrodes.labels},2);
        else
            t.elec_idx = 1:64;
        end

        % extract data
        t.data_ind(:,i_rdk,:,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.idx(1):t.idx(2),t.elec_idx,:,:,pl.sub2sel(i_sub)),[1]));
        t.data_evo(:,i_rdk,:,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.idx(1):t.idx(2),t.elec_idx,:,:,pl.sub2sel(i_sub)),[1]));
        
        % write some data
        t.data_freq(:,i_rdk,:,i_sub) = t.freq;
        [t.dat_position{:,i_rdk,:,:,i_sub}] = deal(TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).poslabel);
        
    end
end

% now extract the data and collapse across RDKs
t.rdkidx = [3 4 5];
t.data_ind_coll = nan(numel(TFA.electrodes),numel(t.rdkidx),numel(pl.con2plot),numel(TFA.ffttimewin),numel(pl.sub2sel));
t.data_evo_coll = t.data_ind_coll;
for i_rdk = 1:numel(t.rdkidx)
    for i_con = 1:numel(pl.con2plot)
        t.conidx = strcmp(F.conRDKcolattended_label(:,t.rdkidx(i_rdk)),pl.con2plot{i_con}) & F.conRDKpresent(:,t.rdkidx(i_rdk));
        t.data_ind_coll(:,i_rdk,i_con,:,:) = mean(t.data_ind(:,t.rdkidx(i_rdk),t.conidx,:,:),3);
        t.data_evo_coll(:,i_rdk,i_con,:,:) = mean(t.data_evo(:,t.rdkidx(i_rdk),t.conidx,:,:),3);
    end
end

% t.data_ind_coll_bc = ((mean(t.data_ind_coll(:,:,:,pl.time2plot,:),4) ./ t.data_ind_coll(:,:,:,pl.time_base,:))-1)*100;
% t.data_evo_coll_bc = ((mean(t.data_evo_coll(:,:,:,pl.time2plot,:),4) ./ t.data_evo_coll(:,:,:,pl.time_base,:))-1)*100;
t.data_ind_coll_bc = mean(t.data_ind_coll(:,:,:,pl.time2plot,:),4) - t.data_ind_coll(:,:,:,pl.time_base,:);
t.data_evo_coll_bc = mean(t.data_evo_coll(:,:,:,pl.time2plot,:),4) - t.data_evo_coll(:,:,:,pl.time_base,:);

% collapse across RDKs
t.data_ind_coll_bc = squeeze(mean(t.data_ind_coll_bc,2,"omitnan"));
t.data_evo_coll_bc = squeeze(mean(t.data_evo_coll_bc,2,"omitnan"));

% difference between conditions
t.contrasts = [1 2; 1 3; 2 3];
for i_cont = 1:size(t.contrasts,1)
    t.data_ind_coll_bc(:,end+1,:) = diff(t.data_ind_coll_bc(:,fliplr(t.contrasts(i_cont,:)),:),1,2);
    t.data_evo_coll_bc(:,end+1,:) = diff(t.data_evo_coll_bc(:,fliplr(t.contrasts(i_cont,:)),:),1,2);
end

t.data_ind_coll_bc_m = mean(t.data_ind_coll_bc,3);
t.data_evo_coll_bc_m = mean(t.data_evo_coll_bc,3);

% do ttests
[t.data_ind_coll_ttp t.data_evo_coll_ttp] = deal(nan(size(t.data_ind_coll_bc_m)));

for i_con = 1:size(t.data_ind_coll_ttp,2)
    [tt.h t.data_ind_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_ind_coll_bc(:,i_con,:))');
    [tt.h t.data_evo_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_evo_coll_bc(:,i_con,:))');
end

% evoked
figure;
set(gcf,'Position',[100 100 1400 400],'PaperPositionMode','auto')
pl.conlabel = pl.con2plot;
for i_cont = 1:size(t.contrasts,1)
    pl.conlabel{end+1} = sprintf('diff %s - %s',pl.con2plot{t.contrasts(i_cont,:)});
end

% modulations
for i_con = 1:size(t.data_evo_coll_bc_m,2)
    h.s(i_con) = subplot(2,size(t.data_evo_coll_bc,2),i_con);
    pl.clim = [-1 1] *max(abs(t.data_evo_coll_bc_m),[],'all');
%     topoplot(t.data_evo_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',fake_parula,...
%         'whitebk','on');
    topoplot(t.data_evo_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',flipud(cbrewer2('RdBu')),...
        'whitebk','on');
   
%     title(sprintf('evoked SSVEP mod | [%1.0f %1.0f]ms\n%s', ...
%         TFA.ffttimewin{pl.time2plot}*1000, pl.conlabel{i_con}), ...
%         'FontSize',6)
    title(sprintf('evoked SSVEP diff | [%1.0f %1.0f]ms\n%s', ...
        TFA.ffttimewin{pl.time2plot}*1000, pl.conlabel{i_con}), ...
        'FontSize',6)
    colorbar
end

% uncorrected t-values
for i_con = 1:size(t.data_evo_coll_ttp,2)
    t.data = abs(log10(t.data_evo_coll_ttp(:,i_con)));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end

    h.s2(i_con)=subplot(2,size(t.data_evo_coll_ttp,2),i_con+size(t.data_ind_coll_ttp,2));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap, 'whitebk','on');
    title(sprintf('respective t-values'), 'FontSize',6)
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))

end


% induced
figure;
set(gcf,'Position',[100 100 1400 400],'PaperPositionMode','auto')
pl.conlabel = pl.con2plot;
for i_cont = 1:size(t.contrasts,1)
    pl.conlabel{end+1} = sprintf('diff %s - %s',pl.con2plot{t.contrasts(i_cont,:)});
end

% modulations
for i_con = 1:size(t.data_evo_coll_bc_m,2)
    h.s(i_con) = subplot(2,size(t.data_ind_coll_bc,2),i_con);
    pl.clim = [-1 1] *max(abs(t.data_ind_coll_bc_m),[],'all');
%     topoplot(t.data_evo_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',fake_parula,...
%         'whitebk','on');
    topoplot(t.data_ind_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',flipud(cbrewer2('RdBu')),...
        'whitebk','on');    
   
    title(sprintf('induced SSVEP mod | [%1.0f %1.0f]ms\n%s', ...
        TFA.ffttimewin{pl.time2plot}*1000, pl.conlabel{i_con}), ...
        'FontSize',6)
    
    colorbar
end

% uncorrected t-values
for i_con = 1:size(t.data_evo_coll_ttp,2)
    t.data = abs(log10(t.data_ind_coll_ttp(:,i_con)));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end

    h.s2(i_con)=subplot(2,size(t.data_ind_coll_ttp,2),i_con+size(t.data_ind_coll_ttp,2));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap, 'whitebk','on');
    title(sprintf('respective t-values'), 'FontSize',6)
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))

end


%% extract amplitude values for FFT
% plotting parameters
pl.elec2plot = {{'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
    {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% smaller for center
% pl.elec2plot = {{'POz';'O1';'Oz';'I2';'Iz'}, 'center';...
%     {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
%     {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% small adapted
pl.elec2plot = {{'Oz';'Iz'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% small adapted II
pl.elec2plot = {{'Oz';'Iz';'O1';'O2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% very large center
% pl.elec2plot = {{'P3';'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P4';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
%     {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
%     {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'}, 'right'};

% large center as in tango | periphery: central and lateral smaller (fits topographies)
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO4';'POz';'Oz';'O2'}, 'right'};

% % large center as in tango | periphery: only central
% pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
%     {'PO3';'POz';'Oz';'O1'}, 'left'; ...
%     {'PO2';'POz';'Oz';'O2'}, 'right'};
% 
% % % large center as in tango | periphery: only lateral 
% pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
%     {'P8';'P10';'PO8'}, 'left'; ...
%     {'P7';'P9';'PO7'}, 'right'};
% 
% 
% % large center as in tango | as in tango periphery: central and lateral 
% pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center';...
%     {'P8';'P10';'PO8';'PO4';'O2';'I2';'P3';'P1';'Pz';'P2';'PO3';'POz'}, 'left'; ...
%     {'P7';'P9';'PO7';'PO3';'O1';'I1';'P4';'P2';'Pz';'P1';'PO4';'POz'}, 'right'};




pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


pl.time2plot = [1:3];
pl.time2baseline = [1];
pl.freqrange=[-0.1 0.1];
% pl.freqrange=[0 0];
pl.sub2plot = 1:numel(F.Subs2use);

pl.colorlum = repmat({'offset_to_bckgrd'}, 1, numel(pl.sub2plot));
pl.colorlum(F.SSVEP_Freqsforsubs(pl.sub2plot) == 2) = {'isolum__to_bckgrd'};


pl.data_ind = []; pl.data_evo = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3';'RDK4';'RDK5'};
pl.RDKposlabel2 = {'center';'center';'peri';'peri';'peri'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);

% extract data
% ([RDK num], con, time, sub)
pl.RDK.data_ind = nan(numel(TFA.RDK(1).RDK),numel(F.conlabel_att),numel(pl.time2plot),numel(pl.sub2plot));
pl.RDK.data_evo = pl.RDK.data_ind;


% the ones that are fixed
pl.RDK.con = permute(repmat((1:numel(F.conlabel_att))',[1 size(pl.RDK.data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.RDK_id = repmat(pl.RDKlabel,[1 size(pl.RDK.data_ind, [2 3 4 ])]);
pl.RDK.RDK_ispresent = permute(repmat(F.conRDKpresent_label,[1 1 size(pl.RDK.data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.RDK_isattended = permute(repmat(F.conRDKcolattended_label,[1 1 size(pl.RDK.data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.RDK_isattended2 = permute(repmat(F.conRDKcolattended_label2,[1 1 size(pl.RDK.data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.timewin = permute(repmat(pl.timelabel,[1 size(pl.RDK.data_ind,[1 2 4])]), [2 3 1 4]);
pl.RDK.sub = permute(repmat(F.Subs2use(pl.sub2plot)',[1 size(pl.RDK.data_ind,[1 2 3])]),[2 3 4 1]);
pl.RDK.RDK_pos2 = repmat(pl.RDKposlabel2,[1 size(pl.RDK.data_ind, [2 3 4 ])]);

% the onses that change
pl.RDK.RDK_freq = pl.RDK.data_ind;
pl.RDK.RDK_color = repmat({''},size(pl.RDK.data_ind));
pl.RDK.RDK_colorlum = repmat({''},size(pl.RDK.data_ind));
pl.RDK.RDK_pos1 = repmat({''},size(pl.RDK.data_ind));
pl.RDK.RDK_electrodes = repmat({''},size(pl.RDK.data_ind));

for i_sub = 1:numel(pl.sub2plot)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        t.freq = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).freq; % which SSVEP frequency?
        t.fidx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
        t.posidx = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).poslabel;
        t.elidx = strcmp(pl.elec2plot(:,2), t.posidx);
        
        % extract data
        pl.RDK.data_ind(i_rdk,:,:,i_sub) = squeeze(mean( ...
            TFA.fftdata_ind(t.fidx(1):t.fidx(2),pl.elec2plot_i{t.elidx},:,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
            ));
        pl.RDK.data_evo(i_rdk,:,:,i_sub) = squeeze(mean( ...
            TFA.fftdata_evo(t.fidx(1):t.fidx(2),pl.elec2plot_i{t.elidx},:,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %evoked
            ));

        % write some data
        pl.RDK.RDK_freq(i_rdk,:,:,i_sub) = t.freq;
        pl.RDK.RDK_color(i_rdk,:,:,i_sub) = {TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).colnames};
        pl.RDK.RDK_colorlum(i_rdk,:,:,i_sub) = pl.colorlum(pl.sub2plot(i_sub));
        pl.RDK.RDK_pos1(i_rdk,:,:,i_sub) = {TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).poslabel};
        pl.RDK.RDK_electrodes(i_rdk,:,:,i_sub) = {vararg2str(pl.elec2plot(t.elidx,1))};
    end
end

% baseline corrected data
pl.RDK.data_ind_bc = 100*(bsxfun(@rdivide, pl.RDK.data_ind, pl.RDK.data_ind(:,:,1,:))-1);
pl.RDK.data_evo_bc = 100*(bsxfun(@rdivide, pl.RDK.data_evo, pl.RDK.data_evo(:,:,1,:))-1);
pl.RDK.data_ind_bc_sub = bsxfun(@minus,  pl.RDK.data_ind, pl.RDK.data_ind(:,:,1,:));
pl.RDK.data_evo_bc_sub = bsxfun(@minus, pl.RDK.data_evo, pl.RDK.data_evo(:,:,1,:));

% % trouble shooting do subtraction and modulation correspond?
% % pl.tdata = reshape(squeeze(pl.RDK.data_evo_bc(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4));
% % pl.tdata(:,:,2) = reshape(squeeze(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4)).*10;
% pl.tdata = reshape(squeeze(mean(pl.RDK.data_evo_bc(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(mean(pl.RDK.data_evo_bc_sub(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4)).*10;
% figure;
% for i_sub = 1:size(pl.tdata,2)
% %     plot([-0.25 +0.25]+i_sub, squeeze(pl.tdata(:,i_sub,:)),'Color',[0.3 0.3 0.3 0.5])
%     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% 
% pl.tdata = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% figure;
% for i_sub = 1:size(pl.tdata,2)
%     plot([-0.25 +0.25]+i_sub, ...
%         squeeze(pl.tdata(:,i_sub,:))+repmat(((randn(size(pl.tdata,1),1)-0.5).*0.1),1,2), ...
%         'Color',[0.3 0.3 0.3 0.5])
% %     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% figure; plot(sign(pl.RDK.data_evo_bc_sub(:)), sign(pl.RDK.data_evo_bc(:)))
% figure; plot(sign(R_Mat.all_table.modulation_evoked), sign(R_Mat.all_table.subtraction_evoked))


% % do some interim plotting for checking everything
% t.tdata = cat(3, ...
%     [squeeze(mean(pl.RDK.data_evo(1,[1 3 5],2,:),2)) squeeze(mean(pl.RDK.data_evo(1,[2 4 6],2,:),2))], ...
%     [squeeze(mean(pl.RDK.data_evo(2,[2 4 6],2,:),2)) squeeze(mean(pl.RDK.data_evo(2,[1 3 5],2,:),2))]);
% figure; boxplot(mean(t.tdata,3))
% t.tidx = strcmp(pl.RDK.RDK_pos1,'center') & strcmp(pl.RDK.timewin,'[0.5 1.5] ') & pl.RDK.RDK_freq == 29 ;
% t.tdata_ev0 =  pl.RDK.data_evo_bc;  t.tdata_ev0(~t.tidx)= nan;
% t.tdata = cat(3, ...
%     [squeeze(mean(t.tdata_ev0(1,[1 3 5],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(1,[2 4 6],2,:),2,'omitnan'))], ...
%     [squeeze(mean(t.tdata_ev0(2,[2 4 6],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(2,[1 3 5],2,:),2,'omitnan'))]);
% figure; boxplot(mean(t.tdata,3,'omitnan'))
% figure; plot(mean(t.tdata,3,'omitnan')')
% t.tdata2 =diff(mean(t.tdata,3,'omitnan'),1,2);
% figure; boxplot(diff(mean(t.tdata,3,'omitnan'),1,2));
% figure; histfit(diff(mean(t.tdata,3,'omitnan'),1,2),20)
% [tt.h,tt.p,tt.ci,tt.stats] = ttest(diff(mean(t.tdata,3,'omitnan'),1,2));

R_Mat.all = [{'amplitude_induced','amplitude_evoked','modulation_induced','modulation_evoked','subtraction_induced','subtraction_evoked', ...
    'subjects', 'condition', 'time', ...
    'RDK_id', 'RDK_position1','RDK_position2', 'RDK_freq', 'RDK_color','RDK_colorlum', 'RDK_ispresented', 'RDK_isattended', 'RDK_isattended2', 'RDK_electrodes'}; ...
    num2cell([pl.RDK.data_ind(:) pl.RDK.data_evo(:) pl.RDK.data_ind_bc(:) pl.RDK.data_evo_bc(:) pl.RDK.data_ind_bc_sub(:) pl.RDK.data_evo_bc_sub(:) ...
    pl.RDK.sub(:) pl.RDK.con(:)]) ...
    pl.RDK.timewin(:) pl.RDK.RDK_id(:) pl.RDK.RDK_pos1(:) pl.RDK.RDK_pos2(:) num2cell(pl.RDK.RDK_freq(:)) pl.RDK.RDK_color(:) ...
    pl.RDK.RDK_colorlum(:) pl.RDK.RDK_ispresent(:) pl.RDK.RDK_isattended(:) pl.RDK.RDK_isattended2(:) pl.RDK.RDK_electrodes(:)
    ];

R_Mat.all_table = cell2table(R_Mat.all(2:end,:), "VariableNames",R_Mat.all(1,:));


t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftperirr\data_in';
% t.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
% xlswrite(fullfile(t.path,sprintf('FFT_Amp_data_largeclust_%s.csv',t.datestr)),R_Mat.all)
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclust_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclust_allsubs_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_smallclust_allsubs_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largecenterclust_allsubs_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_TangoStudy_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_TangoStudyPeriOnlyCenter_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_TangoStudyPeriOnlyLateral_%s.csv',t.datestr)),'Delimiter',';')

% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_CenterTangoStudyPeriSmall_%s.csv',t.datestr)),'Delimiter',';')


%% actual plotting data | TFA Grand Mean timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.flims= TFA.frequency([1 end]);
pl.flims_i=dsearchn(TFA.frequency', pl.flims');

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
pl.base_i = dsearchn(TFA.time', pl.base');

% pl.sub2sel = find(F.Subs2use<22); % luminance offset
pl.sub2sel = find(F.Subs2use>21); % isoluminant to background

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
% pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


pl.data_ind = squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4));
pl.data_evo = squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4));
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, ...
    mean(squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4)),2))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, ...
    mean(squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4)),2))-1);

figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | induced\n for channel [%s]', vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | baseline corrected | induced\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


%% actual plotting data | TFA timecourse | central stimuli
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-5 5];

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
pl.base_i = dsearchn(TFA.time', pl.base');

% pl.sub2plot = find(F.Subs2use<22); % luminance offset
pl.sub2plot = find(F.Subs2use>21); % isoluminant to background

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2'};
pl.RDKidx = [1 2];


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % raw
        pl.data_ind(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),3));
        pl.data_evo(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),3));
        % baseline corrected
        pl.data_ind_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_ind(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),3)),2))-1);
        pl.data_evo_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_evo(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),3)),2))-1);
    end   
end

% collapse across RDKs
t.conidxRDK1 = [1:6]; t.conidxRDK2 = [2 1 4 3 6 5];
t.data = pl.data_ind(:,:,t.conidxRDK1,1,:); t.data(:,:,:,:,:,2)= pl.data_ind(:,:,t.conidxRDK2,2,:);
pl.data_ind_coll = squeeze(mean(t.data,6));
t.data = pl.data_evo(:,:,t.conidxRDK1,1,:); t.data(:,:,:,:,:,2)= pl.data_evo(:,:,t.conidxRDK2,2,:);
pl.data_evo_coll = squeeze(mean(t.data,6));
t.data = pl.data_ind_bc(:,:,t.conidxRDK1,1,:); t.data(:,:,:,:,:,2)= pl.data_ind_bc(:,:,t.conidxRDK2,2,:);
pl.data_ind_bc_coll = squeeze(mean(t.data,6));
t.data = pl.data_evo_bc(:,:,t.conidxRDK1,1,:); t.data(:,:,:,:,:,2)= pl.data_evo_bc(:,:,t.conidxRDK2,2,:);
pl.data_evo_bc_coll = squeeze(mean(t.data,6));


% raw
pl.conunique = unique(F.conRDKcolattended_label(:,1));
pl.data = nan([size(pl.data_ind_coll,[1 2]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_ind_coll(:,:,t.idx,:),[3,4]);
end

pl.clims1=[0 1].*repmat(max(pl.data,[],"all"),numel(pl.conunique),1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | induced | %s\n for channel [%s]',pl.conunique{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conunique{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


pl.data = nan([size(pl.data_ind_coll,[1 2]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_coll(:,:,t.idx,:),[3,4]);
end
pl.clims1=[0 1].*repmat(max(pl.data,[],"all"),numel(pl.conunique),1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | corrected | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conunique{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});



% baseline corrected
pl.conunique = unique(F.conRDKcolattended_label(:,1));
pl.data = nan([size(pl.data_ind_coll,[1 2]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_ind_bc_coll(:,:,t.idx,:),[3,4]);
end
pl.clims1=[-1 1].*repmat(max(abs(pl.data),[],"all"),numel(pl.conunique),1);

figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | induced | %s\n for channel [%s]',pl.conunique{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conunique{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.conunique = unique(F.conRDKcolattended_label(:,1));
pl.data = nan([size(pl.data_ind_coll,[1 2]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_bc_coll(:,:,t.idx,:),[3,4]);
end
pl.clims1=[-1 1].*repmat(max(abs(pl.data),[],"all"),numel(pl.conunique),1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conunique{i_con};'frequency in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});



%% actual plotting data | TFA timecourse lineplot | central stimuli
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.1 0.1];

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-500 -0];
pl.base_i = dsearchn(TFA.time', pl.base');

% pl.sub2plot = find(F.Subs2use<22); % luminance offset
pl.sub2plot = find(F.Subs2use>21); % isoluminant to background

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2'};
pl.RDKidx = [1 2];

% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % raw
        pl.data_ind(:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),[3,1]));
        pl.data_evo(:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),[3,1]));
        % baseline corrected
        pl.data_ind_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_ind(:,:,i_RDK,i_sub), ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),[3,1,2]))')...
            -1);
        pl.data_evo_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_evo(:,:,i_RDK,i_sub), ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot(i_sub)),[3,1,2]))')...
            -1);
    end
end

% collapse across RDKs
t.conidxRDK1 = [1:6]; t.conidxRDK2 = [2 1 4 3 6 5];
t.data = pl.data_ind(:,t.conidxRDK1,1,:); t.data(:,:,:,:,2)= pl.data_ind(:,t.conidxRDK2,2,:);
pl.data_ind_coll = squeeze(mean(t.data,5));
t.data = pl.data_evo(:,t.conidxRDK1,1,:); t.data(:,:,:,:,2)= pl.data_evo(:,t.conidxRDK2,2,:);
pl.data_evo_coll = squeeze(mean(t.data,5));
t.data = pl.data_ind_bc(:,t.conidxRDK1,1,:); t.data(:,:,:,:,2)= pl.data_ind_bc(:,t.conidxRDK2,2,:);
pl.data_ind_bc_coll = squeeze(mean(t.data,5));
t.data = pl.data_evo_bc(:,t.conidxRDK1,1,:); t.data(:,:,:,:,2)= pl.data_evo_bc(:,t.conidxRDK2,2,:);
pl.data_evo_bc_coll = squeeze(mean(t.data,5));



% actual plotting
pl.conunique = unique(F.conRDKcolattended_label(:,1));
pl.data = nan([size(pl.data_ind_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_ind_coll(:,t.idx,:),[2]);
end

pl.data = pl.data(:);
pl.cons = permute(repmat(pl.conunique',[size(pl.data_ind_coll,1),1,size(pl.data_ind_coll,3)]),[1,3,2]); pl.cons = pl.cons(:);
pl.times = repmat(TFA.time(pl.xlims_i(1):pl.xlims_i(2))',[1 numel(pl.sub2plot) numel(pl.conunique)]); pl.times = pl.times(:);

clear g
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (induced analysis)');
g(1,1).axe_property('YLim',[4 6]);


pl.data = nan([size(pl.data_evo_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_coll(:,t.idx,:),[2]);
end
pl.data = pl.data(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(2,1).axe_property('YLim',[1 2.5]);

pl.data = nan([size(pl.data_ind_bc_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_ind_bc_coll(:,t.idx,:),[2]);
end
pl.data = pl.data(:);
g(1,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,2).stat_summary('type','sem');
g(1,2).set_names('x','time in ms','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs from baseline (induced analysis)');
g(1,2).axe_property('YLim',[-5 5]);

pl.data = nan([size(pl.data_evo_bc_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_bc_coll(:,t.idx,:),[2]);
end
pl.data = pl.data(:);
g(2,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,2).stat_summary('type','sem');
g(2,2).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,2).axe_property('YLim',[-20 40]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();


clear g
pl.data = nan([size(pl.data_evo_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_coll(:,t.idx,:),[2]);
end
pl.data = pl.data(:);
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(1,1).axe_property('YLim',[1 2.5]);

pl.data = nan([size(pl.data_evo_bc_coll,[1 3]),numel(pl.conunique)]);
for i_con = 1:numel(pl.conunique)
    t.idx = strcmpi(F.conRDKcolattended_label(:,1),pl.conunique{i_con});
    pl.data(:,:,i_con) = mean(pl.data_evo_bc_coll(:,t.idx,:),[2]);
end
pl.data = pl.data(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV/m²','color','RDK');
g(2,1).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,1).axe_property('YLim',[-20 40]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 700 700]);
g.draw();



%% actual plotting data | TFA timecourse lineplot | peripheral stimuli
% plotting parameters
% large center as in tango | periphery: central and lateral smaller (fits topographies)
pl.elec2plot = {{'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'}, 'left'; ...
    {'P7';'P9';'PO7';'PO4';'POz';'Oz';'O2'}, 'right'};
% cluster analysis
pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.freqrange=[-0.1 0.1];

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-500 -0];
pl.base_i = dsearchn(TFA.time', pl.base');

% pl.sub2plot = find(F.Subs2use<22); % luminance offset
pl.sub2plot = find(F.Subs2use>21); % isoluminant to background

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK3';'RDK4';'RDK5'};
pl.RDKidx = [3 4 5];

% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % index side of RDK and electrode cluster
        t.elecidx = strcmp(pl.elec2plot(:,2),TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx(i_RDK)).poslabel);
        % raw
        pl.data_ind(:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i{t.elecidx},:,pl.sub2plot(i_sub)),[3,1]));
        pl.data_evo(:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i{t.elecidx},:,pl.sub2plot(i_sub)),[3,1]));
        % baseline corrected
        pl.data_ind_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_ind(:,:,i_RDK,i_sub), ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i{t.elecidx},:,pl.sub2plot(i_sub)),[3,1,2]))')...
            -1);
        pl.data_evo_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_evo(:,:,i_RDK,i_sub), ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i{t.elecidx},:,pl.sub2plot(i_sub)),[3,1,2]))')...
            -1);
    end
end
% collapse across RDKs first
pl.conunique = unique(F.conRDKcolattended_label2(F.conRDKpresent & repmat([false false true true true],size(F.conRDKcolattended_label2,1),1)));
[pl.data_ind_coll pl.data_evo_coll pl.data_ind_bc_coll pl.data_evo_bc_coll] = deal(nan([size( pl.data_ind,1) numel(pl.conunique) numel(pl.RDKidx) numel(pl.sub2plot)]));
for i_con = 1:numel(pl.conunique)
    for i_rdk = 1:numel(pl.RDKidx)
        % index condition
        t.cidx = strcmp(F.conRDKcolattended_label2(:,pl.RDKidx(i_rdk)),pl.conunique{i_con});
        % induced
        pl.data_ind_coll(:,i_con,i_rdk,:) = mean(pl.data_ind(:,t.cidx,i_rdk,:),2);
        % evoked
        pl.data_evo_coll(:,i_con,i_rdk,:) = mean(pl.data_evo(:,t.cidx,i_rdk,:),2);
        % induced bc
        pl.data_ind_bc_coll(:,i_con,i_rdk,:) = mean(pl.data_ind_bc(:,t.cidx,i_rdk,:),2);
        % evoked bc
        pl.data_evo_bc_coll(:,i_con,i_rdk,:) = mean(pl.data_evo_bc(:,t.cidx,i_rdk,:),2);
    end
end

% average across RDKs
pl.data_ind_coll=squeeze(mean(pl.data_ind_coll,3,"omitnan"));
pl.data_evo_coll=squeeze(mean(pl.data_evo_coll,3,"omitnan"));
pl.data_ind_bc_coll=squeeze(mean(pl.data_ind_bc_coll,3,"omitnan"));
pl.data_evo_bc_coll=squeeze(mean(pl.data_evo_bc_coll,3,"omitnan"));



pl.data = pl.data_ind_coll(:);
pl.cons = repmat(pl.conunique',[size(pl.data_ind_coll,1),1,size(pl.data_ind_coll,3)]); pl.cons = pl.cons(:);
pl.times = repmat(TFA.time(pl.xlims_i(1):pl.xlims_i(2))',[1  numel(pl.conunique) numel(pl.sub2plot)]); pl.times = pl.times(:);

clear g
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (induced analysis)');
g(1,1).axe_property('YLim',[5 8]);

pl.data = pl.data_evo_coll(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(2,1).axe_property('YLim',[1 2.5]);

pl.data = pl.data_ind_bc_coll(:);
g(1,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,2).stat_summary('type','sem');
g(1,2).set_names('x','time in ms','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs from baseline (induced analysis)');
g(1,2).axe_property('YLim',[-15 5]);

pl.data = pl.data_evo_bc_coll(:);
g(2,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,2).stat_summary('type','sem');
g(2,2).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,2).axe_property('YLim',[-25 25]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

clear g
pl.data = pl.data_evo_coll(:);
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(1,1).axe_property('YLim',[1 2.5]);

pl.data = pl.data_evo_bc_coll(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,1).axe_property('YLim',[-25 25]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 700 700]);
g.draw();



% average across conditions
pl.conunique2 = {'attended';'irrelevant';'unattended'}
pl.data_ind_coll2=mean(pl.data_ind_coll(:,1:2,:),2); pl.data_ind_coll2(:,2,:)=mean(pl.data_ind_coll(:,3:4,:),2); pl.data_ind_coll2(:,3,:)=mean(pl.data_ind_coll(:,5:6,:),2);
pl.data_evo_coll2=mean(pl.data_evo_coll(:,1:2,:),2); pl.data_evo_coll2(:,2,:)=mean(pl.data_evo_coll(:,3:4,:),2); pl.data_evo_coll2(:,3,:)=mean(pl.data_evo_coll(:,5:6,:),2);
pl.data_ind_bc_coll2=mean(pl.data_ind_bc_coll(:,1:2,:),2); pl.data_ind_bc_coll2(:,2,:)=mean(pl.data_ind_bc_coll(:,3:4,:),2); pl.data_ind_bc_coll2(:,3,:)=mean(pl.data_ind_bc_coll(:,5:6,:),2);
pl.data_evo_bc_coll2=mean(pl.data_evo_bc_coll(:,1:2,:),2); pl.data_evo_bc_coll2(:,2,:)=mean(pl.data_evo_bc_coll(:,3:4,:),2); pl.data_evo_bc_coll2(:,3,:)=mean(pl.data_evo_bc_coll(:,5:6,:),2);


pl.data = pl.data_ind_coll2(:);
pl.cons = repmat(pl.conunique2',[size(pl.data_ind_coll2,1),1,size(pl.data_ind_coll2,3)]); pl.cons = pl.cons(:);
pl.times = repmat(TFA.time(pl.xlims_i(1):pl.xlims_i(2))',[1  numel(pl.conunique2) numel(pl.sub2plot)]); pl.times = pl.times(:);

clear g
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (induced analysis)');
g(1,1).axe_property('YLim',[5 8]);

pl.data = pl.data_evo_coll2(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(2,1).axe_property('YLim',[1 2.5]);

pl.data = pl.data_ind_bc_coll2(:);
g(1,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,2).stat_summary('type','sem');
g(1,2).set_names('x','time in ms','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs from baseline (induced analysis)');
g(1,2).axe_property('YLim',[-15 5]);

pl.data = pl.data_evo_bc_coll2(:);
g(2,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,2).stat_summary('type','sem');
g(2,2).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,2).axe_property('YLim',[-25 25]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

clear g
pl.data = pl.data_evo_coll2(:);
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
g(1,1).axe_property('YLim',[1 2.5]);


pl.data = pl.data_evo_bc_coll2(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
g(2,1).axe_property('YLim',[-25 25]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 700 700]);
g.draw();




