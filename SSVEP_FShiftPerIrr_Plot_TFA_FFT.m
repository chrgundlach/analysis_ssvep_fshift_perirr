%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\EEG\TFA'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:50,'UniformOutput',false)';
% F.Subs2use              = [1:13 15:21];
% changed experiment from participant 22 onwards (stimuli isoluminant to
% background and used other frequencies
F.Subs2use              = [1:13 15:34];
                        
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




%% plot grand mean FFT data | spectra | for distinct frequencies (lookup of respective electrode cluster)
% plotting parameters
pl.elec2plot = {{'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
    {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% pl.elec2plot = {{'POz';'O1';'Oz';'I2';'Iz'}, 'center';...
%     {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
%     {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


pl.time2plot = [1];
pl.sub2plot = 1:numel(F.Subs2use);
% pl.sub2plot = find(F.Subs2use<22); % luminance offset
% pl.freq2plot=F.SSVEP_Freqs{1}(2);

pl.sub2plot = find(F.Subs2use>21); % isoluminant to background
pl.freq2plot=F.SSVEP_Freqs{2}(5);

% extract data
pl.data_ind = nan(size(TFA.fftdata_ind,1), numel(pl.sub2plot)); pl.data_evo = pl.data_ind;
for i_sub = 1:numel(pl.sub2plot)
    % index position of frequency
    t.rdkidx = [TFA.RDK(pl.sub2plot(i_sub)).RDK.freq]==pl.freq2plot;
    t.posidx = TFA.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).poslabel;
    t.elidx = strcmp(pl.elec2plot(:,2), t.posidx);
    % index conditions for which the RDK is shown
    t.conidx = any(t.rdkidx & F.conRDKpresent,2);

    % extract data
    pl.data_ind(:,i_sub) = mean(TFA.fftdata_ind(:,pl.elec2plot_i{t.elidx},t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[2,3,4]);
    pl.data_evo(:,i_sub) = mean(TFA.fftdata_evo(:,pl.elec2plot_i{t.elidx},t.conidx,pl.time2plot,pl.sub2plot(i_sub)),[2,3,4]);
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
ylabel('amplitude in \muV')
title(sprintf('induced GrandMean FFT spectra | N = %1.0f | FOI = %1.1f Hz', ...
    numel(pl.sub2plot), pl.freq2plot),'Interpreter','none')
vline(pl.freq2plot,'k:')
% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i))),'o','r',4,1});

subplot(2,1,2)
hold on;
plot(TFA.fftfreqs,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(TFA.fftfreqs,mean(pl.data_evo,2),'Color','k','LineWidth',2)
xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f | FOI = %1.1f Hz', ...
    numel(pl.sub2plot), pl.freq2plot),'Interpreter','none')
vline(pl.freq2plot,'k:')
box on

% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.28 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i))),'o','r',4,1});




%% plot Grand Mean FFT data | topoplot for positions (as frequencies are random)
pl.time2plot = [1:3];
pl.time2plot = [1];
% pl.pos2plot='center';
% pl.pos2plot='right';
pl.pos2plot='left';
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
        t.data_ind(:,i_rdk,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.idx(1):t.idx(2),:,:,pl.time2plot,i_sub),[1 4]));
        t.data_evo(:,i_rdk,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.idx(1):t.idx(2),:,:,pl.time2plot,i_sub),[1 4]));
        
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



%% extract amplitude values for FFT
% plotting parameters
pl.elec2plot = {{'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'}, 'center';...
    {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
    {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% smaller for center
% pl.elec2plot = {{'POz';'O1';'Oz';'I2';'Iz'}, 'center';...
%     {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'}, 'left'; ...
%     {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'}, 'right'};

% % only center for periphery
% pl.elec2plot = {{'POz';'O1';'Oz';'I2';'Iz'}, 'center';...
%     {'POz';'Oz';'O1';'PO3'}, 'left'; ...
%     {'POz';'Oz';'O2';'PO4'}, 'right'};

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
% pl.RDK.data_ind_bc = bsxfun(@minus,  pl.RDK.data_ind, pl.RDK.data_ind(:,:,1,:));
% pl.RDK.data_evo_bc = bsxfun(@minus, pl.RDK.data_evo, pl.RDK.data_evo(:,:,1,:));


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

R_Mat.all = [{'amplitude_induced','amplitude_evoked','modulation_induced','modulation_evoked','subjects', 'condition', 'time', ...
    'RDK_id', 'RDK_position1','RDK_position2', 'RDK_freq', 'RDK_color','RDK_colorlum', 'RDK_ispresented', 'RDK_isattended', 'RDK_isattended2', 'RDK_electrodes'}; ...
    num2cell([pl.RDK.data_ind(:) pl.RDK.data_evo(:) pl.RDK.data_ind_bc(:) pl.RDK.data_evo_bc(:) pl.RDK.sub(:) pl.RDK.con(:)]) ...
    pl.RDK.timewin(:) pl.RDK.RDK_id(:) pl.RDK.RDK_pos1(:) pl.RDK.RDK_pos2(:) num2cell(pl.RDK.RDK_freq(:)) pl.RDK.RDK_color(:) ...
    pl.RDK.RDK_colorlum(:) pl.RDK.RDK_ispresent(:) pl.RDK.RDK_isattended(:) pl.RDK.RDK_isattended2(:) pl.RDK.RDK_electrodes(:)
    ];

R_Mat.all_table = cell2table(R_Mat.all(2:end,:), "VariableNames",R_Mat.all(1,:));


t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftperirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
% xlswrite(fullfile(t.path,sprintf('FFT_Amp_data_largeclust_%s.csv',t.datestr)),R_Mat.all)
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclust_%s.csv',t.datestr)),'Delimiter',';')
writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclust_allsubs_%s.csv',t.datestr)),'Delimiter',';')
% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_pericenter_allsubs_%s.csv',t.datestr)),'Delimiter',';')

%% actual plotting data | TFA Grand Mean timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
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

pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


pl.data_ind = squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4));
pl.data_evo = squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4));
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, ...
    mean(squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4)),2))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, ...
    mean(squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2plot),3),4)),2))-1);

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


%% actual plotting data | TFA timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-5 5];

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(i_sub).RDK.freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(i_sub).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % raw
        pl.data_ind(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3));
        pl.data_evo(:,:,:,i_RDK,i_sub) = ...
            squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3));
        % baseline corrected
        pl.data_ind_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_ind(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3)),2))-1);
        pl.data_evo_bc(:,:,:,i_RDK,i_sub) = ...
            100*(bsxfun(@rdivide, pl.data_evo(:,:,:,i_RDK,i_sub), ...
            mean(squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3)),2))-1);
    end   
end

% collapse across RDKs
t.data = pl.data_ind(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_ind(:,:,[2 1],2,:);
pl.data_ind_coll = squeeze(mean(t.data,6));
pl.data_ind_coll(:,:,3,:)=mean(pl.data_ind(:,:,:,3,:),3);
t.data = pl.data_evo(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_evo(:,:,[2 1],2,:);
pl.data_evo_coll = squeeze(mean(t.data,6));
pl.data_evo_coll(:,:,3,:)=mean(pl.data_evo(:,:,:,3,:),3);
t.data = pl.data_ind_bc(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_ind_bc(:,:,[2 1],2,:);
pl.data_ind_bc_coll = squeeze(mean(t.data,6));
pl.data_ind_bc_coll(:,:,3,:)=mean(pl.data_ind_bc(:,:,:,3,:),3);
t.data = pl.data_evo_bc(:,:,:,1,:); t.data(:,:,:,:,:,2)= pl.data_evo_bc(:,:,[2 1],2,:);
pl.data_evo_bc_coll = squeeze(mean(t.data,6));
pl.data_evo_bc_coll(:,:,3,:)=mean(pl.data_evo_bc(:,:,:,3,:),3);
pl.conlabel2 = {'attended';'unattended';'irrelevant'};


% raw
pl.data = mean(pl.data_ind_coll,4);
pl.clims1=[0 1].*squeeze(max(max(pl.data)));
pl.clims1=[0 1].*repmat(max(max(max(pl.data))),3,1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | induced | %s\n for channel [%s]',pl.conlabel2{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.data = mean(pl.data_evo_coll,4);
pl.clims1=[0 1].*squeeze(max(max(pl.data)));
pl.clims1=[0 1].*repmat(max(max(max(pl.data))),3,1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
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
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
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
pl.data = mean(pl.data_ind_bc_coll,4);
pl.clims1=[-1 1].*squeeze(max(max(abs(pl.data))));
pl.clims1=[-1 1].*repmat(max(max(max(abs(pl.data)))),3,1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | induced | %s\n for channel [%s]',pl.conlabel2{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.data = mean(pl.data_evo_bc_coll,4);
pl.clims1=[-1 1].*squeeze(max(max(abs(pl.data))));
pl.clims1=[-1 1].*repmat(max(max(max(abs(pl.data)))),3,1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.conlabel2{i_con};'frequency in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});



%% actual plotting data | TFA timecourse lineplot
% plotting parameters
pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.1 0.1];

pl.xlims=[-1000 2000]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-500 -0];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(i_sub).RDK.freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(i_sub).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        % raw
        pl.data_ind(:,:,i_RDK,i_sub) = ...
            squeeze(mean(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3),1));
        pl.data_evo(:,:,i_RDK,i_sub) = ...
            squeeze(mean(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,i_sub),3),1));
        % baseline corrected
        pl.data_ind_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_ind(:,:,i_RDK,i_sub), ...
            squeeze(mean(mean(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3),1),2))')...
            -1);
        pl.data_evo_bc(:,:,i_RDK,i_sub) = ...
            100*(...
            bsxfun(@rdivide, pl.data_evo(:,:,i_RDK,i_sub), ...
            squeeze(mean(mean(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,i_sub),3),1),2))')...
            -1);
    end
end

% collapse across RDKs
t.data = pl.data_ind(:,:,1,:); t.data(:,:,:,:,2)= pl.data_ind(:,[2 1],2,:);
pl.data_ind_coll = squeeze(mean(t.data,5));
pl.data_ind_coll(:,3,:)=mean(pl.data_ind(:,:,3,:),2);
t.data = pl.data_evo(:,:,1,:); t.data(:,:,:,:,2)= pl.data_evo(:,[2 1],2,:);
pl.data_evo_coll = squeeze(mean(t.data,5));
pl.data_evo_coll(:,3,:)=mean(pl.data_evo(:,:,3,:),2);
t.data = pl.data_ind_bc(:,:,1,:); t.data(:,:,:,:,2)= pl.data_ind_bc(:,[2 1],2,:);
pl.data_ind_bc_coll = squeeze(mean(t.data,5));
pl.data_ind_bc_coll(:,3,:)=mean(pl.data_ind_bc(:,:,3,:),2);
t.data = pl.data_evo_bc(:,:,1,:); t.data(:,:,:,:,2)= pl.data_evo_bc(:,[2 1],2,:);
pl.data_evo_bc_coll = squeeze(mean(t.data,5));
pl.data_evo_bc_coll(:,3,:)=mean(pl.data_evo_bc(:,:,3,:),2);
pl.conlabel2 = {'attended';'unattended';'irrelevant'};


% actual plotting
pl.data = pl.data_ind_coll(:);
pl.cons = repmat(pl.conlabel2',size(pl.data_ind_coll,1),size(pl.data_ind_coll,3)); pl.cons = pl.cons(:);
pl.times = repmat(TFA.time(pl.xlims_i(1):pl.xlims_i(2))',size(pl.data_ind_coll,[2 3])); pl.times = pl.times(:);

clear g
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (induced analysis)');

pl.data = pl.data_evo_coll(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');

pl.data = pl.data_ind_bc_coll(:);
g(1,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,2).stat_summary('type','sem');
g(1,2).set_names('x','time in ms','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs from baseline (induced analysis)');

pl.data = pl.data_evo_bc_coll(:);
g(2,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,2).stat_summary('type','sem');
g(2,2).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();




