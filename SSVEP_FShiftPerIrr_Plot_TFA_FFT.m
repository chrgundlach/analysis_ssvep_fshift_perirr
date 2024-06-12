%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\EEG\TFA'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:50,'UniformOutput',false)';
F.Subs2use              = [1:13 15:18];
                        
F.TFA.baseline          = [-500 -250];

% F.SSVEP_Freqs           = [20 24 30];
F.RDK_pos               = [0 -255 255];
F.RDK_pos_label         = {'center';'left';'right'};



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




%% plot grand mean FFT data | spectra
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));


pl.time2plot = [2];
pl.freq2plot=F.SSVEP_Freqs;
pl.freqrange=[-0.1 0.1];
pl.sub2plot = F.Subs2use;

% extract data
pl.data = [];
pl.data = mean(mean(mean(TFA.fftdata_ind(:,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot),2),4),3);

% plotting
figure;
set(gcf,'Position',[100 100 600 600],'PaperPositionMode','auto')
subplot(2,1,1)
plot(TFA.fftfreqs,squeeze(pl.data),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,squeeze(mean(pl.data,5)),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('induced GrandMean FFT spectra | N = %1.0f',numel(pl.sub2plot)),'Interpreter','none')
vline(pl.freq2plot,'k:')
% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});

subplot(2,1,2)
pl.data = mean(mean(mean(TFA.fftdata_evo(:,pl.elec2plot_i,:,pl.time2plot,pl.sub2plot),2),4),3);
hold on;
plot(TFA.fftfreqs,squeeze(pl.data),'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(TFA.fftfreqs,squeeze(mean(pl.data,5)),'Color','k','LineWidth',2)
xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f',numel(pl.sub2plot)),'Interpreter','none')
vline(pl.freq2plot,'k:')




%% plot Grand Mean FFT data | topoplot for positions (as frequencies are random)
pl.time2plot = [1:3];
pl.time2plot = [1];
pl.pos2plot='center';
pl.freqrange=[-0.1 0.1];
pl.sub2plot = 1:numel(F.Subs2use);

% extract data
pl.data_ind = []; pl.data_evo = []; pl.title = {}; h.s =[];
for i_sub = 1:numel(pl.sub2plot)
    
end


% extract data
pl.data_ind = []; pl.data_evo = []; pl.title = {}; h.s =[];
for i_freq = 1:numel(pl.freq2plot)
    t.fidx = dsearchn(TFA.fftfreqs',pl.freq2plot(i_freq)+pl.freqrange');
    pl.data_ind(:,:,i_freq) = squeeze(mean(mean(mean(TFA.fftdata_ind(t.fidx(1):t.fidx(2),:,:,pl.time2plot,pl.sub2plot),1),4),3));
    pl.data_evo(:,:,i_freq) = squeeze(mean(mean(mean(TFA.fftdata_evo(t.fidx(1):t.fidx(2),:,:,pl.time2plot,pl.sub2plot),1),4),3));
    pl.title{i_freq} = sprintf('%1.2f [%1.2f %1.2f] Hz', pl.freq2plot(i_freq),pl.freqrange);
end
pl.data_ind(:,:,i_freq+1)=mean(pl.data_ind(:,:,1:i_freq),3);
pl.data_evo(:,:,i_freq+1)=mean(pl.data_evo(:,:,1:i_freq),3);
pl.title{i_freq+1}='average freqs';

figure;
set(gcf,'Position',[100 100 1200 300],'PaperPositionMode','auto')
for i_freq = 1:size(pl.data_ind,3)
    h.s(i_freq,1)=subplot(2,size(pl.data_ind,3),i_freq);
    t.pldata = squeeze(mean( pl.data_ind,2));
%     topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:))],'conv','on','colormap',fake_parula,...
%         'whitebk','on');
    topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:,i_freq))],'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(sprintf('induced %s',pl.title{i_freq}))
    colorbar
    
    h.s(i_freq,2)=subplot(2,size(pl.data_ind,3),i_freq+size(pl.data_ind,3));
    t.pldata = squeeze(mean( pl.data_evo,2));
%     topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
%         'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:))],'conv','on','colormap',fake_parula,...
%         'whitebk','on');
    topoplot( t.pldata(:,i_freq), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(t.pldata(:,i_freq))],'conv','on','colormap',fake_parula,...
        'whitebk','on');
    title(sprintf('evoked %s',pl.title{i_freq}))
    colorbar
end

% % % draw colorbar into new axis
% t.pos2 = get(h.s{i_pl,1},'Position');
% t.pos3 = get(h.s{i_pl,1},'OuterPosition');
% h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
% caxis(pl.clims(1,:));
% colormap(gca,'jet')
% h.cb{i_pl,1}=colorbar;
% t.pos4 = get(h.cb{i_pl,1},'Position');
% set(h.cb{i_pl,1},'Position',[t.pos4(1)+0.04 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3)/2 t.pos4(4)*(2/3)])


%% extract amplitude values for FFT
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'P7';'P9';'P6';'P8';'P10';'PO3';'PO7';'PO4';'PO8';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));


pl.time2plot = [1:4];
pl.freqrange=[-0.1 0.1];
% pl.freqrange=[0 0];
pl.sub2plot = F.Subs2use;
% pl.sub2plot = [7];

pl.data_ind = []; pl.data_evo = [];
pl.conlabel = {'attend RDK1';'attend RDK2'};
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);

% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.fftfreqs', x+pl.freqrange'),[TFA.RDK(i_sub).RDK.freq],'UniformOutput',false));
    for i_RDK = 1:size(t.fidx,2)
        pl.data_ind(i_RDK,:,:,i_sub) = mean(mean(TFA.fftdata_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.elec2plot_i,:,:,i_sub),1),2);
        pl.data_evo(i_RDK,:,:,i_sub) = mean(mean(TFA.fftdata_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.elec2plot_i,:,:,i_sub),1),2);
    end        
end
% into long format
pl.RDK.data_ind_l = pl.data_ind(:);
pl.RDK.data_evo_l = pl.data_evo(:);
% pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, mean(pl.data_ind(:,:,1,:),2))-1);
% pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, mean(pl.data_evo(:,:,1,:),2))-1);
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, pl.data_ind(:,:,1,:))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, pl.data_evo(:,:,1,:))-1);
% pl.data_ind_bc = bsxfun(@minus, pl.data_ind, pl.data_ind(:,:,1,:));
% pl.data_evo_bc = bsxfun(@minus, pl.data_evo, pl.data_evo(:,:,1,:));
pl.RDK.data_ind_bc_l = pl.data_ind_bc(:);
pl.RDK.data_evo_bc_l = pl.data_evo_bc(:);
pl.RDK.subs = repmat(pl.sub2plot,prod(size(pl.data_ind,[1 2 3])),1); pl.RDK.subs = pl.RDK.subs(:);
pl.RDK.RDK = repmat(pl.RDKlabel,1,prod(size(pl.data_ind,[2 3 4]))); pl.RDK.RDK=pl.RDK.RDK(:);
pl.RDK.con = repmat(pl.conlabel',size(pl.data_ind,1),prod(size(pl.data_ind,[3 4]))); pl.RDK.con = pl.RDK.con(:);
pl.RDK.time = repmat(pl.timelabel',prod(size(pl.data_ind,[1 2])),prod(size(pl.data_ind,[4]))); pl.RDK.time = pl.RDK.time(:);


% plot 1
clear g
% plot with gramm
g(1,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_ind_l,'color',pl.RDK.con);
g(1,1).facet_grid([],pl.RDK.RDK);
g(1,1).stat_summary('geom',{'bar','black_errorbar'});
g(1,1).set_names('x','time windows','y','amplitude in \muV','color','attention','column','SSVEP');
g(1,1).set_title('amplitude values of SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

% % troubleshooting
% [pl.RDK.time pl.RDK.con pl.RDK.RDK num2cell(pl.RDK.data_ind_l)]

g(2,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_evo_l,'color',pl.RDK.con);
g(2,1).facet_grid([],pl.RDK.RDK);
g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(2,1).set_names('x','time windows','y','amplitude in \muV','color','attention','column','SSVEP');
g(2,1).set_title('amplitude values of SSVEPs (evoked analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

% plot baseline corrected data (modulation in %)
% plot with gramm
g(1,2)=gramm('x',pl.RDK.time,'y',pl.RDK.data_ind_bc_l,'color',pl.RDK.con);
g(1,2).facet_grid([],pl.RDK.RDK);
g(1,2).stat_summary('geom',{'bar','black_errorbar'});
g(1,2).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(1,2).set_title('amplitude modulations of SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

g(2,2)=gramm('x',pl.RDK.time,'y',pl.RDK.data_evo_bc_l,'color',pl.RDK.con);
g(2,2).facet_grid([],pl.RDK.RDK);
g(2,2).stat_summary('geom',{'bar','black_errorbar'});
g(2,2).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(2,2).set_title('amplitude modulations of SSVEPs (evoked analysis)');
% g.stat_boxplot();
g.set_text_options('base_size',8);
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

% plot 2
% collapse across RDKs
t.data = pl.data_ind(1,:,:,:); t.data(:,:,:,:,2)= pl.data_ind(2,[2 1],:,:);
pl.data_ind_coll = squeeze(mean(t.data,5));
pl.data_ind_coll(3,:,:)=mean(pl.data_ind(3,:,:,:),2);
t.data = pl.data_evo(1,:,:,:); t.data(:,:,:,:,2)= pl.data_evo(2,[2 1],:,:);
pl.data_evo_coll = squeeze(mean(t.data,5));
pl.data_evo_coll(3,:,:)=mean(pl.data_evo(3,:,:,:),2);
pl.conlabel2 = {'attended';'unattended';'irrelevant'};

pl.coll.data_ind_l = pl.data_ind_coll(:);
pl.coll.data_evo_l = pl.data_evo_coll(:);
pl.data_ind_coll_bc = 100*(bsxfun(@rdivide, pl.data_ind_coll, pl.data_ind_coll(:,1,:))-1);
pl.data_evo_coll_bc = 100*(bsxfun(@rdivide, pl.data_evo_coll, pl.data_evo_coll(:,1,:))-1);
pl.coll.data_ind_bc_l = pl.data_ind_coll_bc(:);
pl.coll.data_evo_bc_l = pl.data_evo_coll_bc(:);
pl.coll.subs = repmat(pl.sub2plot,prod(size(pl.data_ind_coll,[1 2])),1); pl.coll.subs = pl.coll.subs(:);
pl.coll.con = repmat(pl.conlabel2,1,prod(size(pl.data_ind_coll,[2 3]))); pl.coll.con = pl.coll.con(:);
pl.coll.time = repmat(pl.timelabel',prod(size(pl.data_ind_coll,[1])),prod(size(pl.data_ind_coll,[3]))); pl.coll.time = pl.coll.time(:);

% plot with gramm
g(1,1)=gramm('x',pl.coll.time,'y',pl.coll.data_ind_l,'color',pl.coll.con);
g(1,1).stat_summary('geom',{'bar','black_errorbar'});
g(1,1).set_names('x','time windows','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitude values of collapsed SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(2,1)=gramm('x',pl.coll.time,'y',pl.coll.data_evo_l,'color',pl.coll.con);
g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(2,1).set_names('x','time windows','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitude values of collapsed SSVEPs (evoked analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(1,2)=gramm('x',pl.coll.time,'y',pl.coll.data_ind_bc_l,'color',pl.coll.con);
g(1,2).stat_summary('geom',{'bar','black_errorbar'});
g(1,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(2,2)=gramm('x',pl.coll.time,'y',pl.coll.data_evo_bc_l,'color',pl.coll.con);
g(2,2).stat_summary('geom',{'bar','black_errorbar'});
g(2,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs (evoked analysis)');
% g.stat_boxplot();
g.set_text_options('base_size',8);
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

% plot 2b
% collapse across RDKs
% different style
clear g
% plot with gramm
g(1,1)=gramm('x',pl.coll.time,'y',pl.coll.data_ind_l,'color',pl.coll.con);
g(1,1).stat_violin('half',true,'normalization','count','width',2,'fill','transparent');
g(1,1).set_names('x','time windows','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitude values of collapsed SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(2,1)=gramm('x',pl.coll.time,'y',pl.coll.data_evo_l,'color',pl.coll.con);
g(2,1).stat_violin('half',true,'normalization','count','width',2,'fill','transparent');
g(2,1).set_names('x','time windows','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitude values of collapsed SSVEPs (evoked analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(1,2)=gramm('x',pl.coll.time,'y',pl.coll.data_ind_bc_l,'color',pl.coll.con);
g(1,2).stat_violin('half',true,'normalization','count','width',2,'fill','transparent');
g(1,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(2,2)=gramm('x',pl.coll.time,'y',pl.coll.data_evo_bc_l,'color',pl.coll.con);
g(2,2).stat_violin('half',true,'normalization','count','width',2,'fill','transparent');
g(2,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs (evoked analysis)');
% g.stat_boxplot();
g.set_text_options('base_size',8);
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

% plot 3 with single subject
clear g


% plot baseline corrected data (modulation in %)
g(1,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_ind_bc_l,'color',pl.RDK.con);
g(1,1).facet_grid([],pl.RDK.RDK);
% g(1,1).stat_summary('geom',{'bar','black_errorbar'});
g(1,1).geom_point('alpha',0.2);
g(1,1).geom_line('alpha',0.2);
g(1,1).stat_summary('geom','line');
g(1,1).stat_summary('type','ci','geom','black_errorbar','width',0.8);
g(1,1).stat_summary('geom','point');
% g(1,1).stat_glm('geom','area','disp_fit',false);
g(1,1).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(1,1).set_title('amplitude modulations of SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 800 350]);
% g.draw();

g(2,1)=gramm('x',pl.RDK.time,'y',pl.RDK.data_evo_bc_l,'color',pl.RDK.con);
g(2,1).facet_grid([],pl.RDK.RDK);
% g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(2,1).geom_point('alpha',0.2);
g(2,1).geom_line('alpha',0.2);
g(2,1).stat_summary('geom','line');
g(2,1).stat_summary('type','ci','geom','black_errorbar','width',0.8);
g(2,1).stat_summary('geom','point');
% g(1,1).stat_glm('geom','area','disp_fit',false);
g(2,1).set_names('x','time windows','y','modulation in %','color','attention','column','SSVEP');
g(2,1).set_title('amplitude modulations of SSVEPs (evoked analysis)');


g(1,2)=gramm('x',pl.coll.time,'y',pl.coll.data_ind_bc_l,'color',pl.coll.con);
% g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(1,2).geom_point('alpha',0.2,'dodge',0.2);
g(1,2).geom_line('alpha',0.2,'dodge',0.2);
g(1,2).stat_summary('geom','line');
g(1,2).stat_summary('type','ci','geom','black_errorbar','width',0.8,'dodge',0.2);
g(1,2).stat_summary('geom','point','dodge',0.2);
% g(1,1).stat_glm('geom','area','disp_fit',false);
g(1,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs (induced analysis)');
% g.stat_boxplot();
% figure('Position',[100 100 600 250]);
% g.draw();

g(2,2)=gramm('x',pl.coll.time,'y',pl.coll.data_evo_bc_l,'color',pl.coll.con);
% g(2,1).stat_summary('geom',{'bar','black_errorbar'});
g(2,2).geom_point('alpha',0.2,'dodge',0.2);
g(2,2).geom_line('alpha',0.2,'dodge',0.2);
g(2,2).stat_summary('geom','line','dodge',0.2);
g(2,2).stat_summary('type','ci','geom','black_errorbar','width',0.8,'dodge',0.2);
g(2,2).stat_summary('geom','point','dodge',0.2);
% g(1,1).stat_glm('geom','area','disp_fit',false);
g(2,2).set_names('x','time windows','y','modulation in %','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs (evoked analysis)');
% g.stat_boxplot();
g.set_text_options('base_size',8);
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();

% save data
% time resolvedcompounddata

R_Mat.all = [{'amplitude_induced','amplitude_evoked','modulation_induced','modulation_evoked','subjects', 'condition', 'time'};...
    num2cell([pl.coll.data_ind_l pl.coll.data_evo_l pl.coll.data_ind_bc_l pl.coll.data_evo_bc_l pl.coll.subs]) ...
    pl.coll.con pl.coll.time];
t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_FShiftBase\data';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
% xlswrite(fullfile(t.path,sprintf('FFT_Amp_data_collapsed_largeclust_%s.csv',t.datestr)),R_Mat.all)


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




