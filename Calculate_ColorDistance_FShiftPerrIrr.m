% script to tranfer RGB colors and calculate distance
clearvars

p.path =              '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\behavior\';
p.subs=                 arrayfun(@(x) sprintf('%02.0f',x),1:60,'UniformOutput',false)';

p.subs2use              = [1:13 15:21];
% changed experiment from participant 22 onwards (stimuli isoluminant to
% background and used other frequencies
% participant 42 has lower trial number
% p.subs2use              = [22:41 43:52];

p.subs2use               = [1:13 15:21 22:41 43:52];
p.subsets               = {[1:13 15:21];[22:41 43:52]};
p.subset_label          = {'luminance offset';'isoluminant'};
p.subset_label_mat      = cellfun(@(x,y) repmat({y},1,numel(x)),p.subsets,p.subset_label,'UniformOutput',false);
p.subset_label_mat      = [p.subset_label_mat{:}];

p.responsewin =         [0.2 1.2]; % according to p.targ_respwin from run_posnalpha

p.cols =                [0.0980 0.0392 0 1; % red
                        0 0.0596 0.1490 1; % blue
                        0 0.0745 0 1]; % green

p.cols = p.cols*6;

p.colnames =            {'redish';'blue';'green'};
p.pos =                 [];




col_RGB = p.cols;

% propixx color values
col.xy_r = [0.665 0.321]; col.xy_g = [0.172 0.726]; col.xy_b = [0.163 0.039]; col.xyz_w = [0.3127 0.3290];
% col_RGB = col_RGB.*[0.5 0.5 0.5; 0.3 0.3 0.3; 0.4 0.4 0.4];
% [xyz] = rgb2xyz_custom([1 0 0],col.xy_r, col.xy_g, col.xy_b, col.xyz_w);
% [RGB] = xyz2rgb_custom([0.6650    0.3210    0.0140],col.xy_r, col.xy_g, col.xy_b, col.xyz_w);


%% extract data for each subject and calculate distances
for i_sub = 1:numel(p.subs2use)
    % load data
    temp.files = dir(sprintf('%sVP%s_timing*.mat',p.path,p.subs{p.subs2use(i_sub)}));
    data_in.resp.experiment = repmat({[nan]},1,8);
    data_in.button_presses.experiment = repmat({[nan]},1,8);
    for i_fi = 1:numel(temp.files)
        temp.data_in{i_fi}=load(sprintf('%s%s',p.path,temp.files(i_fi).name));
        % extract relevant data
        try data_in.conmat.experiment = temp.data_in{i_fi}.randmat.experiment;
            data_in.RDK = temp.data_in{i_fi}.RDK.RDK;
        end
        if any(strcmp(fieldnames(temp.data_in{i_fi}.resp),'experiment'))
            temp.index=cell2mat(cellfun(@(x) ~isempty(x), temp.data_in{i_fi}.resp.experiment,'UniformOutput',false));
            data_in.resp.experiment(temp.index)=temp.data_in{i_fi}.resp.experiment(temp.index);
            data_in.button_presses.experiment(temp.index)=temp.data_in{i_fi}.button_presses.experiment(temp.index);
        end
    end
    
    %% get color values for RDKs
    t.colmat = cell2mat(cellfun(@(x) x(1,1:3),{data_in.RDK.col},'UniformOutput',false)');
    t.colmat_base = cell2mat(cellfun(@(x) x(1,1:3),{data_in.RDK.col_init},'UniformOutput',false)');
    t.colmat_bckgrnd = data_in.RDK(1).col(2,1:3);
    
    t.colmat_name = {data_in.RDK.colnames}';
    t.RDKmat = {'RDK1';'RDK2';'RDK3';'RDK4';'RDK5'};
    t.RDKposmat = {'center';'center';'peri';'peri';'peri'};
    t.groupmat = repmat(find(cellfun(@(x) any(ismember(x,p.subs2use(i_sub))),p.subsets)),numel(t.RDKmat),1);
    t.colmat_XYZ = rgb2xyz_custom(t.colmat, col.xy_r, col.xy_g, col.xy_b, col.xyz_w );
    t.colmat_LAB = xyz2lab(t.colmat_XYZ);
    t.colmat_LAB_deg = rad2deg(atan2(t.colmat_LAB(:,3),t.colmat_LAB(:,2)));
    t.colmat_LAB_deg2 = reshape(repmat(t.colmat_LAB_deg',2,1),numel(t.colmat_name)*2,[]);
    t.colmat_base_XYZ = rgb2xyz_custom(t.colmat_base, col.xy_r, col.xy_g, col.xy_b, col.xyz_w );
    t.colmat_base_LAB = xyz2lab(t.colmat_base_XYZ);
    t.colmat_bckgrnd_XYZ = rgb2xyz_custom(t.colmat_bckgrnd, col.xy_r, col.xy_g, col.xy_b, col.xyz_w );
    t.colmat_bckgrnd_LAB = xyz2lab(t.colmat_bckgrnd_XYZ);
    
    if i_sub == 1
        colmat = cell2struct([repmat({p.subs2use(i_sub)},numel(t.RDKmat),1) num2cell(t.groupmat) t.RDKmat t.RDKposmat t.colmat_name ...
            num2cell(t.colmat_LAB_deg) num2cell(t.colmat_LAB(:,1))]',...
            {'participant','group', 'RDK_id', 'RDK_pos','RDK_color', 'RDK_color_LAB_deg', 'RDK_color_LAB_Lum'});
        colmat_RGB = t.colmat;
        colmat_XYZ = t.colmat_XYZ;
        colmat_lab = t.colmat_LAB;
        colmat_base_rgb = t.colmat_base;
        colmat_base_LAB =t.colmat_base_LAB;
        colmat_name = t.colmat_name;
        colmat_bckgrnd_LAB = t.colmat_bckgrnd_LAB;
    else
        colmat = [colmat; ...
            cell2struct([repmat({p.subs2use(i_sub)},numel(t.RDKmat),1) num2cell(t.groupmat) t.RDKmat t.RDKposmat t.colmat_name ...
            num2cell(t.colmat_LAB_deg) num2cell(t.colmat_LAB(:,1))]',...
            {'participant','group', 'RDK_id', 'RDK_pos','RDK_color', 'RDK_color_LAB_deg', 'RDK_color_LAB_Lum'})];
        colmat_RGB = [colmat_RGB; t.colmat];
        colmat_XYZ = [colmat_XYZ; t.colmat_XYZ];
        colmat_lab = [colmat_lab; t.colmat_LAB];
        colmat_base_rgb = [colmat_base_rgb; t.colmat_base];
        colmat_base_LAB = [colmat_base_LAB; t.colmat_base_LAB];
        colmat_name = [colmat_name; t.colmat_name];
        colmat_bckgrnd_LAB = [colmat_bckgrnd_LAB; t.colmat_bckgrnd_LAB];
    end
    
end

%% graphical checking
% figure; histogram([colmat.dist2attended],20)

% plot lab space
% plot_Lab(1,colmat_lab',1,colmat_base_rgb,15,0)
figure; polarhistogram(atan2(colmat_lab(:,3),colmat_lab(:,2)),90)
% figure; histogram([colmat.dist2attended_deg],150)
% figure; histogram([colmat.dist2unattended_deg],150)


pl.idx = 1;
pl.pos = unique(t.RDKposmat);
figure;
set(gcf,'Position',[100 100 900 600],'PaperPositionMode','auto')
for i_group = 1:numel(p.subsets)
    for i_pos = 1:2
        subplot(numel(p.subsets),2, pl.idx)
        pl.data = [];
        % index of group
        t.idxgr = [colmat.group]==i_group;
        t.idxpos = strcmp({colmat.RDK_pos},pl.pos{i_pos });
        for i_col = 1:3
            t.idxcol = strcmp(colmat_name,p.colnames{ i_col});
            pl.data{i_col} = colmat_lab(t.idxgr&t.idxpos&t.idxcol',1);
        end

        for i_col = 1:3
             histogram(pl.data{i_col},20:0.5:70,'FaceColor',p.cols(i_col,1:3)); hold on
        end
        title(sprintf('L values of CIE L*a*b | %s | %s',pl.pos{i_pos},p.subset_label{i_group}));
        pl.idx = pl.idx +1;
    end
end


% Calculate average CIE LAB values for each color
[pl.data_avg, pl.data_sd, pl.data_avg_rgb, pl.data_sd_rgb, pl.data_avg_lab_bcgr, pl.data_sd_lab_bcgr] = deal([]);
for i_group = 1:numel(p.subsets)
    for i_pos = 1:2
        for i_col = 1:3
            pl.idxcol = strcmp(colmat_name, p.colnames{i_col});
            pl.idxgr = [colmat.group]==i_group;
            pl.idxpos = strcmp({colmat.RDK_pos},pl.pos{i_pos });
            pl.data_avg(i_col,i_group,i_pos, :) = mean(colmat_lab(pl.idxgr&pl.idxpos&pl.idxcol', :), 1);
            pl.data_sd(i_col,i_group,i_pos, :) = std(colmat_lab(pl.idxgr&pl.idxpos&pl.idxcol', :), 1);
            pl.data_avg_rgb(i_col,i_group,i_pos, :) = mean(colmat_RGB(pl.idxgr&pl.idxpos&pl.idxcol', :), 1);
            pl.data_sd_rgb(i_col,i_group,i_pos, :) = std(colmat_RGB(pl.idxgr&pl.idxpos&pl.idxcol', :), 1);
        end
    end
    % background
    pl.idxgr2 = strcmp(p.subset_label_mat, p.subset_label{i_group});
    pl.data_avg_lab_bcgr(i_group, :) = mean(colmat_bckgrnd_LAB(pl.idxgr2, :), 1);
    pl.data_sd_lab_bcgr(i_group, :) = std(colmat_bckgrnd_LAB(pl.idxgr2, :), 1);
end

% Display the average CIE LAB values
fprintf('Average CIE L*a*b values:\n');
for i_group = 1:numel(p.subsets)
    for i_pos = 1:2
        for i_col = 1:3
            pl.coldata = [squeeze(pl.data_avg(i_col,i_group,i_pos, :)) squeeze(pl.data_sd(i_col,i_group,i_pos, :))]';
            fprintf('%s | %s | %s: L = %1.3f(%1.3f), a = %1.3f(%1.3f), b = %1.3f(%1.3f)\n', ...
                p.subset_label{i_group},pl.pos{i_pos},p.colnames{i_col}, ...
                pl.coldata);
        end
    end
    pl.coldata = [squeeze(pl.data_avg_lab_bcgr(i_group, :)) squeeze(pl.data_sd_lab_bcgr(i_group, :))]';
    fprintf('%s | background: L = %1.3f(%1.3f), a = %1.3f(%1.3f), b = %1.3f(%1.3f)\n', ...
        p.subset_label{i_group},pl.coldata);
end


% display three dimensional data
% figure
% plotChromaticity()
% 
% 
% gamut = CIELabGamut();
% PlotRings(gamut);
% figure;
% PlotVolume(gamut);

%%
% figure;
% subplot(1,3,1); scatter(pl.data{:,1},pl.data{:,2}); xlabel(p.colnames{1});ylabel(p.colnames{2})
% [corres.r corres.p] = corr(pl.data(:,1),pl.data(:,2),'Type','Spearman');
% title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,1), corres.r, corres.p))
% subplot(1,3,2); scatter(pl.data(:,1),pl.data(:,3)); xlabel(p.colnames{1});ylabel(p.colnames{3})
% [corres.r, corres.p] = corr(pl.data(:,1),pl.data(:,3),'Type','Spearman');
% title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,1), corres.r, corres.p))
% subplot(1,3,3); scatter(pl.data(:,2),pl.data(:,3)); xlabel(p.colnames{2});ylabel(p.colnames{3})
% [corres.r, corres.p] = corr(pl.data(:,2),pl.data(:,3),'Type','Spearman');
% title(sprintf('L values; Rho(%1.0f) = %1.3f; p = %1.3f', size(pl.data,2), corres.r, corres.p))


figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = rad2deg(atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2))-mean(atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2))));
    histogram(pl.data(:,i_col),linspace(-0.5,0.5,30),'FaceColor',p.cols(i_col,1:3)); hold on
end
title('hue angle differences of CIE L*a*b space for different colors');


figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = atan2(colmat_lab(pl.idx,3),colmat_lab(pl.idx,2));
    polarhistogram(pl.data(:,i_col),'BinWidth',2*pi/90,'FaceColor',p.cols(i_col,:)); hold on
end
title('hue angles of CIE L*a*b space for isoluminant colors');

figure;
pl.data = [];
for i_col = 1:3
    pl.idx = strcmp(colmat_name,p.colnames{ i_col});
    pl.data(:,i_col) = atan2(colmat_base_LAB(pl.idx,3),colmat_base_LAB(pl.idx,2));
    polarhistogram(pl.data(:,i_col),'BinWidth',2*pi/90,'FaceColor',p.cols(i_col,:)); hold on
end
title('hue angles of CIE L*a*b space for base colors');


pl.data_m = []; pl.data_std = [];
for i_col = 1:3
    pl.data_m(i_col,:) = mean(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
    pl.data_std(i_col,:) = std(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
    fprintf('\n%s Mean(Std) of CIE L*a*b values: %1.3f (%1.3f), %1.3f (%1.3f), %1.3f (%1.3f)',...
        p.colnames{i_col}, reshape([pl.data_m(i_col,:);pl.data_std(i_col,:)],1,[]))
end
fprintf('\n')

% calculate color differences
pl.data_m = [];
for i_col = 1:3
    pl.data_m(i_col,:) = mean(colmat_lab(strcmp(colmat_name,p.colnames{ i_col}),:));
end
pl.data_m_hue = atan2(pl.data_m(:,3),pl.data_m(:,2));
pl.data_m_hue(pl.data_m_hue<0) = 2*pi -abs(pl.data_m_hue(pl.data_m_hue<0));
pl.hue_diffs = abs([pl.data_m_hue(1)-pl.data_m_hue(2) pl.data_m_hue(1)-pl.data_m_hue(3) pl.data_m_hue(2)-pl.data_m_hue(3)]);
pl.hue_diffs(pl.hue_diffs>pi)=(2*pi)-pl.hue_diffs(pl.hue_diffs>pi);
pl.hue_diffs = rad2deg(pl.hue_diffs);
fprintf('range of differences between color hues: %1.3f° %1.3f° %1.3f°\n',pl.hue_diffs)


% plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
%     'SavePath','C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftBase\figures',...
%     'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)
% 
plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
    'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100,'disp_LAB_vals',true)
% plot_colorwheel([0    0.0580    0.1451;     0.0941    0.0376         0;          0    0.0784         0],'ColorSpace','propixxrgb',...
%     'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100,'disp_LAB_vals',true)

% plot_colorwheel([0 1 0; 0 0.4 1; 1 0.4 0],'ColorSpace','propixxrgb',...
%     'SavePath','C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShiftBase\figures',...
%     'LAB_L',75,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)

plot_colorwheel(p.cols,'ColorSpace','propixxrgb',...
    'SavePath','C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_FShift_Probabil\figures',...
    'LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)

%% plot color planes
plot_colorplanes

%% saving
colmat_table = struct2table(colmat);
% general summary mean RT and hit rate
% export data
p.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftprobabil\data';
p.filename1 = 'color_values.csv';
writetable(colmat_table, fullfile(p.path,p.filename1))







