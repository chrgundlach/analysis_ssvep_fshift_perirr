%% script to analyze behavioral data
clearvars

p.path =                'N:\AllgPsy\experimental_data\2024_FShiftPerIrr\behavior\';
p.path =                '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\behavior\';
p.subs=                 arrayfun(@(x) sprintf('%02.0f',x),1:50,'UniformOutput',false)';
% p.subs2use=             [1:13 15:21];% 
% changed experiment from participant 22 onwards (stimuli isoluminant to background
p.subs2use=             [1:13 15:34];%
p.responsewin =         [0.2 1]; % according to p.targ_respwin from run_posnalpha

p.eventtype =           {'target';'distractor'};
p.SSVEP_Freqsforsubs =  ones(size(p.subs2use)); p.SSVEP_Freqsforsubs(p.subs2use>21)=2;
p.colorlum =            repmat({'offset_to_bckgrd'}, 1, numel(p.subs2use));
p.colorlum(p.SSVEP_Freqsforsubs == 2) = {'isolum__to_bckgrd'};


% p.con1index =           [1 1 2 2];
% p.con1name =            'pre cue task';
% p.con1label =           {'attend fixation cross';'attend fixation cross';'attend both rings';'attend both rings'};
% p.con1label =           {'cross';'cross';'rings';'rings'};
% p.con2index =           [1 2 1 2];
% p.con2name =            'post cue attention';
% p.con2label =           {'attend left';'attend right';'attend left';'attend right'};


%% actual calculation
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
            temp.index=cell2mat(cellfun(@(x) ~isempty(cell2mat({x(:).trialnumber})), temp.data_in{i_fi}.resp.experiment,'UniformOutput',false));
            data_in.resp.experiment(temp.index)=temp.data_in{i_fi}.resp.experiment(temp.index);
            data_in.button_presses.experiment(temp.index)=temp.data_in{i_fi}.button_presses.experiment(temp.index);
        end
    end
    
    %% setup response matrices
    data_all_struct = horzcat(data_in.resp.experiment{:});
    data_all_table = struct2table(data_all_struct);
    clear data_events
    
    % extract single event data
    t.evidx = ~isnan([data_all_struct.eventtype]);
    t.data_eventtype = repmat({''},size([data_all_struct.eventtype]));
    t.mat = [data_all_struct.eventtype];
    t.data_eventtype(t.evidx) = arrayfun(@(x) p.eventtype{x}, t.mat(~isnan(t.mat)),'UniformOutput',false);
    t.data_trialnum = repmat([data_all_struct.trialnumber],2,1); % trialnumber
    t.data_blocknum = repmat([data_all_struct.blocknumber],2,1); % block number
    t.data_condition = repmat([data_all_struct.condition],2,1); % block number
    t.data_eventRDK = [data_all_struct.eventRDK]; % RDK of event
    t.data_colorlum =  repmat(p.colorlum(i_sub),size(t.data_eventRDK));
    t.data_eventcolorname = repmat({''},size(t.data_eventRDK));
    for i_rdk = 1:2
        t.data_eventcolorname(t.data_eventRDK==i_rdk)=deal({data_in.RDK(i_rdk).colnames});
    end
    t.data_eventfreq = nan(size(t.data_eventRDK));
    t.data_eventfreq(t.evidx)=arrayfun(@(x) data_in.RDK(x).freq, t.data_eventRDK(t.evidx));
    t.data_eventnum = repmat([1;2],1,size(t.evidx,2)).*t.evidx;
    t.data_response = repmat({''},size(t.data_eventRDK));
    t.data_response(t.evidx) = [data_all_struct.event_response_type];
    t.data_response_RT = nan(size(t.data_eventRDK));
    t.data_response_RT(t.evidx) = [data_all_struct.event_response_RT];
    t.data_subject = repmat(p.subs2use(i_sub),size(t.data_eventtype));
    %t.data_irr_color = repmat({data_in.RDK(5).colnames},size(t.data_eventtype));
    %t.data_irr_freq = repmat(data_in.RDK(3).freq,size(t.data_eventtype));
    

   
    
    % write into one file
    data_events = cell2struct([num2cell([t.data_subject(t.evidx) t.data_trialnum(t.evidx) t.data_blocknum(t.evidx) t.data_condition(t.evidx)]) t.data_eventtype(t.evidx)  ...
        num2cell(t.data_eventRDK(t.evidx)) t.data_eventcolorname(t.evidx) num2cell([t.data_eventfreq(t.evidx) t.data_eventnum(t.evidx)]) ...
        t.data_response(t.evidx) num2cell(t.data_response_RT(t.evidx)) t.data_colorlum(t.evidx)]',...
        {'participant','trialnumber','blocknumber','condition','eventtype','RDKnumber','eventcolor','eventfrequency','eventnumber',...
        'response','RT', 'colorlum'});
    
    if i_sub == 1
        response_events = data_events;
    else
        response_events = vertcat(response_events,data_events);
    end
    
    %% check for false alarms?
    response_FA(i_sub).subject = p.subs2use(i_sub);
    response_FA(i_sub).FA_proper = sum(strcmpi({data_events.response},'FA_proper'));
    response_FA(i_sub).FA = ...
        sum(cellfun(@(x) sum(cellfun(@(y) sum(y(:,8)~=0),x.presses_t)),data_in.button_presses.experiment))-...
        sum(strcmpi({data_events.response},'hit')|strcmpi({data_events.response},'FA_proper'));
    
end

%% interim statistics
response_events_table = struct2table(response_events);
t.dat = grpstats(response_events_table, ["participant","response"],'numel');
t.dat2 = grpstats(response_events_table, ["participant","response"],'mean',"DataVars",["RT"]);

t.idx = strcmp(t.dat.response,'FA_proper');
[t.dat.participant(1:4:end) t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')) t.dat.GroupCount(strcmp(t.dat.response,'CR')) ...
    100*[t.dat.GroupCount(strcmp(t.dat.response,'FA_proper'))./(t.dat.GroupCount(strcmp(t.dat.response,'CR')) + ...
    t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')))]]
mean(100*[t.dat.GroupCount(strcmp(t.dat.response,'FA_proper'))./(t.dat.GroupCount(strcmp(t.dat.response,'CR')) + ...
    t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')))])

t.dat.GroupCount(strcmp(t.dat.response,'hit'))./(t.dat.GroupCount(strcmp(t.dat.response,'hit'))+t.dat.GroupCount(strcmp(t.dat.response,'miss')))

mean(t.dat2.mean_RT(strcmp(t.dat.response,'hit')))

%% summary statistics/descriptives
response_events_table = struct2table(response_events);
response_FA_table = struct2table(response_FA);
% head(response_events_table)

% general summary mean RT and hit rate
% export data
p.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftperirr\data_in';
p.filename1 = 'behavior_events.csv';
p.filename2 = 'behavior_FAs.csv';
writetable(response_events_table, fullfile(p.path,p.filename1))
writetable(response_FA_table, fullfile(p.path,p.filename2))

