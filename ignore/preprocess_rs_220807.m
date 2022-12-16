%% Script description
% Script to preprocess all the given raw EEG data.
%
% Author: Fabian Mager

%%
clc, clear all, close all,

addpath('/Volumes/Promise/Shared/Microstates/eeglab/eeglab2021.1');
addpath('eeglab/eeglab2019_1/plugins/dipfit/standard_BEM/elec/');

eeglab

%% Parameters

% Visualization
vis = 1;

% session
visit = {'visit_0','visit_6','visit_104'};

% resample
sr = 128;

% load data in minutes
load_s = 8;

% each epoch has 4 seconds for eeg analysis
epoch = 4;

% number of channels
n_ch = 64;

% filter
filter = [1 40];

% ICA
ica = 1;

% init table
run_sub = 0;
CleanTab = [];

process_date = datestr(date,'yymmdd');

%%

for ic = 1:length(visit)
    
    % path to results
    path_to_results = ['/Volumes/Promise/Shared/Microstates/RS/Preprocessed/'...
        'preprocessed_' process_date '_Filt_'  num2str(filter(1)) '_' num2str(filter(2))...
        '/' visit{ic} '/'];
    
    path_to_results_final = ['/Volumes/Promise/Shared/Microstates/RS/Preprocessed_final/'...
        'preprocessed_' process_date '_Filt_'  num2str(filter(1)) '_' num2str(filter(2))...
        '/' visit{ic} '/'];
    
    % make directory
    mkdir(path_to_results);
    mkdir(path_to_results_final);
    
    diary([path_to_results 'diary.txt']);
    date_now = datetime('now');
    
    fprintf('\n----------------- %s -----------------\n',string(date_now))
    fprintf('Starting preprocessing:\n')
    fprintf('sampling rate: %.0f Hz\n', sr)
    fprintf('epoching: %.0f s\n', epoch)
    fprintf('filter: %.0f to %.0f Hz\n', filter(1), filter(2))
    fprintf('ICA: %.0f\n', ica)
    
    % List of all subjects   
    ID_list = dir(['/Volumes/Promise/Shared/Microstates/RS/RAW/' visit{ic} '/*.bdf']);
    
    % Remove artifacted files
    ID_list([ID_list.bytes]<1e8) = [];
    
    proc_time = zeros(1, numel(ID_list));
    filtorder = []; revfilt = 1;
    
    % init table
    varNames = {'name','cohort','visit',...
        'group','num_epochs','num_bad_raw_chan',...
        'bad_raw_list','num_bad_chan','bad_chan_list',...
        'num_all_chan','all_chan_list','num_rem_epochs',...
        'num_rem_epoch_chan','num_rem_ic','process_time',...
        'perc_dat_int'};
    varTypes = {'string','string','string',...
        'string','double','double',...
        'string','double','string',...
        'double','string','double',...
        'double','double','double',...
        'double'};
    CleanTab = [CleanTab; table('Size',[numel(ID_list) numel(varNames)],'VariableTypes',varTypes,'VariableNames',varNames)];
    
    %% for each subject ...
    for sub = 1:numel(ID_list)

        run_sub = run_sub + 1;
        
        % Start timer
        tstart = tic;
        
        % Get subject ID
        ID = ID_list(sub).name(1:end-4);
        
        % Print subject ID
        fprintf('Processing subject %.0f of %.0f\n',sub,numel(ID_list))
        
        %% Reading the data in
        
        bdfFile = [ID_list(sub).folder '/' ID_list(sub).name];
        
        EEGnoref = pop_biosig(bdfFile,...
            'blockrange', [0 load_s*60]);
                
        EEG = EEGnoref;
               
        EEG = pop_chanedit(EEG, 'lookup','64-4_Biosemi.xyz');
                
        EEG = pop_select(EEG,'channel',1:n_ch);
        
        EEG = eeg_checkset(EEG);
        
        chanlocs = eeg_mergelocs(EEG.chanlocs);
        
        CleanTab.name(run_sub) = ID_list(sub).name(1:end-4);
                              
        %% Filtering and Resampling
        
        % Resample to 256 Hz per second
        EEG = pop_resample(EEG, sr); 
        
        % Bandpass filter - Filter data using Hamming windowed sinc FIR filter
        EEG = pop_eegfiltnew(EEG, filter(1), filter(2));
       
        % Notch filtering
        EEG = pop_eegfiltnew(EEG,48,51,filtorder, revfilt);
 
        %% Bad portions of data
        EEGprebadport = EEG;
        
        EEG.etc.clean_channel_mask = true(EEG.nbchan,1);
        
        % Here we can add channels we know are bad
        num_raw_bad_chan = {'Iz'};
        
        options = {'FlatlineCriterion',4,...
            'LineNoiseCriterion',6,...
            'ChannelCriterion',0.65,...
            'ChannelCriterionMaxBadTime',0.6,...
            'BurstRejection','off',...
            'BurstCriterion','off',...
            'WindowCriterion','off',...
            'Highpass', 'off',...
            'WindowCriterionTolerances','off',...
            };
       
        
        fprintf('Artifact Removal options:\n')
        disp((options(:)));
        
        
        EEG = clean_artifacts(EEG,options{:});
        
        num_raw_bad_chan = unique([num_raw_bad_chan   chanlocs(~EEG.etc.clean_channel_mask).labels]);
        
        CleanTab.num_bad_raw_chan(run_sub) = length(num_raw_bad_chan);
         
        if length(num_raw_bad_chan) > 60
            continue
        end
        
        if vis == 1
            vis_artifacts(EEG,EEGprebadport); end

        %% Load data again, this time average refence, ignore bad channels
        
        EEGnoref = EEG;
        
        [ref_ch, ref_ch_idx] = setdiff({chanlocs.labels},num_raw_bad_chan,'stable');
                
        EEGref = pop_biosig(bdfFile,...
            'ref',ref_ch_idx,...
            'refoptions',{'keepref' 'on'},...
            'blockrange', [0 load_s*60]);
        
        EEG = EEGref;
        
        EEG.subject = ID;
        
        EEG.filename = ID_list(sub).name;
        
        EEG.filepath = ID_list(sub).folder;
        
        % Add channel locations
        EEG = pop_chanedit(EEG, 'lookup','64-4_Biosemi.xyz');
        
        EEG = pop_select(EEG,'channel',1:n_ch);
        
        EEG.etc.clean_channel_mask = false(n_ch,1);
        
        EEG.etc.clean_channel_mask(ref_ch_idx) = true;

        EEG = pop_select(EEG,'nochannel',num_raw_bad_chan);
        
        EEG = eeg_checkset(EEG);
        
        
        %% Filtering and Resampling
        
        % Resample to 256 Hz per second
        EEG = pop_resample(EEG, sr); 

        % Bandpass filter - Filter data using Hamming windowed sinc FIR filter
        EEG = pop_eegfiltnew(EEG, filter(1), filter(2));
   
        % Notch filtering
        EEG = pop_eegfiltnew(EEG,48,51,filtorder, revfilt);
    
        EEGpreart = EEG;
        
        %% Artefact Removal
        
        % Clean data from artifacts
        %   ChannelCriterion : Minimum channel correlation.
        %   LineNoiseCriterion : If a channel has more line noise relative
        %   to its signal (default 0.85)
        %   BurstCriterion : Standard deviation cutoff for removal of
        %   bursts (via ASR) (default 4).
        %   WindowCriterion : Criterion for removing time windows that were
        %   not repaired completely (default 5).
        
        
        options = {'FlatlineCriterion',4,...
            'ChannelCriterion',0.75,...
            'ChannelCriterionMaxBadTime',0.5,...
            'LineNoiseCriterion',4,...
            'Highpass','off',...
            'BurstCriterion',20,...
            'WindowCriterion',0.1,...
            'BurstRejection','off',...
            'Distance','Euclidian',...
            'WindowCriterionTolerances', [-Inf 7]};
        
        
        fprintf('Artifact Removal options:\n')
        disp((options(:)));
        
        [EEGclean,~,BUR]  = clean_artifacts(EEG, options{:});
        
        if vis == 1
            vis_artifacts(EEGclean,EEG); end
               
        %% Find clean segments of the data
        % Segment the remaining epochs such that we get the mostn amount of data 
        
        EEG.etc.clean_sample_mask = false(1,EEGpreart.pnts*EEGpreart.trials);
        
        % run length incoding of clean and noisy segments
        [val, val_length] = my_RLE(EEGclean.etc.clean_sample_mask);
        
        % starting sample of segments
        sample = cumsum(val_length);
        
        % epoch length in samples
        ep_length = EEG.srate*epoch;
        
        % number of clean epochs in each segment
        ep_num_clean = floor(val_length.*val./ep_length);
        
        % clean idx
        clean_idx = find(ep_num_clean);
        
        % remove noisy segments of idx
        ep_num_clean = ep_num_clean(clean_idx);
             
        for i = 1:numel(clean_idx)
            
            % add sample 1 if first sample is clean 
            if clean_idx(i) == 1
                start_idx(i) = 1;
            else
                start_idx(i) = sample(clean_idx(i)-1)+1;
            end
            
            end_idx(i) = start_idx(i) + ep_length*ep_num_clean(i) - 1;
            
            
            EEG.etc.clean_sample_mask(start_idx(i):end_idx(i)) = true;
            
       end
        
        EEG.data = EEG.data(:,EEG.etc.clean_sample_mask);    
                  
        %% Epoch remaining EEG data
        EEG = eeg_regepochs( EEG, 'limits', [0 epoch], 'rmbase', NaN, 'recurrence', epoch, 'extractepochs','on');
        
        EEG.event = [];
        EEG.epoch = [];
                     
        % Channels to remove from EEG data
        bad_chan = setdiff({EEGpreart.chanlocs.labels},{EEGclean.chanlocs.labels});
        EEG.etc.clean_channel_mask = EEGclean.etc.clean_channel_mask;
        EEG = pop_select(EEG,'nochannel',bad_chan);
        
        % Check set
        EEG = eeg_checkset(EEG);
        
        if vis == 1
            vis_artifacts(EEG,EEGclean); end
        
        
        % Write into results table
        rem_epoch = [];
        num_rem_epoch = numel(rem_epoch);
              
        CleanTab.num_bad_chan(run_sub) = length(bad_chan);
        CleanTab.num_epochs(run_sub) = EEG.trials;
 
        if EEG.nbchan < 2 || EEG.trials < 2
            continue
        end

        %% Check each channel in each epoch
        rem_epoch_chan_mat = zeros(size(EEG.data,1), size(EEG.data,3));
                
        %% ICA
        
        % check set
        EEG = eeg_checkset(EEG);
        
        EEGpreica = EEG;
        
        if ica == 1            
            % Run ICA
            EEG = pop_runica(EEG, 'icatype', 'runica', 'dataset',1, 'options',{});
            
            % save the dataset without ICA
            % Avergare reference
            EEGout = pop_reref(EEG, []);
            filename = [ID '_preprocessed_noica'];
            mkdir([path_to_results '../noICA/' visit{ic} '/']);
            pop_saveset(EEGout, 'filepath', [path_to_results '../noICA/' visit{ic} '/'],'filename',filename);
            
            % ICA label
            EEG = iclabel(EEG);
      
            % Reconstruct
            % Choose components where brain has the highest probability
            Lprob = EEG.etc.ic_classification.ICLabel.classifications;
            
            % evtl end-1 to include "Others"
            [~, maxIdx] = max(Lprob(:,1:end-1), [], 2);
            ic_remove = find(maxIdx ~= 1);
            
            if numel(ic_remove)==EEG.nbchan
                continue
            end
            
            fprintf('Removing %.0f components\n',length(ic_remove));
            
            EEG = pop_subcomp(EEG, ic_remove);
            
        else
            ic_remove = [];
        end
            
        if vis == 1
            vis_artifacts(EEG,EEGpreica); end
        
        %% Interpolate all missing channels
        if size(EEG.data,1) < n_ch
            EEG = pop_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
            EEG.etc.clean_channel_mask(EEG.etc.clean_channel_mask==false)=true;
        end
        
        EEGinterpol = EEG;
        
        % Avergare reference
        EEG = pop_reref(EEG, []);
                       
        % Check set
        EEGout = eeg_checkset(EEG);
        
        % Get other info
        EEGout.visit = visit{ic};
        EEGout = get_rs_info(EEGout);
        
        %% End timer
        proc_time = toc(tstart);
        
        %% Write Cleaning Table
        CleanTab.visit(run_sub) = string(EEGout.visit);
        CleanTab.group(run_sub) = string(EEGout.group);
        CleanTab.cohort(run_sub) = string(EEGout.cohort);
        CleanTab.bad_raw_list(run_sub) = join(string(num_raw_bad_chan));
        CleanTab.bad_chan_list(run_sub) = join(string(bad_chan));
        CleanTab.all_chan_list(run_sub) = join(string(unique([num_raw_bad_chan bad_chan])));
        CleanTab.num_all_chan(run_sub) = length(unique([num_raw_bad_chan bad_chan]));
        CleanTab.num_rem_epochs(run_sub) = num_rem_epoch;
        CleanTab.num_rem_epoch_chan(run_sub) = sum(rem_epoch_chan_mat,[1 2]);
        CleanTab.num_rem_ic(run_sub) = length(ic_remove);
        CleanTab.perc_dat_int(run_sub) = CleanTab.num_all_chan(run_sub)/EEG.nbchan + CleanTab.num_rem_epoch_chan(run_sub)/EEG.nbchan/EEG.trials;
        CleanTab.process_time(run_sub) = proc_time;
        
        %% Save the preprocessed data
        filename = [ID '_preprocessed'];
        pop_saveset(EEGout, 'filepath', path_to_results,'filename',filename);
        save([path_to_results '/' filename],'EEGout');
        
        %%
        if ~ (CleanTab.num_epochs(run_sub)<15 || CleanTab.perc_dat_int(run_sub)>0.25)
            pop_saveset(EEGout, 'filepath', path_to_results_final,'filename',filename);
            save([path_to_results '/' filename],'EEGout');
        end
        
        %% End timer
        close all;
        
        save([path_to_results '../CleanTab_' strrep([visit{:}],'/','_')],'CleanTab');
        writetable(CleanTab,[path_to_results '../CleanTab_' strrep([visit{:}],'/','') '.txt']);

    end
end     
