%% Script discription

% Script is accessed via run_microstate_rs
% DISCRIPTION: Script to do microtaste segmentation using a individual approach.
% AUTHOR: Fabian Mager (fabian.martin.mager@regionh.dk)

%% Clean up
% clc, clear all, close all

%% Parameter

% analysis:
analysis = 'analysis3';

% number of microstates
%Nmicrostates = [4,7];

% load data from visit:
visit = {'visit_0','visit_6','visit_104'};

% clustering algorithm:
algorithm = 'modkmeans';

% Template for sorting
%tempfile = 'Global';


% Proprocessing folder
% preprocessing = 'preprocessed_210920_Filt_2_20';

% sorting loss measure and procedure:
label_fit = 'Competitive'; %or
% label_fit = 'Overall';

%% Adding relevant paths

% path to eeg files:
path_to_files = ['/Volumes/Promise/Shared/Microstates/RS/Preprocessed_final/' preprocessing '/'];

date_today = datestr(now,'yymmdd');

% path to result:
path_to_results = ['/Volumes/Promise/Shared/Microstates/Results_RS/' preprocessing '/' date_today '/Prototypes/' tempfile '_sorting/' label_fit '_fit/' algorithm '/'];

% figures
mkdir([path_to_results 'Figures/']);

% visit path
for i=1:numel(visit)
    mkdir([path_to_results visit{i} '/']);
end

% EEG lab
addpath('/Volumes/Promise/Shared/Microstates/eeglab/eeglab2021.1'); %path to toolbox
eeglab;

% Other
addpath([pwd '/export_fig'])

close all

%% Retrieve files
[EEG_files] = load_rs_EEG_from(visit,path_to_files);

if strcmp(tempfile,'Literature')
    % Da Cruz file   
    template_file = ['/Volumes/Promise/Shared/Microstates/DC_Proto/Proto_group_CONSCZ_cohort_PECANS1PECANS2_dc.set'];
elseif strcmp(tempfile,'Global')
    template_file = ['/Volumes/Promise/Shared/Microstates/Results_Agg/' preprocessing '/Prototypes/modkmeans/MicroGFPpeakData.set'];
end
TempEEG = pop_loadset('filename',template_file);

%% Load data into eeglab
diary([path_to_results 'diary.txt']);
date_today = datetime('now');
fprintf('\n----------------- %s -----------------\n',string(date_today))

fprintf('----- Starting with the following parameter: -----\n')
fprintf('Included subjects from: %s\n',analysis)
fprintf('Path to files: %s \n',path_to_files)
fprintf('Cohort: %s \n',visit{:})
fprintf('Clustering Algorithm: %s \n',algorithm)
fprintf('Prototype label measure: %s \n','Spatial Correlation')
fprintf('Path to results: %s \n',path_to_results)

diary off

ALLEEG = [];

for i = 1:numel(EEG_files)
    
    % Printing processs
    fprintf('Loading subject %.0f of %.0f\n',i,numel(EEG_files))
    
    EEG = pop_loadset('filename',EEG_files(i).name,'filepath',EEG_files(i).folder);
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, i);
    
    fprintf('Subject %s included in analysis.\n', EEG.subject)
    
end

diary on

% Print subject informations
N = numel(ALLEEG);
fprintf('I have loaded %.0f subjects now.\n',N)


%% Segmentation

% Select data for microstate analysis
Npeaks = 1000;
MinPeakDist = 20;
avgref = 0;
gfp_normalise = 0;
GFPthresh = 10;

% Select clustering parameter
sorting = 'Global explained variance';
verbose = 1;
clust_normalise = 0;
Nrepetitions = 100;
max_iterations = 1000;
threshold = 1e-08;
fitmeas = 'GEV';
optimised = 1;


% results table
optClustTab = table();

% Printing status
fprintf('------- Aggregating data using: --------\n')
fprintf('Number of peaks per subject: %.0f\n',Npeaks)
fprintf('Minimum peak distance: %.0f\n',MinPeakDist)
fprintf('Average referencing: %.0f\n',avgref)
fprintf('Normalizing GFP peaks: %.0f\n',gfp_normalise)
fprintf('Excluding GFP peaks higher than GFP std: %.0f\n',GFPthresh)
fprintf('------- Clustering using: --------\n')
fprintf('Sorting: %.s\n',sorting)
fprintf('Number of microstates: %.0f\n',Nmicrostates)
fprintf('Normalizing EEG: %.0f\n',clust_normalise)
fprintf('Nrepetitions: %.0f\n',Nrepetitions)
fprintf('Fit measure: %s\n',fitmeas)
fprintf('Optimised: %.0f\n',optimised)

diary off

for i=1:N
    
    
    %% Load data into eeglab
    fprintf('Processing EEG file %.0f of %.0f\n', i, N);

    %% Select EEG data set
    [ALLEEG, EEG, CURRENTSET] = pop_newset( ALLEEG, EEG, CURRENTSET,'retrieve', i);
        
    %% Select data for microstate analysis
    [EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG,...
        'datatype', 'spontaneous',...
        'avgref', avgref, ...
        'normalise', gfp_normalise, ...
        'Npeaks',Npeaks,...
        'MinPeakDist', MinPeakDist, ...
        'GFPthresh', GFPthresh);

    %% Store data in a new EEG structure
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, i,'overwrite','on');

    %% Perform the microstate segmentation
    fprintf('----- Starting segmentation -----\n')
     
    switch algorithm
        case 'modkmeans'
            EEG = pop_micro_segment(EEG,...
                'algorithm', algorithm, ...
                'sorting', sorting, ...
                'Nmicrostates', Nmicrostates, ...
                'verbose', verbose, ...
                'normalise', clust_normalise, ...
                'Nrepetitions', Nrepetitions, ...
                'max_iterations', max_iterations, ...
                'threshold', threshold, ...
                'fitmeas', fitmeas,...
                'optimised', optimised);
            
        case 'aahc'
            EEG = pop_micro_segment(EEG,...
                'algorithm', algorithm,...
                'sorting', sorting, ...
                'Nmicrostates', Nmicrostates, ...
                'verbose', verbose, ...
                'normalise', clust_normalise);
    end
    
    %% Clustering Results
    optClustTab.subject(i) = string(EEG.subject);
    [~, optIdx] = (min(EEG.microstate.Res.CV));
    optClustTab.CVopt(i) = Nmicrostates(optIdx);
    [~, optIdx] = (min(EEG.microstate.Res.W));
    optClustTab.Wopt(i) = Nmicrostates(optIdx);
    [~, optIdx] = (min(EEG.microstate.Res.KL));
    optClustTab.KLopt(i) = Nmicrostates(optIdx);
    
    optClustTab.GEVoptCV(i,:) = EEG.microstate.Res.GEV(optIdx);
    
    save([path_to_results 'optClustTab'],'optClustTab');
    writetable(optClustTab,[path_to_results 'optClustTab.txt']);
    
    EEG.microstate.Res_unsorted = EEG.microstate.Res;

    %% Sort Microstates
    for k = Nmicrostates
        % Select active number of microstates
        EEG = pop_micro_selectNmicro(EEG, 'Nmicro', k);
        if strcmp(tempfile,'Global')
            TempEEG = pop_micro_selectNmicro(TempEEG, 'Nmicro', k);
        end
        % Rearanging microstates
        EEG = sort_microstates_spCorr(EEG,TempEEG,label_fit);
    end
    
    %% Store data in a new EEG structure
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, i,'overwrite','on');

    %% Save EEG
    pop_saveset(EEG, 'filepath', [path_to_results EEG.visit '/'], 'filename', [EEG.subject '_prototyped']);
    
    %% plot microstate prototype topographies
    MicroPlotTopo(EEG, 'plot_range', [] );
    fig = gcf;
    set(fig,'visible','off','WindowState','maximized','color','w')
    export_fig(fig, [path_to_results 'Figures/' EEG.subject '_protorange' '.fig']);
    export_fig(fig, [path_to_results 'Figures/' EEG.subject '_protorange' '.png']);
    close all
    
    %% Spatial correlation

    for k = Nmicrostates
        
        % Select active number of microstates
        EEG = pop_micro_selectNmicro(EEG, 'Nmicro', k);
        
        k_idx = find(EEG.microstate.algorithm_settings.Nmicrostates==k);
        
        spat_corr_diag{i,k_idx} = diag(EEG.microstate.Res.Spat_corr{k_idx});
        
    end
    
    close all
end

%% Plot Spatial Correlations of Microstates
fig = figure();

for k = 1:numel(Nmicrostates)
    spat_corr_mat = [spat_corr_diag{:,k}]';
    spat_corr_m(k,:) = mean(spat_corr_mat,1,'omitnan');
    spat_corr_se(k,:) = std(spat_corr_mat,0,1,'omitnan')./sqrt(N);
    
    subplot(1,numel(Nmicrostates),k)
    boxplot(spat_corr_mat,repmat(string({ALLEEG.group}),1,4),'Plotstyle','compact')
    ylabel('Spatial Correlation')
    xlabel('Group')
    xtickangle(45)
        
    spat_corr_tab = table(string({ALLEEG.subject})',spat_corr_mat);
    writetable(spat_corr_tab,[path_to_results 'spat_corr_tab_' num2str(Nmicrostates(k)) '.txt'])
end

set(fig,'color','w')
saveas(fig, [path_to_results 'Figures/Spat_corr_boxplot.fig']);
saveas(fig, [path_to_results 'Figures/Spat_corr_boxplot.png']);


fig = figure();
plot(Nmicrostates,spat_corr_m,'LineWidth',2)
hold on
grid on
errorbar(repmat(Nmicrostates,4,1)',spat_corr_m,spat_corr_se,'ok')
legend(['A':'D']','Location','southeast')
ylim([0.5, 1])
xlim([Nmicrostates(1)-0.5, Nmicrostates(end)+0.5])
ylabel('Average Spatial Correlation and SE for each MS')
xlabel('Number of Microstates')


set(fig,'color','w')
saveas(fig, [path_to_results 'Figures/Spat_corr_over_each_MS.fig']);
saveas(fig, [path_to_results 'Figures/Spat_corr_over_each_MS.png']);

%% Plot mean microstates
A_temp = TempEEG.microstate.prototypes;
for k = Nmicrostates
    A = zeros(64,k,N);
    for i = 1:N
        EEG = pop_micro_selectNmicro(EEG, 'Nmicro', k);
        A(:,:,i) = EEG.microstate.prototypes;
    end
    A_mean=mean(A,3);
    fig = figure();
    for t = 1:4
        subplot(2,k,t)
        topoplot(A_temp(:,t),EEG.chanlocs,'style','map','electrodes','off');
    end
    for kk = 1:k
        subplot(2,k,k+kk)
        topoplot(A(:,kk),EEG.chanlocs,'style','map','electrodes','off');
    end
    set(fig,'color','w')
    saveas(fig, [path_to_results 'Figures/Mean_MS' num2str(k) '.fig']);
    saveas(fig, [path_to_results 'Figures/Mean_MS' num2str(k) '.png']);
end
%%
close all
fig = figure();
N_row = 1+numel(Nmicrostates);
N_col = Nmicrostates(end);
for t = 1:4
    subplot(N_row,N_col,t)
    topoplot(A_temp(:,t),EEG.chanlocs,'style','map','electrodes','off');
end
n_row = 2;
for k = Nmicrostates
    A = zeros(64,k,N);
    for i = 1:N
        EEG = pop_micro_selectNmicro(EEG, 'Nmicro', k);
        A(:,:,i) = EEG.microstate.prototypes;
    end
    A_mean=mean(A,3);
    
    for kk = 1:k
        subplot(N_row,N_col,((n_row-1)*N_col)+kk)
        topoplot(A(:,kk),EEG.chanlocs,'style','map','electrodes','off');
    end
    n_row = n_row + 1;
end
set(fig,'color','w')
saveas(fig, [path_to_results 'Figures/Mean_MS_all.fig']);
saveas(fig, [path_to_results 'Figures/Mean_MS_all.png']);
        


%% Backfitting

minFrame = [1];
sWdth = [3];
sWght = [10];

label_type = 'backfit';
% smooth_type = 'reject segments';
smooth_type = 'windowed';

% path to result:
%date_today = datestr(now,'yymmdd');
date_today = '220812';
path_to_files = ['/Volumes/Promise/Shared/Microstates/Results_RS/' preprocessing '/' date_today '/Prototypes/' tempfile '_sorting/' label_fit '_fit/' algorithm '/'];

EEG_files = [];

for i = 1:length(visit)    
    EEG_files = [EEG_files;dir([path_to_files visit{i} '/*.set'])];  
end

idx = 0;

ALLEEG = [];

for i = 1:numel(EEG_files)
    
    % Printing processs
    fprintf('Loading subject %.0f of %.0f\n',i,numel(EEG_files))
    
    EEG = pop_loadset('filename',EEG_files(i).name,'filepath',EEG_files(i).folder);
    
    EEG.group = [EEG.type '_' extractAfter(EEG.visit,'_')];
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, i);
    
    fprintf('Subject %s included in analysis.\n', EEG.subject)
    
end

% Print subject informations
N = numel(ALLEEG);
fprintf('I have loaded %.0f subjects now.\n',N)


for k = Nmicrostates
    
    % path to results:
    path_to_results = ['/Volumes/Promise/Shared/Microstates/Results_RS/' preprocessing '/' date_today '/Backfitted/' tempfile '_sorting/' label_fit '_fit/' algorithm '/' num2str(k) '_micro/'];
    
    % Setting up results folder
    mkdir([path_to_results 'Statistics/Scatter/']);
    for i=1:numel(visit)
        mkdir([path_to_results visit{i} '/']);
    end
    
    %  Prepare results table
    ResTab = setup_res_tab(N,k);

    
    for j = 1:length(sWdth)
        
        diary on
        
        % Print status
        fprintf('----- Starting Backfitting with the following parameter: -----\n')
        fprintf('Number of microstates: %.0f\n',k)
        fprintf('Labeltype: %s\n',label_type)
        fprintf('Minimum time of microstate: %.0f\n',minFrame(j))
        fprintf('Smoothing width: %.0f\n',sWdth(j))
        fprintf('Smoothing weight: %.0f\n',sWght(j))
        
        diary off
        
        filename = ['statistics_minTime' num2str(minFrame(j))...
                '_width' num2str(sWdth(j))...
                '_weight' num2str(sWght(j))];
        
        for i=1:N
            %% Load data into eeglab
            fprintf('Processing EEG file %.0f of %.0f\n', i, N);
            
            % Select EEG data set
            [ALLEEG, EEG, CURRENTSET] = pop_newset( ALLEEG, EEG, CURRENTSET,'retrieve', i);
            % Select active number of microstates
            EEG = pop_micro_selectNmicro(EEG, 'Nmicro', k);
            
            % back-fit microstates on EEG, do not account for polarity
            fprintf('I am backfitting prototypes from %s to %s...',ALLEEG(i).filename(1:4), EEG.filename(1:4))
            EEG = pop_micro_fit( EEG, 'polarity', 0 );
            
            % temporally smooth microstates labels
            fprintf('smoothing...')
            
            if strcmp(smooth_type,'windowed')
            % windowed smoothing
            EEG = pop_micro_smooth( EEG,...
                'label_type', label_type, ...
                'smooth_type', smooth_type, ...
                'smooth_width', sWdth(j),...
                'smooth_weight',sWght(j));
            
            elseif strcmp(smooth_type,'reject segments')
            % reject segments
            EEG = pop_micro_smooth( EEG,...
                'label_type', label_type, ...
                'smooth_type', smooth_type, ...
                'minTime', minFrame(j)/EEG.srate*1000 );
            
            else
                error('check label type')
            end
            
            EEG = my_micro_stats(EEG,minFrame(j));
            
            fprintf('done!\n')
            
            %% Store data in a new EEG structure
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, i,'overwrite','on');

            %% Save results to table
            ResTab = write_to_res_tab(ResTab,EEG,i);
            
            %% Save EEG
            % Saving epoched data
            full_path_to_results = [path_to_results EEG.visit...
                '/' EEG.visit '/'...
                'smoothing_minFrame' num2str(minFrame(j))...
                '_width' num2str(sWdth(j))...
                '_weight' num2str(sWght(j)) '/'];
            mkdir(full_path_to_results);
            pop_saveset(EEG, 'filename', [full_path_to_results strrep(EEG.filename,'preprocessed','backfitted')]);
            
            %% Save table
            writetable(ResTab,[path_to_results 'Statistics/' filename '.txt'])
        end
        
        %% Plot Results
        
        % find Group idx
        gIdx = find(string(ResTab.Properties.VariableNames)=='Group');

        % plot figure
        settings.showplot = 'on';
        settings.color = 'group';

%        [settings.h, settings.p] = run_ttest(ResTab,0.05);

        fig1 = plotScatter_rs(ResTab,ResTab{:,gIdx},settings);

        set(fig1,'Units','Normalized','Outerposition',[0 0 1 1]);
        fprintf('saving...')
        saveas(fig1,[path_to_results 'Statistics/Scatter/' filename '_scat.png']);
        savefig(fig1,[path_to_results 'Statistics/Scatter/' filename '_scat.fig']);
    end
    
    close all
end