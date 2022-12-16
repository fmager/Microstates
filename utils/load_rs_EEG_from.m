function [EEG_files] = load_rs_EEG_from(cohort,path)
%LOAD_EEG Function, that loads the data into EEG_files structure.
%   COHORT = (string) cell array, from which cohort group to load.
%   GROUP = (integer) array indicating group, where 0 = conotrol, 1 =
%   schizophrenic. Example group = [0, 1] loads control and schizophrenic

%% Retrieve files
EEG_files = [];
    for i = 1:length(cohort)
        eeg_dir = fullfile([path cohort{i} '/']);
        EEG_files = [EEG_files; dir([eeg_dir,'*.set'])];
    end

end

