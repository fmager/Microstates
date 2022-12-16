function [EEG] = my_micro_stats(EEG,n_single)

% settings
settings.polarity = 0;
% settings.label_type = label_type;

Ntrials = size(EEG.data,3);
K = size(EEG.microstate.prototypes,2);


GEVtotal = zeros(Ntrials,1);
MGFP = nan(Ntrials,K);
MDur = zeros(Ntrials,K);
MOcc = zeros(Ntrials,K);
TCov = zeros(Ntrials,K);
GEV = nan(Ntrials,K);
MspatCorr = nan(Ntrials,K);
TP = zeros(K,K,Ntrials);
Nsingletons = zeros(Ntrials,1);

for trial=1:Ntrials
    
    % truncate first and last microstate
    [data, labels] = remove_first_last(EEG.data(:,:,trial),EEG.microstate.fit.labels(trial,:));
    if isempty(data) || isempty(labels)
        warning('Ignoring one epoch.')
        continue
    end
       
    % calculate statistics
    Mstats = MicroStats_ignoreSingletons(data,EEG.microstate.prototypes,labels,n_single,settings.polarity,EEG.srate);

    GEVtotal(trial) = Mstats.GEVtotal;
    MGFP(trial,1:K) = Mstats.Gfp;
    MOcc(trial,1:K) = Mstats.Occurence;
    MDur(trial,1:K) = Mstats.Duration;
    TCov(trial,1:K) = Mstats.Coverage;
    GEV(trial,1:K) = Mstats.GEV;
    MspatCorr(trial,1:K) = Mstats.MspatCorr;
    TP(:,:,trial) = Mstats.TP;
    Nsingletons(trial) = Mstats.Nsingletons;
        
end

Mstats.polarity = settings.polarity;

if Ntrials > 1
    %sum
    Mstats.sums.Nsingletons = sum(Nsingletons);
    % mean parameters
    Mstats.avgs.GEVtotal = nanmean(GEVtotal,1);
    Mstats.avgs.Gfp = nanmean(MGFP,1);
    Mstats.avgs.Occurence = nanmean(MOcc,1);
    Mstats.avgs.Duration = nanmean(MDur,1);
    Mstats.avgs.Coverage = nanmean(TCov,1);
    Mstats.avgs.GEV = nanmean(GEV,1);
    Mstats.avgs.MspatCorr = nanmean(MspatCorr,1);
    Mstats.avgs.TP = nanmean(TP,3);
    
    % standard deviation of parameters
    Mstats.avgs.stdGEVtotal = nanstd(GEVtotal,1);
    Mstats.avgs.stdGfp = nanstd(MGFP,1);
    Mstats.avgs.stdOccurence = nanstd(MOcc,1);
    Mstats.avgs.stdDuration = nanstd(MDur,1);
    Mstats.avgs.stdCoverage = nanstd(TCov,1);
    Mstats.avgs.stdGEV = nanstd(GEV,1);
    Mstats.avgs.stdMspatCorr = nanstd(MspatCorr,1);
    Mstats.avgs.stdTP = nanstd(TP,[],3);
    
    EEG.microstate.stats.avgs = Mstats.avgs;
    EEG.microstate.stats.sums = Mstats.sums;
else
    EEG.microstate.stats = Mstats;
end

end

