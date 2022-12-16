function EEGout = get_rs_info(EEG)

    EEGout = EEG;
    
    if strcmp(EEGout.subject(1),'D')
        EEGout.cohort = 'pecans';        
    elseif strcmp(EEGout.subject(1),'E')
        EEGout.cohort = 'pecans2';
    elseif strcmp(EEGout.subject(1),'P')
        EEGout.cohort = 'prodromal';
    else
        EEGout.cohort = '';
    end
    
    if strcmp(EEGout.cohort,'pecans') || strcmp(EEGout.cohort,'pecans2')
        if strcmp(EEGout.subject(2),'1')
            EEGout.type = 'patient';        
        elseif strcmp(EEGout.subject(2),'2') || strcmp(EEGout.subject(2),'9')
            EEGout.type = 'control';
        end
    elseif strcmp(EEGout.cohort,'prodromal')
        EEGout.type = 'uhr';
    else
        EEGout.type = '';
    end
    
    EEGout.group = [EEGout.type '_' extractAfter(EEGout.visit,'_')];

end