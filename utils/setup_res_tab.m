function [ResTab] = setup_res_tab(N,num_micro)

letters = 'A':'Z';

for nm = 1:num_micro
    dur{nm} = ['Dur' letters(nm)];
    occ{nm}= ['Occ' letters(nm)];
    cov{nm} = ['Cov' letters(nm)];
    gev{nm} = ['Gev' letters(nm)];
    gfp{nm} = ['Gfp' letters(nm)];
    corr{nm} = ['Corr' letters(nm)];
end

varNames = {'Subject','Session','Group','Gender', 'Num_micro','GEVtotal',...
            dur{:}, occ{:}, cov{:}, gev{:}, gfp{:}, corr{:} };

temp_var_types = repmat({'double'},1,6*num_micro);

varTypes = {'string','string','string','string','double','double',...
    temp_var_types{:}};

varNum = numel(varNames);

ResTab = table('Size',[N varNum],'VariableTypes',varTypes,'VariableNames',varNames);


end

