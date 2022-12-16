function fig = plotScatter_rs(T,Grouping,settings)
%PLOTPARAMETER Plots the statistics in table T as a scatter plot

showplot = settings.showplot;
groupstr = unique(string(Grouping));
% setup parameter
N = height(T);
letters = 'A':'Z';
Nm = T.Num_micro(1); % number of microstates
Ng = length(unique(Grouping));
%groupstr(1:Ng) = groupstr([2,6,4,5,7,1,3]);

Dur = [];
Occ = [];
Cov = [];
Gev = [];
Corr = [];

for nm = 1:Nm
    % extract duration
    Dur = [Dur, T.(['Dur' letters(nm)])];    
    % extract occurence
    Occ = [Occ, T.(['Occ' letters(nm)])];
    % extract coverage
    Cov = [Cov, T.(['Cov' letters(nm)])];
    % extract gev
    Gev = [Gev, T.(['Gev' letters(nm)])];
    % extract correlation
    %Corr = [Corr, T.(['Corr' letters(nm)])];
end

Cov = Cov.*100;
Gev = Gev.*100;

% setup data matrix
Xall = {Dur,Occ,Cov,Gev};

% ylabel
yfeat = {'Duration in ms','Occurence in 1/s','Coverage in \%','GEV in \%'};
ylimits = {[0 200],[0 7],[0 60],[0 50]};

% setup group matrix
%G = repmat(string(Grouping),1,Nm);

% setup axis
offset = 0:1:(Nm-1); % offset
for i = 1:Ng
    x_grid(i,:) = i:Ng:(Ng*Nm);
end
x_grid = x_grid + offset;
%x = [1:Ng:(Ng*Nm-1);2:Ng:(Ng*Nm);3:Ng:(Ng*Nm)]+s;
for i = 1:Ng
    X_scat{i} = repmat(linspace(-.3,.3,sum(Grouping==groupstr{i}))',1,Nm)+x_grid(i,:);
end
% setup x indices for each group
%xCon = repmat(linspace(-.3,.3,sum(Grouping==groupstr{1}))',1,Nm)+x(1,:);
%xScz = repmat(linspace(-.3,.3,sum(Grouping==groupstr{2}))',1,Nm)+x(2,:);

% setup figure
fig = figure('visible',showplot);
set(fig,'color','w');

% setup color
if strcmp(settings.color,'group')
    colormap lines;
    cmap = colormap;
    cmap(3,:) = [];
    for i = 1:Ng
        for j = 1:Nm
            color{i,j} = cmap(i,:);
        end
    end
elseif strcmp(settings.color,'corr')
    colormap parula
    for i = 1:Ng
        for j = 1:Nm
            %color{i,j} = Corr(Grouping==groupstr{i},j);
        end
    end
end


%%%%%%%% HARDCODING OF DA CRUZ RESULTS %%%%%%%%%
DC_con = zeros(numel(Xall),Nm);
DC_scz = zeros(numel(Xall),Nm);

DC_con(:,1:4) = [72 71 90 82 ; 1.8 1.9 2.5 2.3; 18.8 19.3 34.0 27.9; 0 0 0 0];
DC_scz(:,1:4) = [69 66 107 70; 1.7 1.7 2.8 1.6; 18.2 16.2 46.4 19.3; 0 0 0 0];

% DC_con(2,:) = DC_con(3,:)./100.*1000./DC_con(1,:);
% DC_scz(2,:) = DC_scz(3,:)./100.*1000./DC_scz(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% HARDCODING OF TOMESCU RESULTS %%%%%%%%%
% TOM_con = [44 47 60 51; 4.2 4.7 5.2 4.5; 19 22.5 33 25.4; 0.07 0.1 0.19 0.14];
% TOM_scz = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% HARDCODING OF De Bock RESULTS %%%%%%%%%
% DB_con = [66.7 65.1 77.3 66.7; 3.47 3.31 4.33 3.69; 22.74 20.85 32.23 24.18; 0 0 0 0];
% DB_scz = [73.7 59.6 70.8 60.2; 4.3 3 4.21 3.7; 30.8 18 29.2 22; 0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% HARDCODING OF STOFFER %%%%%%%%%%%%%%%%%
% sub_stoffer = ['D150';'D154';'D156';'E104'];
% com_stoffer = intersect(sub_stoffer,T.Subject);
% for i = 1:length(com_stoffer)
%     x_stoffer(i) = find(strcmp(com_stoffer(i),string(T.Subject(Grouping==groupstr{2}))));
%     idx_stoffer(i) = find(strcmp(com_stoffer(i),string(T.Subject)));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:length(Xall)
    % setup subplot
    subplot(1,length(Xall),d)
    hold on
    
    % extract data
    X = Xall{d};
    
    % extract mean and std
    %mX   = [mean(X(Grouping==groupstr{1},:)); mean(X(Grouping==groupstr{2},:))];
    %stdX = [std(X(Grouping==groupstr{1},:)); std(X(Grouping==groupstr{2},:))];   
    
    % plotting
    for i = 1:size(X,2) %for each microstate
        
        for j = 1:Ng %for each group
            
            scat = scatter(X_scat{j}(:,i),X(Grouping==groupstr{j},i),35,color{j,i},'.');
            %scat.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Sub',T.Subject(Grouping==groupstr{j}));
            
            % mean control
            plot(x_grid(j,i),mean(X(Grouping==groupstr{j},i)),'dk','LineWidth',2)
            % SE control
            errorbar(x_grid(j,i),mean(X(Grouping==groupstr{j},i)),std(X(Grouping==groupstr{j},i))/sqrt(sum(Grouping==groupstr{j})),'k','LineWidth',1)
            
            if strcmp(groupstr{j},'FUS')
                plot(X_scat{j}(:,i),DC_scz(d,i),'*g','LineWidth',2','HandleVisibility','off')
            elseif strcmp(groupstr{j},'CON')
                plot(X_scat{j}(:,i),DC_con(d,i),'*g','LineWidth',2','HandleVisibility','off')
            end
            
        end
        
        % plot control
        %sCon = scatter(xCon(:,i),X(Grouping==groupstr{1},i),35,color{1,i},'.');
        %sCon.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Sub',T.Subject(Grouping==groupstr{1}));
        
        
        % plot schizophrenic
        %sScz = scatter(xScz(:,i),X(Grouping==groupstr{2},i),35,color{2,i},'.');
        %sScz.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Sub',T.Subject(Grouping==groupstr{2}));
        
        % mean control
        %plot(X_scat(1,i),mX(1,i),'dk','LineWidth',2)
        % SE control
        %errorbar(X_scat(1,i),mX(1,i),stdX(1,i)/sqrt(Ncon),'k','LineWidth',1)
        
        % statistics
        %if isfield(settings,'p')
         %   text(sum(X_scat(:,i))/2, ylimits{d}(2)*0.05, sprintf('p=%.3f',settings.p(d,i)),'FontSize',9,'HorizontalAlignment','center')
        %end
        
        %%%% CON DA CRUZ %%%%
        %plot(X_scat(1,i),DC_con(d,i),'*g','LineWidth',2')
        %%%%%%%%%%%%%%%%%

        %%%% Stoffer %%%%
        %scatter(xScz(x_stoffer,i),X(idx_stoffer,i),35,color{2,i},'ob');
        %%%%%%%%%%%%%%%%%

        % mean schizophrenic
        %plot(X_scat(2,i),mX(2,i),'dk','LineWidth',2)
        % std schizophrenix
        %errorbar(X_scat(2,i),mX(2,i),stdX(2,i)/sqrt(Nscz),'k','LineWidth',1)

        %%%% SCZ DA CRUZ %%%%
        %plot(X_scat(2,i),DC_scz(d,i),'*g','LineWidth',2')
        %%%%%%%%%%%%%%%%%
         
        %%%% DB %%%%
%         plot(x(1,i),DB_con(d,i),'*m','LineWidth',2')
%         plot(x(2,i),DB_scz(d,i),'*m','LineWidth',2')
        %%%%%%%%%%%%%%%%%%         
    end
    
    % label
    xticks(median(x_grid,1))
    xticklabels({letters(1:Nm)'})
    xlabel('Microstate')
    ylabel(yfeat(d))
    ylim(ylimits{d})
        
    
end
    % legend
    lgnd = [groupstr'; repmat("mean",[1, Ng]); repmat("se",[1, Ng])];
    lgnd_pos = [0.37 .7 0.1 0.2];
    set(legend(lgnd(:)),'Position',lgnd_pos,'Units','Normalized','Orientation','Horizontal','NumColumns',3,'AutoUpdate','off')
        
end

