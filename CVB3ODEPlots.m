function [] = CVB3ODEPlots(MedianTable,UpperQuantileTable,LowerQuantileTable,MeanTable,varargin)
%function [] = CVB3ODEPlots(MedianTable,UpperQuantileTable,LowerQuantileTable,varargin)
%Function: Plots specified figures
%
%INPUTS:
    %MedianTable: Median results table of CVB3 infection simulation
    %UpperQuantileTable: Upper quantile results table of CVB3 infection simulation
    %LowerTable: Lower quantile results table of CVB3 infection simulation
    %MeanTable: Mean results table of CVB3 infection simulation
    %PlotType: linear or semilog specification for the plots
    %PlotOutputs: vector of specific plots desired
    
%% Parse Inputs

% Set defaults
defaultPlotType = 'semilog';
defaultPlotOutputs = 'All';
defaultPopulationSetting = 'Single Cell';

CVB3ODEPlotsparser = inputParser;

% Set general criteria for inputs
validInput = @(x) istable(x);
validPlotType = @(x) strcmpi(x,'linear') || strcmpi(x,'semilog');
validPlotOutputs = @(x) iscell(x);
validPopulationString = @(x) strcmpi(x,'Single Cell') || strcmpi(x,'Population');

% Set parser options and valid input criteria
addRequired(CVB3ODEPlotsparser,'MedianTable',validInput);
addRequired(CVB3ODEPlotsparser,'UpperQuantileTable',validInput);
addRequired(CVB3ODEPlotsparser,'LowerQuantileTable',validInput);
addRequired(CVB3ODEPlotsparser,'MeanTable',validInput);
addOptional(CVB3ODEPlotsparser,'PlotType',defaultPlotType,validPlotType);
addOptional(CVB3ODEPlotsparser,'PlotOutputs',defaultPlotOutputs,validPlotOutputs);
addOptional(CVB3ODEPlotsparser,'PopulationSetting',defaultPopulationSetting,validPopulationString);

% Extract parsed inputs
parse(CVB3ODEPlotsparser,MedianTable,UpperQuantileTable,LowerQuantileTable,MeanTable,varargin{:});
MedianTable = CVB3ODEPlotsparser.Results.MedianTable;
UpperQuantileTable = CVB3ODEPlotsparser.Results.UpperQuantileTable;
LowerQuantileTable = CVB3ODEPlotsparser.Results.LowerQuantileTable;
MeanTable = CVB3ODEPlotsparser.Results.MeanTable;
PlotType = CVB3ODEPlotsparser.Results.PlotType;
PlotOutputs = CVB3ODEPlotsparser.Results.PlotOutputs;
PopulationSetting = CVB3ODEPlotsparser.Results.PopulationSetting;

%% Creates function for more concise plot settings

function [] = font_ax(t,xlab,ylab,fsize_ax,fweight,bwidth,cbar_flag)
    %Function: Sets features of the plot display
    %
    %INPUTS:
        %t: Plot title
        %xlab: X-axis label
        %ylab: Y-axis label
        %fsize_ax: Sets the figure's font size
        %fweight: Sets font to bold if desired
        %bwidth: Thickness of the lines in the display
        %cbar_flag: 1 adds a color bar, 0 does not

    title(t,'FontSize',fsize_ax+5,'FontWeight',fweight);
    % title(['Selected Slice (Depth ',num2str(d_ind),'mm)'],'FontSize',15,'FontWeight','bold');
    % xh=xlab;
    % yh=ylab;
    xh=xlabel(xlab);
    yh=ylabel(ylab);
    set([xh,yh],'FontSize',fsize_ax,'FontWeight','bold');
    set(gca,'FontSize',fsize_ax,'FontWeight','bold','LineWidth',bwidth);
    if(cbar_flag==1)
        h=colorbar;
        set(h,'FontSize',fsize_ax,'FontWeight','bold','LineWidth',3);
    end
end

%% Plotting

%Plotting and caption variables
Time = table2array(MedianTable(:,{'Time'}));
MOI = cell2mat(extractBetween(MedianTable.Properties.Description,'MOI: ',','));
VarNames = MedianTable.Properties.VariableNames;

if strcmpi(PlotOutputs,'All')
    %Default Plots
    
    %Delivery
    figure
    sgtitle(sprintf('CVB3 Receptor Binding and Intracellular Delivery at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    for i = 2:16
        subplot(3,5,i-1)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
    %Translation and Cytoplasmic Species
    figure
    sgtitle(sprintf('Viral Translation and Cytoplasmic Species at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    for i = 17:24
        subplot(2,4,i-16)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
    %Replication and VRO Species
    figure
    sgtitle(sprintf('Viral Replication Organelle Species at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    for i = 25:31
        subplot(2,4,i-24)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
    %Virion Assembly
    figure
    sgtitle(sprintf('Virion Capsid Assembly Species at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    for i = 32:42
        subplot(2,11,i-31)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    for i = 44:53
        subplot(2,11,i-31)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
    %Empty vs Filled Virions
    figure
    MedianFilledOutput = table2array(MedianTable(:,{'Virions'}));
    MedianEmptyOutput = table2array(MedianTable(:,{'Empty Provirions'}));
    UQFilledOutput = table2array(UpperQuantileTable(:,{'Virions'}));
    UQEmptyOutput = table2array(UpperQuantileTable(:,{'Empty Provirions'}));
    LQFilledOutput = table2array(LowerQuantileTable(:,{'Virions'}));
    LQEmptyOutput = table2array(LowerQuantileTable(:,{'Empty Provirions'}));
    MeanFilledOutput = table2array(MeanTable(:,{'Virions'}));
    MeanEmptyOutput = table2array(MeanTable(:,{'Empty Provirions'}));
    
    if strcmpi(PlotType,'linear') == 1
        hold on
        filledPlot = plot(Time,MedianFilledOutput,'r-','LineWidth',2);
        emptyPlot = plot(Time,MedianEmptyOutput,'b-','LineWidth',2);
        
        if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
            plot(Time,MeanFilledOutput,'r:','MarkerSize',5)
            plot(Time,MeanEmptyOutput,'b:','MarkerSize',5)
        elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
            plot(Time,UQFilledOutput,'r--','LineWidth',0.75)
            plot(Time,LQFilledOutput,'r--','LineWidth',0.75)
            
            plot(Time,UQEmptyOutput,'b--','LineWidth',0.75)
            plot(Time,LQEmptyOutput,'b--','LineWidth',0.75)
        end
        
        legend([filledPlot,emptyPlot],'Filled Virions','Empty Virions')
        hold off
    elseif strcmpi(PlotType,'semilog') == 1
        filledPlot = semilogy(Time,max(MedianFilledOutput,0.000001),'r-','LineWidth',2);
        hold on
        emptyPlot = semilogy(Time,max(MedianEmptyOutput,0.000001),'b-','LineWidth',2);
        
        if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
            semilogy(Time,max(MeanFilledOutput,0.000001),'r:','MarkerSize',5)
            semilogy(Time,max(MeanEmptyOutput,0.000001),'b:','MarkerSize',5)
        elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
            semilogy(Time,max(UQFilledOutput,0.000001),'r--','LineWidth',0.75)
            semilogy(Time,max(LQFilledOutput,0.000001),'r--','LineWidth',0.75)
            
            semilogy(Time,max(UQEmptyOutput,0.000001),'b--','LineWidth',0.75)
            semilogy(Time,max(LQEmptyOutput,0.000001),'b--','LineWidth',0.75)
        end
        
        legend([filledPlot,emptyPlot],'Filled Virions','Empty Virions')
        hold off
    end
    font_ax(sprintf('Infectious Viral Progeny vs Empty Provions at an MOI of %s',MOI),'Time (h)','',10,'bold',2.0,0)
    
    %Host-viral interaction
    figure
    sgtitle(sprintf('Host Anti-viral Response at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    for i = 55:56
        subplot(1,2,i-54)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
    %Summary species
    figure
    sgtitle(sprintf('Summary Species at an MOI of %s',MOI),'FontSize',12,'FontWeight','bold')
    subplot(1,3,1)
    MedianTotPosOutput = table2array(MedianTable(:,{'Total +ssRNA'}));
    MedianTotNegOutput = table2array(MedianTable(:,{'Total -ssRNA'}));
    
    UQTotPosOutput = table2array(UpperQuantileTable(:,{'Total +ssRNA'}));
    UQTotNegOutput = table2array(UpperQuantileTable(:,{'Total -ssRNA'}));
    
    LQTotPosOutput = table2array(LowerQuantileTable(:,{'Total +ssRNA'}));
    LQTotNegOutput = table2array(LowerQuantileTable(:,{'Total -ssRNA'}));
    
    MeanTotPosOutput = table2array(MeanTable(:,{'Total +ssRNA'}));
    MeanTotNegOutput = table2array(MeanTable(:,{'Total -ssRNA'}));
    
    if strcmpi(PlotType,'linear') == 1
        hold on
        TotPosPlot = plot(Time,MedianTotPosOutput,'r-','LineWidth',2);
        TotNegPlot = plot(Time,MedianTotNegOutput,'b-','LineWidth',2);
        
        if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
            plot(Time,MeanTotPosOutput,'r:','MarkerSize',5)
            plot(Time,MeanTotNegOutput,'b:','MarkerSize',5)
        elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
            plot(Time,UQTotPosOutput,'r--','LineWidth',0.75)
            plot(Time,LQTotPosOutput,'r--','LineWidth',0.75)
            
            plot(Time,UQTotNegOutput,'b--','LineWidth',0.75)
            plot(Time,LQTotNegOutput,'b--','LineWidth',0.75)
        end
        
        legend([TotPosPlot,TotNegPlot],'+ssRNA','-ssRNA')
        hold off
    elseif strcmpi(PlotType,'semilog') == 1
        TotPosPlot=semilogy(Time,max(MedianTotPosOutput,0.000001),'r-','LineWidth',2);
        hold on
        TotNegPlot=semilogy(Time,max(MedianTotNegOutput,0.000001),'b-','LineWidth',2);
        
        if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
            semilogy(Time,max(MeanTotPosOutput,0.000001),'r:','MarkerSize',5)
            semilogy(Time,max(MeanTotNegOutput,0.000001),'b:','MarkerSize',5)
        elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
            semilogy(Time,max(UQTotPosOutput,0.000001),'r--','LineWidth',0.75)
            semilogy(Time,max(LQTotPosOutput,0.000001),'r--','LineWidth',0.75)
            
            semilogy(Time,max(UQTotNegOutput,0.000001),'b--','LineWidth',0.75)
            semilogy(Time,max(LQTotNegOutput,0.000001),'b--','LineWidth',0.75)
        end
        
        legend([TotPosPlot,TotNegPlot],'+ssRNA','-ssRNA')
        hold off
    end
    font_ax('Total +ssRNA vs Total -ssRNA','Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    
    for i = 59:60
        subplot(1,3,i-57)
        Median = table2array(MedianTable(:,{VarNames{i}}));
        UQ = table2array(UpperQuantileTable(:,{VarNames{i}}));
        LQ = table2array(LowerQuantileTable(:,{VarNames{i}}));
        Mean = table2array(MeanTable(:,{VarNames{i}}));
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            hold on
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
            end
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            hold off
        end
        font_ax(sprintf('%s',VarNames{i}),'Time (h)',MedianTable.Properties.VariableUnits{i},10,'bold',2.0,0)
    end
    
else
    %Individual plots
    for j = 1:length(PlotOutputs)
        CustomTable = MedianTable(:,PlotOutputs(j));
        Median = table2array(MedianTable(:,PlotOutputs(j)));
        UQ = table2array(UpperQuantileTable(:,PlotOutputs(j)));
        LQ = table2array(LowerQuantileTable(:,PlotOutputs(j)));
        Mean = table2array(MeanTable(:,PlotOutputs(j)));
        figure
        if strcmpi(PlotType,'linear') == 1
            hold on
            plot(Time,Median,'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                plot(Time,Mean,'r:','LineWidth',2)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                plot(Time,UQ,'r--','LineWidth',0.75)
                plot(Time,LQ,'r--','LineWidth',0.75)
            end
            hold off
        elseif strcmpi(PlotType,'semilog') == 1
            hold on
            semilogy(Time,max(Median,0.000001),'r-','LineWidth',2)
            if strcmpi(PopulationSetting,'Population') %Population mode plots the mean as well
                semilogy(Time,max(Mean,0.000001),'r:','LineWidth',0.75)
            elseif strcmpi(PopulationSetting,'Single Cell') %Only plot quantiles in Single-Cell mode
                semilogy(Time,max(UQ,0.000001),'r--','LineWidth',0.75)
                semilogy(Time,max(LQ,0.000001),'r--','LineWidth',0.75)
            end
            hold off
        end
        font_ax(sprintf('%s at an MOI of %s',PlotOutputs{j},MOI),'Time (h)',CustomTable.Properties.VariableUnits,10,'bold',2.0,0)
    end
end

end
