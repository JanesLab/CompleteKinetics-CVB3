function [SensitivityAnalysis] = SensitivityAnalysisfunc(CVB3ODE,LabeledConstants,LabeledInitVals,SensitivityAnalysisOutput,OtherSpecies,ScalingFactor,MaxScalingOrder,IFNSwitch,VirResponse,IFNStimulation,IFNStimulationTime,MaxTime,Options)
%Function: Runs sensitivity analysis on the model
%
%INPUTS:
    %LabeledConstants: Cell array of constant labels and values
    %LabeledInitVals: Cell array of initial condition labels and values
    %SensitivityAnalysisOutput: Indicates the model species that is evaluated during sensitivity analysis.
    %OtherSpecies: Vector of strings of species to study if the 'Other' SensitivityAnalysisOutput is selected.
    %ScalingFactor: The factor by which each parameter is scaled during sensitivity analysis
    %MaxScalingOrder: The integer absolute value of the maximum order of magnitude change to test
    %IFNSwitch: 'on' to enable IFN response, otherwise it is disabled.
    %IFNStimulation: Indicates if IFN pre-stimulation is simulated.
    %IFNStimulationTime: If IFNStimulation is enabled, indicates the time at which exogenous IFN stimulation occurs. Negative in the case of pre-stimulation.
    %MaxTime: Maximum time for the solver in hours. Used for plotting.
    %Options: Desired optional settings for the ODE solver.
%  
%OUTPUTS:
    %SensitivityAnalysis: A table containing the final values from the model run
        
%% Set-up for the sensitivity analysis    

%Sets the scaling factor that we vary each parameter by:
Scalar = ScalingFactor;

%Pulls Constants and InitVals from their labeled cell arrays
Constants = LabeledConstants{:,:};
InitVals = LabeledInitVals{:,:};
ConstantLabels = LabeledConstants.Properties.RowNames;
InitValLabels = LabeledInitVals.Properties.RowNames;

%These are the labels of potentially variable initial conditions to test:
InitValLabelVec = {'uCVB3' 'Defective uCVB3' 'uDAF' 'uCAR' 'Viral Ribosomes' 'Host Ribosomes' 'ISG Protein'};

%Pre-defining quick options for OutputIndex:
Positive_Strands = {'R p endo' 'Defective R p endo' 'R p cyt' 'Defective R p cyt' 'T c' 'R p VRO' 'RNA-Bound Pentamer' 'Virion' 'P2Filled' 'P3Filled' 'P4Filled' 'P5Filled' 'P6Filled' 'P7Filled' 'P8Filled' 'P9Filled' 'P10Filled' 'P11Filled'};
Negative_Strands = {'R n cyt' 'R n VRO'};
dsRNA = {'R Ip VRO' 'R In VRO'};
Polyproteins = {'Protease'};
Virion = {'Virion'};
EmptyProvirion = {'Empty Provirion'};
    
%Assigns OutputIndex based on the SensAnalysis input. This is a vector of indices indicating the desired species of study.
if strcmpi(SensitivityAnalysisOutput,'Plus RNA')
    OutputIndex = Positive_Strands;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Positive Strand RNA Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'Minus RNA')
    OutputIndex = Negative_Strands;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Negative Strand RNA Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'dsRNA')
    OutputIndex = dsRNA;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Double Strand RNA Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'RNA Ratio')
    OutputIndex = Positive_Strands;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Positive/Negative Strand RNA Concentration Ratio at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'Polyprotein')
    OutputIndex = Polyproteins;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Polyprotein Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'Virion')
    OutputIndex = Virion;    
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Virion Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'Empty Provirion')
    OutputIndex = EmptyProvirion;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes in Empty Capsid Concentration at %0.1f hpi';
elseif strcmpi(SensitivityAnalysisOutput,'Other')
    OutputIndex = OtherSpecies;
    Title = 'Comprehensive Sensitivity Analysis: Log2 Fold Changes of Custom Input Concentration at %0.1f hpi';
else
    fprintf('\nError: Invalid SensitivityAnalysisOutput Input.\n')
    return
end

%% Sensitivity analysis:

SensMat = zeros(1+2*MaxScalingOrder,length(Constants)+length(InitValLabelVec)); %Initializes the matrix (Excludes volume conversion constants)

for j = 1:length(Constants) + length(InitValLabelVec)
    
    fprintf('\n %2d \n',j) %Displays index for tracking of error messages
    
    for i = 1:(1+2*MaxScalingOrder)
        
        fprintf('%1d',i) %Displays index for tracking
        
        %Conditional to alter hill coefficient analysis. Order of magnitude changes are unrealistic and result in solver failure or large delays, so integer/scaler changes are used
        if j <= length(Constants) && j == find(strcmpi(LabeledConstants.Properties.RowNames,'Hill_Constant'))
            HillConstantIndex = j;
            if i == 1
                OriginalValue = Constants(j); %Saves original value for resetting later
                Constants(j) = OriginalValue / 2^MaxScalingOrder; %Brings down to minimum value
            elseif i <= MaxScalingOrder
                Constants(j) = Constants(j) * 2; %Scales each number below the median by 2 up to the original value
            else
                Constants(j) = Constants(j) + 1; %Increases by one integer value
            end 
        elseif j <= length(Constants) %This condition is to manipulate the values of all constants other than the Hill coefficient
            if i == 1
                OriginalValue = Constants(j); %Saves original value for resetting later
                Constants(j) = OriginalValue / Scalar^(MaxScalingOrder); %Brings down to minimum value
            else
                Constants(j) = Constants(j) * Scalar; %Increases by an order of magnitude
            end
        else %This condition is to manipulate the value of initial conditions
            if i == 1
                OriginalValue = LabeledInitVals{InitValLabelVec{j-length(Constants)},:}; %Saves original value for resetting later
                ValueIndex = find(strcmpi(LabeledInitVals.Properties.RowNames,InitValLabelVec{j-length(Constants)})); %Saves original value index for resetting later
                InitVals(ValueIndex) = OriginalValue / Scalar^(MaxScalingOrder); %Brings down to minimum value
            else
                InitVals(ValueIndex) = InitVals(ValueIndex) * Scalar; %Increases by an order of magnitude
            end
        end
        
        [~,Solutions] = ode15s(@(t,y)CVB3ODE(t,y,Constants,IFNSwitch,VirResponse,IFNStimulation,IFNStimulationTime),[0 MaxTime],InitVals,Options); %Solving the ODE
       
        LabeledSolutions = table(Solutions(end,:)','RowNames',InitValLabels,'VariableNames',{'Value'}); %Adding labels to the solutions matrix. Only uses end values.
        
        for k = 1:length(OutputIndex) %Takes output indices and sums end values in each species
            Value = LabeledSolutions{OutputIndex(k),:};
            SensMat(i,j) = SensMat(i,j) + Value;
        end        
        if strcmpi(SensitivityAnalysisOutput,'RNA Ratio') %Necessary to do Sensitivity analysis on the RNA ratio
            Denominator = 0;
            for l = 1:length(Negative_Strands)
                Value = LabeledSolutions{Negative_Strands(l),:};
                Denominator = Denominator + Value;
            end
            SensMat(i,j) = SensMat(i,j)/Denominator; %Creates the ratio
        end
        
        if i == (1+2*MaxScalingOrder) %For resetting original value when max value has been reached
            if j <= length(Constants) %Resets original value for constants
                Constants(j) = OriginalValue;
            else %Resets original value for initial conditions
                InitVals(ValueIndex) = OriginalValue;
            end
        end
 
    end
end

%Min and Max values of the solution in question. Used for colormap axis.
MaxSolution = max(max(SensMat)); 
MinSolution = min(min(SensMat));
UnchangedSolution = SensMat(1+MaxScalingOrder,1); %Value if the model is run without changing any parameters.

SensMat = flipud(real(SensMat)); %Eliminates imaginary numbers and flips the matrix for the later heatmap function
SensMatNorm = log2(SensMat./UnchangedSolution); %Sets the value at original conditions as 0 and puts others relative to that

fprintf('\nThe Minimum Value is %3.2e\n',MinSolution) %Displays minimum of the matrix
fprintf('The Maximum Value is %3.2e\n',MaxSolution) %Displays maximum of the matrix
fprintf('The Average Value is %3.2e\n',UnchangedSolution) %Displays the average value of the matrix

%% Heat Map Generation

%Generates string labels for the heatmap rows based on MaxScalingOrder
RowLabels = cell(1,1+2*MaxScalingOrder); RowMinValue = Scalar^-MaxScalingOrder;
for i = 1:1+2*MaxScalingOrder
    RowLabels{i}= num2str(RowMinValue);
    RowMinValue = Scalar * RowMinValue;
end

%Generates string labels for the rows of the Hill Constant heatmap based on MaxScalingOrder
RowLabelsHillConstant = cell(1,1+2*MaxScalingOrder); RowHillConstantValue = 1/2^MaxScalingOrder;
for i = 1:1+2*MaxScalingOrder
    RowLabelsHillConstant{i}= num2str(RowHillConstantValue);
    if i <= MaxScalingOrder
        RowHillConstantValue = RowHillConstantValue * 2;
    else
        RowHillConstantValue = RowHillConstantValue + 1;
    end
end

%Determines the color axis settings- selects which relative difference is greater for caxislimit
if abs(MaxSolution/UnchangedSolution) >= abs(UnchangedSolution/MinSolution)
    caxislimit = MaxSolution/UnchangedSolution;
else
    caxislimit = UnchangedSolution/MinSolution;
end

%Creates row and column labels for the heatmap
RowLabels = fliplr(RowLabels); %Flips the row label vector horizontally so it aligns with the heatmap data
RowLabelsHillConstant = fliplr(RowLabelsHillConstant); %Flips the row label vector horizontally so it aligns with the heatmap data
ColLabels = strrep([ ConstantLabels' InitValLabelVec ],'_',' '); %Creates vector of string labels for the heatmap columns

%Vector for the heatmap color bar
Colors = [.3137 .3647 .4549; .3333 .4706 .6235; .3922 .5686 .7098; .4902 .6353 .7686; .6157 .7255 .8314; ...
    .9647 .9647 .9686; .9412 .7059 .5373; .8902 .6157 .4667; .8353 .5333 .4039; .7176 .3843 .3294; .5020 .3373 .3255];

%Creates a trimmed version without the middle row of all unchanged solutions
SensMatNorm(1+MaxScalingOrder,:) = [];
RowLabels(:,1+MaxScalingOrder) = [];
RowLabelsHillConstant(:,1+MaxScalingOrder) = [];

%Alters main heatmap and labels by removing Hill Constant for its own subplot with a separate y-axis
HillConstantSensMatNorm = SensMatNorm(:,HillConstantIndex);
SensMatNorm(:,HillConstantIndex) = [];
ColLabelHillConstant = ColLabels(HillConstantIndex);
ColLabels(HillConstantIndex) = [];

%Heatmap plot (Without center row of 1x)
figure
subplot(1,10,1:9)
    heatmap(ColLabels,RowLabels,SensMatNorm,'Colormap',Colors,'Title',sprintf(Title,MaxTime),'XLabel','Constant Name','YLabel','Multiplication Factor','MissingDataColor',[0 0 0],'MissingDataLabel','No Data')
    caxis(round([-log2(caxislimit) log2(caxislimit)],0))
    colorbar('off')
subplot(1,10,10)
    heatmap(ColLabelHillConstant,RowLabelsHillConstant,HillConstantSensMatNorm,'Colormap',Colors,'XLabel','Constant Name','YLabel','Multiplication Factor (< Median), Summand (> Median)','MissingDataColor',[0 0 0],'MissingDataLabel','No Data')
    caxis(round([-log2(caxislimit) log2(caxislimit)],0))

%Assign heatmap data to results table
SensitivityAnalysis = array2table(SensMatNorm,'VariableNames',ColLabels,'RowNames',RowLabels);

end