%% I. Running basic infection simulations

%Simulate an infection at a multiplicity of infection (MOI) of 5
%Without additional options, this will simulate a single cell infection for 16 hours 
%with both host anti-viral interferon response and viral antagonism of host response.
%The model will output the median, mean, 95th quantile, and 5th quantile of 100
%simulations as well as species plots showing the median and 90% CI of the simulations. For more details see methods.

[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ]= CVB3ODEEval(5);

%The outputted results tables contain metadata describing the run parameters, variable names, and variable units.

%Metadata can be accessed with the following commands:
medianResultsTable.Properties.Description
medianResultsTable.Properties.VariableNames
medianResultsTable.Properties.VariableUnits

%Simulate an infection of a single cell at an MOI of 10 for 24 hours
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24);
summary(medianResultsTable)

%Simulate an infection of a population of cells at an MOI of 10 for 24 hours
%Plots generated in Population mode show the median and mean model traces
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'PopulationSetting','Population','MaxTime',24);
summary(medianResultsTable)

%Change the number of model simulations from which the median and 90% confidence interval are calculated.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'RunCount',10);

%Change the confidence interval that is calculated from the simulation results.
[medianResultsTable,meanResultsTable,Q90ResultsTable,Q10ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'UpperQuantile',0.9,'LowerQuantile',0.1);

%Change the CV of the parameters. This will alter the log-normal distribution from which the parameter is randomly chosen.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'CV',0.3);

%Plot the results on a linear y-axis.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'PlotType','linear');

%% II. Running infection simulations where host-viral interactions are perturbed

%Simulate an infection at an MOI of 10 for 24 hours with the host anti-viral interferon response disabled
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24, 'IFNSwitch','off');
summary(medianResultsTable)

%Simulate an infection at an MOI of 10 for 24 hours with viral antagonism of host anti-viral response disabled
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'VirResponse','off');
summary(medianResultsTable)

%Simulate an infection at an MOI of 10 for 24 hours with exogenous interferon stimulation 2 hours after infection.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'IFNStimulation','on','IFNStimulationTime',2);
summary(medianResultsTable)

%Simulate an infection at an MOI of 10 for 24 hours with exogenous interferon pre-stimulation 2 hours before infection.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'IFNStimulation','on','IFNStimulationTime',-2);
summary(medianResultsTable)

%% III. Exporting simulation results

%Each results table will be exported to the current working directory in a separate file.
%The run parameters for the simulation will be outputted as a separate file named RunParameters_<myfile> for reference.

%Export simulation results to a csv file.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'ExportData','on','DataFile','myResults.csv');

%Export simulation results to a text file.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'ExportData','on','DataFile','myResults.txt');

%Export simulation results to a xls file.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'ExportData','on','DataFile','myResults.xls');

%% IV. Running model sensitivity analysis

%Model parameter sensitivity in Single Cell mode at the specifed MOI and MaxTime can be assessed using the following command:
%By default a heatmap of the log2 fold change in total +ssRNA ('Plus RNA') from the base model results will be generated.
%By default base model results are also plotted and outputted.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, sensitivityTable] = CVB3ODEEval(10,'MaxTime',8,'SensitivityAnalysis','on');

%Generate sensitivity heatmap without base model plots
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, sensitivityTable] = CVB3ODEEval(10,'MaxTime',8,'SensitivityAnalysis','on','PlotResults','off');

%Assess model -ssRNA sensitivity (SensitivityAnalysisOutput = 'Minus RNA') to 10-fold changes in parameter values (ScalingFactor = 10).
%Generate data for 5 orders of magnitude above and below base model parameter values (MaxScalingOrder = 5).
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, sensitivityTable] = CVB3ODEEval(10,'MaxTime',24,'PlotResults','off','SensitivityAnalysis','on','SensitivityAnalysisOutput','Minus RNA','ScalingFactor',10,'MaxScalingOrder',5);

%Export heatmap data to csv file.
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, sensitivityTable] = CVB3ODEEval(10,'MaxTime',24,'PlotResults','off','SensitivityAnalysis','on','SensitivityAnalysisOutput','Minus RNA','ScalingFactor',10,'MaxScalingOrder',2,...
'ExportData','on','DataFile','mySensitivityAnalysis.csv');

%% V. Custom Plotting

%Turn off default plotting
[medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable, ~ ] = CVB3ODEEval(10,'MaxTime',24,'PlotResults','off');

%Plot only the desired species by providing the CVB3ODEPlots function with a cell array containing the desired variables.
%Variable names available for plotting can be found by looking that the resultsTable variable names property.
medianResultsTable.Properties.VariableNames
%Plots will default to semilog scale.
CVB3ODEPlots(medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable,'PlotOutputs',{'Total +ssRNA' 'Virions'});

%Plots can be generated on a linear scale using the 'PlotType' option
CVB3ODEPlots(medianResultsTable,meanResultsTable,Q95ResultsTable,Q05ResultsTable,'PlotOutputs',{'Total +ssRNA' 'Virions'},'PlotType','linear');