function [MedianSolutions,MeanSolutions,UpperQuantileSolutions,LowerQuantileSolutions,SensitivitySolutions] = CVB3ODEEval(MOI,varargin)
%Function: Stores all constants, initial values, runs the ODE, and calls other relevant functions
%
%INPUTS:
    %MOI: Simulated multiplicity of infection or the number of plaque forming units of virus. No default value.
    %PopulationSetting: Indicates whether or not to use a Poisson distribution for MOI or just a single cell. 'Single Cell' or 'Population' mode. Default = 'Single Cell'.
    %MaxTime (h): Maximum time for the solver in hours. Default = 16 h.
    %IFNSwitch: Indicates whether the host-cell anti-viral interferon (IFN) response is enabled. Default = 'on'.
    %VirResponse: Indicates whether the viral anti-interferon response is enabled. Default = 'on'.
    %IFNStimulation: Indicates if exogenous IFN stimulation is simulated. Default = 'off'.
    %IFNStimulationTime: If IFNStimulation is enabled, indicates the time at which exogenous IFN stimulation occurs. Negative in the case of pre-stimulation. Default = 0 h.
    %EC50_RNAdeg (nM): Half-effective concentration of ISG protein for viral RNA degradation. Default = 5 nM.
    %EC50_DetectorDeg (nM): Half-effective concentration of viral Protease for viral RNA detector degradation. Default = 1 pM.
    %EC50_Protease (nM): Half-effective concentration of ISG protein for inhibition of 2Apro by ISG15. Default = 20 nM.
    %EC50_Translate (nM): Half-effective concentration of ISG protein for translation inhibition. Default = 10 nM.
    %PlotResults: Indicates whether results plots are generated after simulations have completed. Default = 'on'.
    %PlotType: Specifies the type of axes to plot results on ('linear' or 'semilog'). Default = 'semilog'.
    %SensitivityAnalysis: Indicates whether sensitivity analysis is conducted on the model parameters. Default = 'off'.
    %SensitivityAnalysisOutput: Indicates the model species ('Plus RNA','Minus RNA','dsRNA','RNA Ratio','Polyprotein','Virion','Empty Provirion','Other') that is evaluated during sensitivity analysis. Default = 'Plus RNA'.
    %ScalingFactor: The factor by which each parameter is scaled during sensitivity analysis. Default = 10.
    %MaxScalingOrder: The integer absolute value of the maximum change to test via sensitivty analysis. Default = 2.
    %ExportData: Indicates whether the data should be exported as DataFile in the current working directory. Default = 'off'.
    %DataFile: A string specifying the name and type of file to be exported. The corresponding file will be written to the current working directory. Default = 'results.csv'.
    %RunCount: The number of times the model is simulated with parameters stochastically sampled from independent log-normal distributions. Default = 100.
    %UpperQuantile: The upper quantile outputted for model species. Default = 0.95.
    %LowerQuantile: The lower quantile outputted for model species. Default = 0.05.
    %CV: The coefficient of variation of the log-normal distribution the from which the parameters are stochastically determined. Default = 0.05.
%
%OUTPUTS:
    %MedianSolutions: A vector containing the median value of all model simulations from the ODE at a given time point.
    %MeanSolutons: A vector containing the mean value of all model simulations from the ODE at a given time point.
    %UpperQuantileSolutions: A vector containing the upper quantile value of all model simulations from the ODE at a given time point.
    %LowerQuantileSolutions: A vector containing the lower quantile value of all model simulations from the ODE at a given time point.
    %SensitivitySolutions: A vector containing the log2 fold change in SensitivityAnalysisOutput of all model simulations from the ODE at a given time point.
    
%% Set option defaults and parse inputs

%Set simulation parameter defaults
defaultMaxTime = 16; %hours
defaultPopulationSetting = 'Single Cell';
defaultIFNSwitch = 'on';
defaultVirResponse = 'on';
defaultIFNStimulation = 'off'; 
defaultIFNStimulationTime = 0; %hours
defaultEC50_RNAdeg = 5; %nM
defaultEC50_DetectorDeg = 0.001; %nM
defaultEC50_Protease = 20; %nM
defaultEC50_Translate = 10; %nM
defaultPlotResults = 'on';
defaultPlotType = 'semilog';
defaultSensitivityAnalysis = 'off'; 
defaultSensitivityAnalysisOutput = 'Plus RNA'; 
defaultScalingFactor = 10; 
defaultMaxScalingOrder = 2;
defaultExportData = 'off';
defaultDataFile = 'results.csv';
defaultRunCount = 100;
defaultUpperQuantile = 0.95;
defaultLowerQuantile = 0.05;
defaultCV = 0.05;

CVB3ODEEvalparser = inputParser;

%set general criteria for inputs
   validNum = @(x) isnumeric(x);
   validPosNum = @(x) isnumeric(x) && (x > 0);
   validPosNumOrZero = @(x) isnumeric(x) && (x >= 0);
   validSwitchString = @(x) strcmpi(x,'on') || strcmpi(x,'') || strcmpi(x,'off');
   validPopulationString = @(x) strcmpi(x,'Single Cell') || strcmpi(x,'Population');
   validSensString = @(x) strcmpi(x,'Plus RNA') || strcmpi(x,'Minus RNA') || strcmpi(x,'dsRNA') || strcmpi(x,'RNA Ratio') ||...
       strcmpi(x,'Polyprotein') || strcmpi(x,'Virion') || strcmpi(x,'Empty Provirion') || strcmpi(x,'Other') || strcmpi(x,'');
   validChar = @(x) ischar(x);
   validPercent = @(x) isnumeric(x) && (x > 0) && (x <= 1);
   validPlotType = @(x) strcmpi(x,'linear') || strcmpi(x,'semilog');
   
%set parser options and valid input criteria
   addRequired(CVB3ODEEvalparser,'MOI',validPosNumOrZero);
   
   addOptional(CVB3ODEEvalparser,'PopulationSetting',defaultPopulationSetting,validPopulationString);
   addOptional(CVB3ODEEvalparser,'MaxTime',defaultMaxTime,validPosNum);
   addOptional(CVB3ODEEvalparser,'IFNSwitch',defaultIFNSwitch,validSwitchString);
   addOptional(CVB3ODEEvalparser,'VirResponse',defaultVirResponse,validSwitchString);
   addOptional(CVB3ODEEvalparser,'IFNStimulation',defaultIFNStimulation,validSwitchString);
   addOptional(CVB3ODEEvalparser,'IFNStimulationTime',defaultIFNStimulationTime,validNum);
   addOptional(CVB3ODEEvalparser,'EC50_RNAdeg',defaultEC50_RNAdeg,validPosNumOrZero);
   addOptional(CVB3ODEEvalparser,'EC50_DetectorDeg',defaultEC50_DetectorDeg,validPosNumOrZero);
   addOptional(CVB3ODEEvalparser,'EC50_Protease',defaultEC50_Protease,validPosNumOrZero);
   addOptional(CVB3ODEEvalparser,'EC50_Translate',defaultEC50_Translate,validPosNumOrZero);
   addOptional(CVB3ODEEvalparser,'PlotResults',defaultPlotResults,validSwitchString);
   addOptional(CVB3ODEEvalparser,'PlotType',defaultPlotType,validPlotType);
   addOptional(CVB3ODEEvalparser,'SensitivityAnalysis',defaultSensitivityAnalysis,validSwitchString);
   addOptional(CVB3ODEEvalparser,'SensitivityAnalysisOutput',defaultSensitivityAnalysisOutput,validSensString);
   addOptional(CVB3ODEEvalparser,'ScalingFactor',defaultScalingFactor,validPosNum);
   addOptional(CVB3ODEEvalparser,'MaxScalingOrder',defaultMaxScalingOrder,validPosNum);
   addOptional(CVB3ODEEvalparser,'ExportData',defaultExportData,validSwitchString);
   addOptional(CVB3ODEEvalparser,'DataFile',defaultDataFile,validChar);
   addOptional(CVB3ODEEvalparser,'RunCount',defaultRunCount,validPosNum);
   addOptional(CVB3ODEEvalparser,'UpperQuantile',defaultUpperQuantile,validPercent);
   addOptional(CVB3ODEEvalparser,'LowerQuantile',defaultLowerQuantile,validPercent);
   addOptional(CVB3ODEEvalparser,'CV',defaultCV,validPosNumOrZero);

   parse(CVB3ODEEvalparser,MOI,varargin{:});

%Extract parsed inputs
MOI = CVB3ODEEvalparser.Results.MOI;
PopulationSetting = CVB3ODEEvalparser.Results.PopulationSetting;
MaxTime = CVB3ODEEvalparser.Results.MaxTime;
IFNSwitch = CVB3ODEEvalparser.Results.IFNSwitch;
VirResponse = CVB3ODEEvalparser.Results.VirResponse;
IFNStimulation = CVB3ODEEvalparser.Results.IFNStimulation;
IFNStimulationTime = CVB3ODEEvalparser.Results.IFNStimulationTime;
EC50_RNAdeg = CVB3ODEEvalparser.Results.EC50_RNAdeg;
EC50_DetectorDeg = CVB3ODEEvalparser.Results.EC50_DetectorDeg;
EC50_Protease = CVB3ODEEvalparser.Results.EC50_Protease;
EC50_Translate = CVB3ODEEvalparser.Results.EC50_Translate;
PlotResults = CVB3ODEEvalparser.Results.PlotResults;
SensitivityAnalysis = CVB3ODEEvalparser.Results.SensitivityAnalysis;
SensitivityAnalysisOutput = CVB3ODEEvalparser.Results.SensitivityAnalysisOutput;
ScalingFactor = CVB3ODEEvalparser.Results.ScalingFactor;
MaxScalingOrder = CVB3ODEEvalparser.Results.MaxScalingOrder;
ExportData = CVB3ODEEvalparser.Results.ExportData;
DataFile = CVB3ODEEvalparser.Results.DataFile;

%If senstivity analysis is specified, override the RunCount to be 1 and the standard deviation to be 0.
%This ensures that we are only doing sensitivity analysis on the average parameter values.
if strcmpi(SensitivityAnalysis,'on')
   RunCount = 1;
   CV=0;
else
   RunCount = CVB3ODEEvalparser.Results.RunCount;
   CV = CVB3ODEEvalparser.Results.CV;
end

UpperQuantile = CVB3ODEEvalparser.Results.UpperQuantile;
LowerQuantile = CVB3ODEEvalparser.Results.LowerQuantile;
plotType=CVB3ODEEvalparser.Results.PlotType;

%% Process Inputs
%Prompt the user to enter additional information for custom sensitivity analysis    
if strcmpi(SensitivityAnalysisOutput,'Other')
    fprintf(['The options are: uCVB3, Defective uCVB3, uDAF, bDAF, Defective bDAF, uDAF TJ, bDAF TJ, Defective bDAF TJ, uCVB3 TJ, Defective uCVB3 TJ, uCAR, bCAR, Defective bCAR, R p endo, Defective R p endo,\n' ...
    'R p cyt, Defective R p cyt, T c, R n cyt, VRO Avail, R p VRO, R n VRO, Pol3D VRO, R Ip VRO, R In VRO, ATPase2C,\n' ...
    'Pentamer cyt, Pentamer VRO, RNA-Bound Pentamer, P2Filled, P3Filled, P4Filled, P5Filled, P6Filled, P7Filled, P8Filled, P9Filled, P10Filled, P11Filled, Virion,\n' ...
    'P2Empty, P3Empty, P4Empty, P5Empty, P6Empty, P7Empty, P8Empty, P9Empty, P10Empty, P11Empty, Empty Provirion,\n'
    'Viral Ribosomes, Host Ribosomes, Protease, ISG Protein\n'])
    OtherSpecies = input('Enter the name(s) of the species you wish to study as a vector of their names as strings separated by spaces.');   
else
    OtherSpecies = [];
end

%Ensure IFNStimulationTime is set to 0 if IFNStimulation is off 
if strcmpi(IFNStimulation,'off')||strcmpi(IFNStimulation,'')
    IFNStimulationTime = 0;
else
end

%Convert parameter means to mu and CV to sigma for use in log-normal random number generation for
%constants.

function mu = muFunc(M,V)
    mu = log((M^2)/sqrt(V+M^2));
end
sigma = sqrt(log(CV^2+1));

disp('Run Count: ') %Display the label of the Run Count progress tracker

%% Begin repeated simulations 
for i = 1:RunCount
%% Constants for the ODE (Justification for these constants can be found in the parameter table)
%Cell and virus non-rate constant parameters
CellVolume = 3700; %(um^3) 
NucVolume = 960; %(um^3) 
CellSurfaceArea = 2200; %(µm^2) 
CellCytVolume = CellVolume - NucVolume; %(um^3) 
CVB3diameter = 30; %(nm) 
Pol3DLength = 7.1; %(nm) 
MaxCVB3Conc = (6.022e23)^-1 / (4/3*pi*(CVB3diameter/2)^3) * 1e33; %(nM)
ParticlePFUConv = 800; %Number
VROSurfaceArea = 120; %(µm^2) 
VROSurfaceVolume = VROSurfaceArea * Pol3DLength/1000; %(µm^3) 
ISGbpLength = 0.495; %(kb) 
k_transcribeISG = lognrnd(muFunc(170,(CV*170)^2),sigma); %(kb/h) 
RiboActive = 0.8; %Fraction 
RiboTotal = 1e5; %Number
PolysomeSize = 2.5; %Average polysome size. 
TranslationRate = lognrnd(muFunc(3.5*3600,(CV*3.5*3600)^2),sigma); %(aa/h) 
PolymeraseRate = lognrnd(muFunc(50*3600,(CV*50*3600)^2),sigma); %(bp/h)  
CVB3genomelength = 7397; %(bp) 
CVB3polyprolength = 2185; %(aa) 
EncapsidatedNegStrandRatio = 1/1790; %Fraction

%Receptor binding constants
k_on_DAF = lognrnd(muFunc(1.47e5*3600/1e9,(CV*3600/1e9)^2),sigma); %(nM^-1 * h^-1) 
k_off_DAF = lognrnd(muFunc(0.3*3600,(CV*3600)^2),sigma); %(h^-1)
k_transloc = lognrnd(muFunc(60,(CV*60)^2),sigma); %(h^-1) 
k_on_CAR = lognrnd(muFunc(3.35e3*3600/1e9,(CV*3.35e3*3600/1e9)^2),sigma); %(nM^-1 * h^-1)
k_off_CAR = lognrnd(muFunc(8.21e-4*3600,(CV*8.21e-4*3600)^2),sigma); %(h^-1) 
k_internal = lognrnd(muFunc(0.5,(CV*0.5)^2),sigma); %(h^-1) 
k_endo_escape = lognrnd(muFunc(420,(CV*420)^2),sigma); %(h^-1)  
InternalVolConv = CellSurfaceArea / (6.022e23 * CVB3diameter^2 * MaxCVB3Conc * CellCytVolume) * 1e30; %(dimensionless)

%Replication constants 
VROFormThreshold = 25; %(molecules)
k_T_c_Form = lognrnd(muFunc(25/PolysomeSize,(CV*25/PolysomeSize)^2),sigma); %(h^-1 * nM^-1)
k_Translate = TranslationRate * PolysomeSize / CVB3polyprolength; %(h^-1) 
k_P_on = lognrnd(muFunc(1,(CV*1)^2),sigma); %(h^-1) 
k_P_off = lognrnd(muFunc(1,(CV*1)^2),sigma); %(h^-1) 
k_N_on = lognrnd(muFunc(1,(CV*1)^2),sigma); %(h^-1) 
k_N_off = lognrnd(muFunc(1,(CV*1)^2),sigma); %(h^-1) 
VROVolConv = CellCytVolume/VROSurfaceVolume; %(dimensionless)
k_R_Ip_Form = lognrnd(muFunc(0.36,(CV*0.36)^2),sigma); %(nM^-1 h^-1) 
k_N_Transcript = PolymeraseRate / CVB3genomelength; %(h^-1) 
k_P_Transcript = PolymeraseRate / CVB3genomelength; %(h^-1) 
k_R_In_Form = lognrnd(muFunc(0.36,(CV*0.36)^2),sigma); %(nM^-1 h^-1) 
u_P_cyt = lognrnd(muFunc(0.23,(CV*0.23)^2),sigma); %(h^-1) 
u_P_VRO = lognrnd(muFunc(0.23/VROVolConv,(CV*0.23/VROVolConv)^2),sigma); %(h^-1) 
u_N_cyt = lognrnd(muFunc(0.23,(CV*0.23)^2),sigma); %(h^-1) 
u_N_VRO = lognrnd(muFunc(0.23/VROVolConv,(CV*0.23/VROVolConv)^2),sigma); %(h^-1) 
u_T_c = 0; %(h^-1) 
u_VirProt_cyt = lognrnd(muFunc(0.05,(CV*0.05)^2),sigma); %(h^-1)
u_VirProt_VRO = lognrnd(muFunc(0.05/VROVolConv,(CV*0.05/VROVolConv)^2),sigma); %(h^-1) 

%Capsid assembly constants:
k1f_cap = lognrnd(muFunc(3.6,(CV*3.6)^2),sigma); %(nM^-1 h^-1) 
k1b_cap = lognrnd(muFunc(3.6e3,(CV*3.6e3)^2),sigma); %(h^-1) 
k_Pentamer_on = lognrnd(muFunc(0.36,(CV*0.36)^2),sigma); %(nM^-1 h^-1)
k_Pentamer_off = lognrnd(muFunc(36,(CV*36)^2),sigma); %(h^-1) 
kRNACapBind = lognrnd(muFunc(3.6,(CV*3.6)^2),sigma); %(nM^-1 h^-1)
kRNACapUnbind = lognrnd(muFunc(3.6e3,(CV*3.6e3)^2),sigma); %(h^-1) 
k2f_cap = lognrnd(muFunc(3.6,(CV*3.6)^2),sigma); %(nM^-1 h^-1) 
k2b_cap = lognrnd(muFunc(3.6e3,(CV*3.6e3)^2),sigma); %(h^-1) 
u_cap_cyt = lognrnd(muFunc(0.05,(CV*0.05)^2),sigma); %(h^-1) 
u_cap_VRO = lognrnd(muFunc(0.05/VROVolConv,(CV*0.05/VROVolConv)^2),sigma); %(h^-1) 

%Feedback loop constants
kcat_cleave = lognrnd(muFunc(0.3216*3600,(CV*0.3216*3600)^2),sigma); %(h^-1)
Km_cleave = lognrnd(muFunc(960e3,(CV*960e3)^2),sigma); %(nM)
InitRibAvail = (1 - RiboActive) * RiboTotal / PolysomeSize / (6.022e23 * CellCytVolume) * 1e24; %(nM)
kD_Hill = lognrnd(muFunc(1,(CV*1)^2),sigma); %(nM)
n_Hill = lognrnd(muFunc(1.36,(CV*1.36)^2),sigma); %(dimensionless)
ISGformRate = k_transcribeISG/ISGbpLength; %(h^-1)
ISGBasal = 0; %(nM) 
u_ISG = lognrnd(muFunc(0.1,(CV*0.1)^2),sigma); %(h^-1)
EC50_Protease = lognrnd(muFunc(EC50_Protease,(CV*EC50_Protease)^2),sigma); %(nM)
EC50_Translate = lognrnd(muFunc(EC50_Translate,(CV*EC50_Translate)^2),sigma); %(nM)
OAS_RNAdeg = lognrnd(muFunc(5,(CV*5)^2),sigma); %(nM)
EC50_RNAdeg = lognrnd(muFunc(EC50_RNAdeg,(CV*EC50_RNAdeg)^2),sigma); %(nM)
EC50_DetectorDeg = lognrnd(muFunc(EC50_DetectorDeg,(CV*EC50_DetectorDeg)^2),sigma); %(nM)

%Creates vector of all constants needed for rate equations. 
Constants = [k_on_DAF,k_off_DAF,k_transloc,k_on_CAR,k_off_CAR,k_internal,InternalVolConv,k_endo_escape ...
            PolysomeSize,k_T_c_Form,k_Translate,k_P_on,k_P_off,k_N_on,k_N_off,VROVolConv,k_R_Ip_Form,k_N_Transcript,k_P_Transcript,k_R_In_Form ...
            u_P_cyt,u_P_VRO,u_N_cyt,u_N_VRO,u_T_c,u_VirProt_cyt,u_VirProt_VRO ...
            k1f_cap,k2f_cap,k1b_cap,k2b_cap,kRNACapBind,kRNACapUnbind,k_Pentamer_on,k_Pentamer_off,u_cap_cyt,u_cap_VRO ...
            kcat_cleave,Km_cleave,InitRibAvail,u_ISG,ISGformRate,ISGBasal,kD_Hill,n_Hill,EC50_Protease,EC50_Translate,OAS_RNAdeg,EC50_RNAdeg,EC50_DetectorDeg,VROFormThreshold];
 
%Creates vector of labels for the above constants.      
ConstantLabels = {'k_on_DAF' 'k_off_DAF' 'k_transloc' 'k_on_CAR' 'k_off_CAR' 'k_internal' 'Cyt_Vol._Conv.'  'k_endo_escape'  ...
    'PolysomeSize' 'k_Tc_Form' 'k_Translate' 'k_P_on' 'k_P_off' 'k_N_on' 'k_N_off' 'VRO_Vol._Conv.' 'k_RIp_Form' 'k_N_Transcribe' 'k_P_Transcribe' 'k_RIn_Form'  ...
    'u_P_cyt' 'u_P_VRO' 'u_N_cyt' 'u_N_VRO' 'u_T_c' 'u_VirProt_cyt' 'u_VirProt_VRO' ...
    'k1f_cap' 'k2f_cap' 'k1b_cap' 'k2b_cap' 'k_RNACapBind' 'k_RNACapUnbind' 'k_Pentamer_on' 'k_Pentamer_off' 'u_cap_cyt' 'u_cap_VRO' ...
    'kcat_cleave' 'Km_cleave' 'Initial_Viral_Ribosomes' 'u_ISG' 'ISGformRate' 'ISGBasal' 'kD_Hill' 'Hill_Constant' 'EC50_Protease' 'EC50_Translate' 'OAS_RNAdeg' 'EC50_RNAdeg' 'EC50_DetectorDeg' 'VRO_Formation_Threshold'};

%Creates a labeled table of constants
LabeledConstants = table(Constants','RowNames',ConstantLabels,'VariableNames',{'Value'});

%% Initial conditions for the ODE

%Viral Particles
if strcmpi(PopulationSetting,'Population') %Population mode
    uCVB3 = poissrnd(MOI) * (CVB3diameter/1000)^2 / CellSurfaceArea * MaxCVB3Conc; %(nM) Unbound infectious CVB3 virions at cell surface
else %Defaults to single-cell Mode
    uCVB3 = MOI * (CVB3diameter/1000)^2 / CellSurfaceArea * MaxCVB3Conc; %(nM) Unbound infectious CVB3 virions at cell surface
end
uCVB3_Defective = ParticlePFUConv * uCVB3; %(nM) Unbound defective CVB3 virions at cell surface

%Receptors and Internalization
uDAF = 4.4e4 / (6.022e23) / (CVB3diameter/1000 * CellSurfaceArea) * 1e24; %(nM) Unbound DAF molecules
bDAF = 0; %(nM) CVB3-bound DAF molecules
bDAF_Defective = 0; %(nM) Defective CVB3-Bound DAF not in tight junctions
uDAF_TJ = 0; %(nM) Unbound DAF in tight junctions
bDAF_TJ = 0; %(nM) CVB3-Bound DAF in tight junctions
bDAF_Defective_TJ = 0; %(nM) Defective CVB3-Bound DAF in tight junctions
uCVB3_TJ = 0; %(nM) Unbound CVB3 in tight junctions
uCVB3_Defective_TJ = 0; %(nM) Unbound defective CVB3 virions in tight junctions
uCAR = 5.5e6 / (6.022e23) / (CVB3diameter/1000 * CellSurfaceArea) * 1e24; %(nM) Unbound CAR in tight junctions
bCAR = 0; %(nM) CVB3-Bound CAR in tight junctions
bCAR_Defective = 0; %(nM) Defective CVB3-Bound CAR in tight junctions
R_p_endo = 0; %(nM) Internalized virus that has not yet been unpackaged
R_p_endo_Defective = 0; %(nM) Interalized defective virus that has not yet been unpackaged

%Replication
R_p_cyt = 0; %(nM) Positive-strand RNA in the cytoplasm
R_p_cyt_Defective = 0; %(nM) Defective positive-strand RNA in the cytoplasm
T_c = 0; %(nM) Translation complexes (Positive strand + Ribosome) in the cytoplasm
R_n_cyt = 0; %(nM) Negative-strand RNA in the cytoplasm
R_p_VRO = 0; %(nM) Positive-strand RNA on the VRO
R_n_VRO = 0; %(nM) Negative-strand RNA on the VRO
Pol3D_VRO = 0; %(nM) Polymerases and associated proteins needed for transcription on the VRO
R_Ip_VRO = 0; %(nM) Positive-strand RNA transcription complexes (Positive strand + Polymerase) on the VRO
R_In_VRO = 0; %(nM) Negative-strand RNA transcription complexes (Negative strand + Polymerase) on the VRO
ATPase2C = 0; %(nM) 2C ATPase molecules on the VRO used for capsid pentamer interaction

%Capsid Assembly
Pentamer_cyt = 0; %(nM) 5 protomers (20 VPs) assembled into much more stable pentamer form in the cytoplasm
Pentamer_VRO = 0; %(nM) 5 protomers (20 VPs) assembled into much more stable pentamer form on the VRO
RNAPentamer = 0; %(nM) Pentamer in cytoplasm that has bound a positive RNA strand
Virion = 0; %(nM) Mature, infectious virions in the cytoplasm
EmptyProvirion = 0; %(nM) Self-assembled 12-pentamer capsid molecules that lack +ssRNA in the cytoplasm
P2Filled = 0; P2Empty = 0; P3Filled = 0; P3Empty = 0; P4Filled = 0; P4Empty = 0; P5Filled = 0; P5Empty = 0; P6Filled = 0; P6Empty = 0; %(nM) Partially assembled capsid stages
P7Filled = 0; P7Empty = 0; P8Filled = 0; P8Empty = 0; P9Filled = 0; P9Empty = 0; P10Filled = 0; P10Empty = 0; P11Filled = 0; P11Empty = 0; %(nM) Partially assembled capsid stages

RibAvail = InitRibAvail; %(nM) Polysome complexes (of 2-3 ribosomes each) initially inactive and accessible to the virus
RibUnavail = RiboActive * RiboTotal / PolysomeSize / (6.022e23 * CellCytVolume) * 1e24; %(nM) Polysome complexes (of 2-3 ribosomes each) initially in-use by host and inaccessible to the virus
Protease = 0; %(nM) Viral proteases in the cytoplasm
ISGProtein = ISGBasal; %(nM) ISG protein in the cytoplasm

%Set the initial ISGProtein concentration in the event of pre-stimulation
if IFNStimulationTime < 0
    VirDetection = 1; %Starts Interferon response at the start of the simulation.
    
    %Runs the ODE to calculate ISG levels. Note: Viral proteases do not need to be accounted for since we're running this simulation before infection begins.
    [~,PrimedISGProtein] = ode15s(@(t,PrimedISGProtein) VirDetection * ISGformRate * RibUnavail - u_ISG * PrimedISGProtein,[0 abs(IFNStimulationTime)],ISGProtein,odeset('NonNegative',1));
    
    ISGProtein = ISGProtein + PrimedISGProtein(length(PrimedISGProtein)); %Sets the initial value of ISGProtein to the final value (and adds back in basal levels)
end   

%Create initial values vector for ODE solver.
InitVals = [uCVB3,uCVB3_Defective,uDAF,bDAF,bDAF_Defective,uDAF_TJ,bDAF_TJ,bDAF_Defective_TJ,uCVB3_TJ,uCVB3_Defective_TJ,uCAR,bCAR,bCAR_Defective,R_p_endo,R_p_endo_Defective ...
 R_p_cyt,R_p_cyt_Defective,T_c,R_n_cyt,R_p_VRO,R_n_VRO,Pol3D_VRO,R_Ip_VRO,R_In_VRO,ATPase2C ...
 Pentamer_cyt,Pentamer_VRO,RNAPentamer,P2Filled,P3Filled,P4Filled,P5Filled,P6Filled,P7Filled,P8Filled,P9Filled,P10Filled,P11Filled,Virion ...
 P2Empty,P3Empty,P4Empty,P5Empty,P6Empty,P7Empty,P8Empty,P9Empty,P10Empty,P11Empty,EmptyProvirion ...
 RibAvail,RibUnavail,Protease,ISGProtein];

%Create labels vector for the above initial values.
InitValLabels = {'uCVB3' 'Defective uCVB3' 'uDAF' 'bDAF' 'Defective bDAF' 'uDAF TJ' 'bDAF TJ' 'Defective bDAF TJ' 'uCVB3 TJ' 'Defective uCVB3 TJ' 'uCAR' 'bCAR' 'Defective bCAR' 'R p endo' 'Defective R p endo' ...
    'R p cyt' 'Defective R p cyt' 'T c' 'R n cyt' 'R p VRO' 'R n VRO' 'Pol3D VRO' 'R Ip VRO' 'R In VRO' 'ATPase2C' ...
    'Pentamer cyt' 'Pentamer VRO' 'RNA-Bound Pentamer' 'P2Filled' 'P3Filled' 'P4Filled' 'P5Filled' 'P6Filled' 'P7Filled' 'P8Filled' 'P9Filled' 'P10Filled' 'P11Filled' 'Virion' ... 
    'P2Empty' 'P3Empty' 'P4Empty' 'P5Empty' 'P6Empty' 'P7Empty' 'P8Empty' 'P9Empty' 'P10Empty' 'P11Empty' 'Empty Provirion' ...
    'Viral Ribosomes' 'Host Ribosomes' 'Protease' 'ISG Protein'};

%Create a labeled table of initial conditions
LabeledInitVals = table(InitVals','RowNames',InitValLabels,'VariableNames',{'Value'});

%% Runs the ODE and Sensitivity analysis and Assigns Outputs

%Additional ODE solver inputs and options
TSPAN = linspace(0,MaxTime,60*MaxTime);  %Range of hours for the simulation
Options = odeset('RelTol',10e-6,'AbsTol',10e-6,'NonNegative',1:length(InitVals));

[t,SolutionsMatrix] = ode15s(@(t,y)CVB3ODEfunc(t,y,Constants,IFNSwitch,VirResponse,IFNStimulation,IFNStimulationTime),TSPAN,InitVals,Options); %Solving the ODE

%Store all results
SolutionsTensor(:,:,i) = SolutionsMatrix(:,:);

%Runs sensitivity analysis if toggled on
if strcmpi(SensitivityAnalysis,'on')
    [SensTable]= SensitivityAnalysisfunc(@CVB3ODEfunc,LabeledConstants,LabeledInitVals,SensitivityAnalysisOutput,OtherSpecies,ScalingFactor,MaxScalingOrder,IFNSwitch,VirResponse,IFNStimulation,IFNStimulationTime,MaxTime,Options);
end

fprintf('%2d ',i) %Display current run count

end

%% Prepares Solutions for export and plotting if desired

%Convert VRO localized concentrations back to cytoplasmic concentrations when VROs are present
for k = 1:size(SolutionsTensor,3)
    i = 1; j = size(SolutionsTensor,1);
    while i < j && SolutionsTensor(i,22,k) < VROFormThreshold * 6e-4 %Pol3D_VRO; 6e-4 is concentration of 1 molecule in nM
        i = i + 1;
    end
    if i ~= j
        SolutionsTensor(i:j,20:50,k) = SolutionsTensor(i:j,20:50,k)./VROVolConv; %Scale VRO-resident species
    end
end

%Permute SolutionsTensor for summary species
PermuteSolutionsTensor = permute(SolutionsTensor,[1,3,2]);

%Total VP1- Note: Surface-bound species not counted because experimentally cells are washed before lysis
TotalP1 = 0*InternalVolConv*60*(sum(PermuteSolutionsTensor(:,:,4:5),3) ... %bDAF, bDAF_Defective
    + sum(PermuteSolutionsTensor(:,:,7:8),3) ... %bDAF_TJ, bDAF_Defective_TJ
    + sum(PermuteSolutionsTensor(:,:,12:13),3)) ... %bCAR, bCAR_Defective
    + 60*sum(PermuteSolutionsTensor(:,:,14:15),3) ...  %R_p_endo, R_p_endo_Defective
    + 5*(sum(PermuteSolutionsTensor(:,:,26:28),3) ...  %Pentamer_cyt, Pentamer_VRO, RNAPentamer
    + 2*sum(PermuteSolutionsTensor(:,:,[29 40]),3) ...  %P2Filled, P2Empty
    + 3*sum(PermuteSolutionsTensor(:,:,[30 41]),3) ...  %P3Filled, P3Empty
    + 4*sum(PermuteSolutionsTensor(:,:,[31 42]),3) ...  %P4Filled, P4Empty
    + 5*sum(PermuteSolutionsTensor(:,:,[32 43]),3) ...  %P5Filled, P5Empty
    + 6*sum(PermuteSolutionsTensor(:,:,[33 44]),3) ...  %P6Filled, P6Empty
    + 7*sum(PermuteSolutionsTensor(:,:,[34 45]),3) ...  %P7Filled, P7Empty
    + 8*sum(PermuteSolutionsTensor(:,:,[35 46]),3) ...  %P8Filled, P8Empty
    + 9*sum(PermuteSolutionsTensor(:,:,[36 47]),3) ...  %P9Filled, P9Empty
    + 10*sum(PermuteSolutionsTensor(:,:,[37 48]),3) ... %P10Filled, P10Empty
    + 11*sum(PermuteSolutionsTensor(:,:,[38 49]),3) ... %P11Filled, P11Empty
    + 12*sum(PermuteSolutionsTensor(:,:,[39 50]),3));   %Virion, EmptyProvirion

%Total viral RdRp 3D polymerase
TotalP2 = sum(PermuteSolutionsTensor(:,:,22:24),3); %Pol3D_VRO, R_Ip_VRO, R_In_VRO;

%Total viral protease
TotalP3 = PermuteSolutionsTensor(:,:,53); %Protease

%Total VRO associated viral proteins
TotalP4 = sum(PermuteSolutionsTensor(:,:,[25 27]),3) ... %ATPase2C, Pentamer_VRO
    + sum(PermuteSolutionsTensor(:,:,40:50),3); %Empty pentamer species; filled pentamers are bound by RNA, not ATPase2C
    
%Creates summary species
TotalPosStrands = InternalVolConv*(0*sum(PermuteSolutionsTensor(:,:,4:5),3) ... %bDAF, bDAF_Defective; assumed instantly removed before lysis
    + 0*sum(PermuteSolutionsTensor(:,:,7:8),3) ... %bDAF_TJ, bDAF_Defective_TJ; assumed instantly removed before lysis
    + sum(PermuteSolutionsTensor(:,:,12:13),3)) ... %bCAR, bCAR_Defective
    + sum(PermuteSolutionsTensor(:,:,14:15),3) ...  %R_p_endo, R_p_endo_Defective
    + sum(PermuteSolutionsTensor(:,:,16:18),3) ...  %R_p_cyt, R_p_cyt_Defective, T_c
    + PermuteSolutionsTensor(:,:,20) ...            %R_p_VRO
    + 0*PermuteSolutionsTensor(:,:,23) ... %R_Ip_VRO
    + sum(PermuteSolutionsTensor(:,:,28:39),3); %RNAPentamer and filled pentamer species
TotalNegStrands = EncapsidatedNegStrandRatio*(InternalVolConv*(0*sum(PermuteSolutionsTensor(:,:,4:5),3) ... %bDAF, bDAF_Defective; assumed instantly removed before lysis
    + 0*sum(PermuteSolutionsTensor(:,:,7:8),3) ... %bDAF_TJ, bDAF_Defective_TJ; assumed instantly removed before lysis
    + sum(PermuteSolutionsTensor(:,:,12:13),3)) ... %bCAR, bCAR_Defective
    + sum(PermuteSolutionsTensor(:,:,14:15),3) ... %R_p_endo, R_p_endo_Defective
    + sum(PermuteSolutionsTensor(:,:,17),3)) ... %R_p_cyt_Defective
    + PermuteSolutionsTensor(:,:,19) ... %R_n_cyt
    + PermuteSolutionsTensor(:,:,21) ... %R_n_VRO
    + 0*PermuteSolutionsTensor(:,:,24); %R_In_VRO
TotalRNARatio = TotalPosStrands./TotalNegStrands;
TotalDoubleStrands = sum(PermuteSolutionsTensor(:,:,23:24),3); %R_Ip_VRO, R_In_VRO;
TotalDoubleStrandsVRO = TotalDoubleStrands.*VROVolConv;

%Create host-virus interaction species for later plotting
if strcmpi(IFNSwitch,'on') && (strcmpi(IFNStimulation,'on') == 0) %If endogenous response is active but there is no exogenous IFN stimulation
    VirDetection = ((TotalDoubleStrandsVRO).^n_Hill)./(kD_Hill^n_Hill + (TotalDoubleStrandsVRO).^n_Hill); %Re-creates the pattern of viral detection in CVB3ODEfunc
elseif strcmpi(IFNStimulation,'on') && strcmpi(IFNSwitch,'on') %If endogenous response is active and the cell is exogenously stimulated at or before infection begins
    VirDetection = ones(length(t),size(PermuteSolutionsTensor,2)); %The Hill equation is at the max value of 1 from the start of infection
    if IFNStimulationTime > 0 %If stimulation occurs after infection has started
        for i = 1:length(t)
            if t(i) >= IFNStimulationTime %Finds the index of the vector at the time when IFN stimulation occurs
                t_stim = i; break
            end
        end
        VirDetection(1:t_stim,:) = ((TotalDoubleStrandsVRO(1:t_stim,:)).^n_Hill)./(kD_Hill^n_Hill + (TotalDoubleStrandsVRO(1:t_stim,:)).^n_Hill); %Set the values before exogenous stimulation to the endogenous response curve
        %Once exogenous stimulation occurs, values are set to 1 (maximum response) rather than endogenous response (already set to 1 above)
    end
elseif strcmpi(IFNSwitch,'on') == 0 %If interferon response is off
    VirDetection = zeros(length(t),size(PermuteSolutionsTensor,2)); %No viral detection occurs in this situation
end
if strcmpi(VirResponse,'on')
    VirDetection = VirDetection.*(1 - Protease./(Protease + EC50_DetectorDeg)); %Scales down viral detection for IFN response if the virus's response to the IFN response is active.
end

%Determine the medians of the model species solutions
MediansMatrix = median(SolutionsTensor,3);
VirDetectionSummary(:,1) = median(VirDetection,2);
TotalPosStrandsSummary(:,1) = median(TotalPosStrands,2);
TotalNegStrandsSummary(:,1) = median(TotalNegStrands,2);
TotalDoubleStrandsSummary(:,1) = median(TotalDoubleStrands,2);
TotalP1Summary(:,1) = median(TotalP1,2);

%Determine the quartiles of the model species solutions
UpperQuantileMatrix = quantile(SolutionsTensor,UpperQuantile,3);
VirDetectionSummary(:,2) = quantile(VirDetection,UpperQuantile,2);
TotalPosStrandsSummary(:,2) = quantile(TotalPosStrands,UpperQuantile,2);
TotalNegStrandsSummary(:,2) = quantile(TotalNegStrands,UpperQuantile,2);
TotalDoubleStrandsSummary(:,2) = quantile(TotalDoubleStrands,UpperQuantile,2);
TotalP1Summary(:,2) = quantile(TotalP1,UpperQuantile,2);

LowerQuantileMatrix = quantile(SolutionsTensor,LowerQuantile,3);
VirDetectionSummary(:,3) = quantile(VirDetection,LowerQuantile,2);
TotalPosStrandsSummary(:,3) = quantile(TotalPosStrands,LowerQuantile,2);
TotalNegStrandsSummary(:,3) = quantile(TotalNegStrands,LowerQuantile,2);
TotalDoubleStrandsSummary(:,3) = quantile(TotalDoubleStrands,LowerQuantile,2);
TotalP1Summary(:,3) = quantile(TotalP1,LowerQuantile,2);

%Determine the means of the model species solutions
MeansMatrix = mean(SolutionsTensor,3);
VirDetectionSummary(:,4) = mean(VirDetection,2);
TotalPosStrandsSummary(:,4) = mean(TotalPosStrands,2);
TotalNegStrandsSummary(:,4) = mean(TotalNegStrands,2);
TotalDoubleStrandsSummary(:,4) = mean(TotalDoubleStrands,2);
TotalP1Summary(:,4) = mean(TotalP1,2);

%Create new solutions tensor
SummarySolutionsTensor(:,:,1) = MediansMatrix;
SummarySolutionsTensor(:,:,2) = UpperQuantileMatrix;
SummarySolutionsTensor(:,:,3) = LowerQuantileMatrix;
SummarySolutionsTensor(:,:,4) = MeansMatrix;

%Rearrange tensor
dimensions = size(SummarySolutionsTensor);
FinalSolutionsTensor = permute(SummarySolutionsTensor,[1,3,2]);

%Assign Solutions matrices
FinalSolutionsCell = mat2cell(FinalSolutionsTensor,dimensions(1),dimensions(3),ones([1,dimensions(2)]));
[uCVB3,uCVB3_Defective,uDAF,bDAF,bDAF_Defective,uDAF_TJ,bDAF_TJ,bDAF_Defective_TJ,uCVB3_TJ,uCVB3_Defective_TJ,uCAR,bCAR,bCAR_Defective,R_p_endo,R_p_endo_Defective, ...
 R_p_cyt,R_p_cyt_Defective,T_c,R_n_cyt,R_p_VRO,R_n_VRO,Pol3D_VRO,R_Ip_VRO,R_In_VRO,ATPase2C, ...
 Pentamer_cyt,Pentamer_VRO,RNAPentamer,P2Filled,P3Filled,P4Filled,P5Filled,P6Filled,P7Filled,P8Filled,P9Filled,P10Filled,P11Filled,Virion, ...
 P2Empty,P3Empty,P4Empty,P5Empty,P6Empty,P7Empty,P8Empty,P9Empty,P10Empty,P11Empty,EmptyProvirion, ...
 RibAvail,RibUnavail,Protease,ISGProtein] = FinalSolutionsCell{:};

%Get the parsed run parameters to output
if strcmpi(SensitivityAnalysis,'on')
    RunParameters = {sprintf('MOI: %d', MOI)...
        sprintf('MaxTime (h): %d', MaxTime)...
        strcat('IFNSwitch: ',IFNSwitch)...
        strcat('VirResponse: ',VirResponse)...
        strcat('IFNStimulation: ',IFNStimulation)...
        sprintf('IFNStimulationTime (h): %d', IFNStimulationTime)...
        strcat('SensitivityAnalysis: ',SensitivityAnalysis)...
        strcat('SensitivityAnalysisOutput: ' ,SensitivityAnalysisOutput)...
        sprintf('ScalingFactor: %d', ScalingFactor)...
        sprintf('MaxScalingOrder: %d',MaxScalingOrder)...
        strcat('ExportData: ',ExportData)...
        strcat('DataFile: ',DataFile)...
        strcat('RunCount',RunCount)...
        strcat('UpperQuantile',UpperQuantile)...
        strcat('LowerQuantile',LowerQuantile)...
        sprintf('CV: %d',CV)};
else
    RunParameters = {sprintf('MOI: %d', MOI)...
        sprintf('MaxTime (h): %d', MaxTime)...
        strcat('IFNSwitch: ',IFNSwitch)...
        strcat('VirResponse: ',VirResponse)...
        strcat('IFNStimulation: ',IFNStimulation)...
        sprintf('IFNStimulationTime (h): %d', IFNStimulationTime)...
        strcat('SensitivityAnalysis: ',SensitivityAnalysis)...
        strcat('ExportData: ',ExportData)...
        strcat('DataFile: ',DataFile)...
        strcat('RunCount',RunCount)...
        strcat('UpperQuantile',UpperQuantile)...
        strcat('LowerQuantile',LowerQuantile)...
        sprintf('CV: %d',CV)};
end

ResultsLabels = {'Time'...
    'Unbound CVB3' 'Defective Unbound CVB3' 'Unbound DAF' 'CVB3-DAF' 'Defective CVB3-DAF' 'Unbound DAF in TJ' 'CVB3-DAF in TJ' 'Defective CVB3-DAF in TJ' 'Unbound CVB3 in TJ' 'Defective Unbound CVB3 in TJ' 'Unbound CAR' 'Bound CAR' 'Defective Bound CAR','Endosomal +ssRNA' 'Defective Endosomal +ssRNA'...
    'Cytoplasmic +ssRNA' 'Defective Cytoplasmic +ssRNA' 'Cytoplasmic -ssRNA' 'Translation Complex' 'Viral Ribosomes' 'Host Ribosomes' 'Protease' 'Cytoplasmic Pentamer'...
    'VRO +ssRNA' 'VRO -ssRNA' 'VRO +ssRNA Replication Complex' 'VRO -ssRNA Replication Complex' 'VRO Viral Polymerase 3D' 'VRO Viral 2C' 'VRO Pentamer'...
    'RNA-Bound Pentamer' '+ssRNA-2xPentamer' '+ssRNA-3xPentamer' '+ssRNA-4xPentamer' '+ssRNA-5xPentamer' '+ssRNA-6xPentamer' '+ssRNA-7xPentamer' '+ssRNA-8xPentamer' '+ssRNA-9xPentamer' '+ssRNA-10xPentamer' '+ssRNA-11xPentamer' 'Virions'...
    '2xPentamer' '3xPentamer' '4xPentamer' '5xPentamer' '6xPentamer' '7xPentamer' '8xPentamer' '9xPentamer' '10xPentamer' '11xPentamer' 'Empty Provirions'...
    'Viral Detection' 'Interferon Stimulated Proteins'...
    'Total +ssRNA' 'Total -ssRNA' 'Total dsRNA' 'Total Viral Protein'};

MedianTable = table(t,...%Time
    uCVB3(:,1),uCVB3_Defective(:,1),uDAF(:,1),bDAF(:,1),bDAF_Defective(:,1),uDAF_TJ(:,1),bDAF_TJ(:,1),bDAF_Defective_TJ(:,1),uCVB3_TJ(:,1),uCVB3_Defective_TJ(:,1),uCAR(:,1),bCAR(:,1),bCAR_Defective(:,1),R_p_endo(:,1),R_p_endo_Defective(:,1),...%Delivery
    R_p_cyt(:,1),R_p_cyt_Defective(:,1),R_n_cyt(:,1),T_c(:,1),RibAvail(:,1),RibUnavail(:,1),Protease(:,1),Pentamer_cyt(:,1),...%Translation and cytoplasmic species
    R_p_VRO(:,1),R_n_VRO(:,1),R_Ip_VRO(:,1),R_In_VRO(:,1),Pol3D_VRO(:,1),ATPase2C(:,1),Pentamer_VRO(:,1),...%Viral Replication Organelle 
    RNAPentamer(:,1),P2Filled(:,1),P3Filled(:,1),P4Filled(:,1),P5Filled(:,1),P6Filled(:,1),P7Filled(:,1),P8Filled(:,1),P9Filled(:,1),P10Filled(:,1),P11Filled(:,1),Virion(:,1),...%Virion Assembly
    P2Empty(:,1),P3Empty(:,1),P4Empty(:,1),P5Empty(:,1),P6Empty(:,1),P7Empty(:,1),P8Empty(:,1),P9Empty(:,1),P10Empty(:,1),P11Empty(:,1),EmptyProvirion(:,1),...
    VirDetectionSummary(:,1),ISGProtein(:,1),...%Host-virus Interaction
    TotalPosStrandsSummary(:,1), TotalNegStrandsSummary(:,1), TotalDoubleStrandsSummary(:,1), TotalP1Summary(:,1),...%Summary species
    'VariableNames',ResultsLabels);

UpperQuantileTable = table(t,...%Time
    uCVB3(:,2),uCVB3_Defective(:,2),uDAF(:,2),bDAF(:,2),bDAF_Defective(:,2),uDAF_TJ(:,2),bDAF_TJ(:,2),bDAF_Defective_TJ(:,2),uCVB3_TJ(:,2),uCVB3_Defective_TJ(:,2),uCAR(:,2),bCAR(:,2),bCAR_Defective(:,2),R_p_endo(:,2),R_p_endo_Defective(:,2),...%Delivery
    R_p_cyt(:,2),R_p_cyt_Defective(:,2),R_n_cyt(:,2),T_c(:,2),RibAvail(:,2),RibUnavail(:,2),Protease(:,2),Pentamer_cyt(:,2),...%Translation and cytoplasmic species
    R_p_VRO(:,2),R_n_VRO(:,2),R_Ip_VRO(:,2),R_In_VRO(:,2),Pol3D_VRO(:,2),ATPase2C(:,2),Pentamer_VRO(:,2),...%Viral Replication Organelle 
    RNAPentamer(:,2),P2Filled(:,2),P3Filled(:,2),P4Filled(:,2),P5Filled(:,2),P6Filled(:,2),P7Filled(:,2),P8Filled(:,2),P9Filled(:,2),P10Filled(:,2),P11Filled(:,2),Virion(:,2),...%Virion Assembly
    P2Empty(:,2),P3Empty(:,2),P4Empty(:,2),P5Empty(:,2),P6Empty(:,2),P7Empty(:,2),P8Empty(:,2),P9Empty(:,2),P10Empty(:,2),P11Empty(:,2),EmptyProvirion(:,2),...
    VirDetectionSummary(:,2),ISGProtein(:,2),...%Host-virus Interaction
    TotalPosStrandsSummary(:,2), TotalNegStrandsSummary(:,2), TotalDoubleStrandsSummary(:,2), TotalP1Summary(:,2),...%Summary species
    'VariableNames',ResultsLabels);

LowerQuantileTable = table(t,...%Time
    uCVB3(:,3),uCVB3_Defective(:,3),uDAF(:,3),bDAF(:,3),bDAF_Defective(:,3),uDAF_TJ(:,3),bDAF_TJ(:,3),bDAF_Defective_TJ(:,3),uCVB3_TJ(:,3),uCVB3_Defective_TJ(:,3),uCAR(:,3),bCAR(:,3),bCAR_Defective(:,3),R_p_endo(:,3),R_p_endo_Defective(:,3),...%Delivery
    R_p_cyt(:,3),R_p_cyt_Defective(:,3),R_n_cyt(:,3),T_c(:,3),RibAvail(:,3),RibUnavail(:,3),Protease(:,3),Pentamer_cyt(:,3),...%Translation and cytoplasmic species
    R_p_VRO(:,3),R_n_VRO(:,3),R_Ip_VRO(:,3),R_In_VRO(:,3),Pol3D_VRO(:,3),ATPase2C(:,3),Pentamer_VRO(:,3),...%Viral Replication Organelle 
    RNAPentamer(:,3),P2Filled(:,3),P3Filled(:,3),P4Filled(:,3),P5Filled(:,3),P6Filled(:,3),P7Filled(:,3),P8Filled(:,3),P9Filled(:,3),P10Filled(:,3),P11Filled(:,3),Virion(:,3),...%Virion Assembly
    P2Empty(:,3),P3Empty(:,3),P4Empty(:,3),P5Empty(:,3),P6Empty(:,3),P7Empty(:,3),P8Empty(:,3),P9Empty(:,3),P10Empty(:,3),P11Empty(:,3),EmptyProvirion(:,3),...
    VirDetectionSummary(:,3),ISGProtein(:,3),...%Host-virus Interaction
    TotalPosStrandsSummary(:,3), TotalNegStrandsSummary(:,3), TotalDoubleStrandsSummary(:,3), TotalP1Summary(:,3),...%Summary species
    'VariableNames',ResultsLabels);

MeanTable = table(t,...%Time
    uCVB3(:,4),uCVB3_Defective(:,4),uDAF(:,4),bDAF(:,4),bDAF_Defective(:,4),uDAF_TJ(:,4),bDAF_TJ(:,4),bDAF_Defective_TJ(:,4),uCVB3_TJ(:,4),uCVB3_Defective_TJ(:,4),uCAR(:,4),bCAR(:,4),bCAR_Defective(:,4),R_p_endo(:,4),R_p_endo_Defective(:,4),...%Delivery
    R_p_cyt(:,4),R_p_cyt_Defective(:,4),R_n_cyt(:,4),T_c(:,4),RibAvail(:,4),RibUnavail(:,4),Protease(:,4),Pentamer_cyt(:,4),...%Translation and cytoplasmic species
    R_p_VRO(:,4),R_n_VRO(:,4),R_Ip_VRO(:,4),R_In_VRO(:,4),Pol3D_VRO(:,4),ATPase2C(:,4),Pentamer_VRO(:,4),...%Viral Replication Organelle 
    RNAPentamer(:,4),P2Filled(:,4),P3Filled(:,4),P4Filled(:,4),P5Filled(:,4),P6Filled(:,4),P7Filled(:,4),P8Filled(:,4),P9Filled(:,4),P10Filled(:,4),P11Filled(:,4),Virion(:,4),...%Virion Assembly
    P2Empty(:,4),P3Empty(:,4),P4Empty(:,4),P5Empty(:,4),P6Empty(:,4),P7Empty(:,4),P8Empty(:,4),P9Empty(:,4),P10Empty(:,4),P11Empty(:,4),EmptyProvirion(:,4),...
    VirDetectionSummary(:,4),ISGProtein(:,4),...%Host-virus Interaction
    TotalPosStrandsSummary(:,4), TotalNegStrandsSummary(:,4), TotalDoubleStrandsSummary(:,4), TotalP1Summary(:,4),...%Summary species
    'VariableNames',ResultsLabels);

%MatLab outputs tables to the floating point precision of the workspace so the results need to be rounded back to the precision of the odesolver.
MedianTable.Variables = round(MedianTable.Variables, 6);
UpperQuantileTable.Variables = round(UpperQuantileTable.Variables, 6);
LowerQuantileTable.Variables = round(LowerQuantileTable.Variables, 6);
MeanTable.Variables = round(MeanTable.Variables, 6);

%Assign ResultsTable variable units
unitArray = cell(1,59);
unitArray(:) = {'nM'};
MedianTable.Properties.VariableUnits = cat(2,{'Hours'},unitArray);
UpperQuantileTable.Properties.VariableUnits = cat(2,{'Hours'},unitArray);
LowerQuantileTable.Properties.VariableUnits = cat(2,{'Hours'},unitArray);
MeanTable.Properties.VariableUnits = cat(2,{'Hours'},unitArray);

%Assign run parameters to table description
MedianTable.Properties.Description = strjoin(string(RunParameters),', ');
UpperQuantileTable.Properties.Description = strjoin(string(RunParameters),', ');
LowerQuantileTable.Properties.Description = strjoin(string(RunParameters),', ');
SensTable.Properties.Description = strjoin(string(RunParameters),', ');
MeanTable.Properties.Description = strjoin(string(RunParameters),', ');

%Assign the output
if strcmpi(SensitivityAnalysis,'on')
    SensitivitySolutions = SensTable;
else
    SensitivitySolutions = [];
end

MedianSolutions = MedianTable;
UpperQuantileSolutions = UpperQuantileTable;
LowerQuantileSolutions = LowerQuantileTable;
MeanSolutions = MeanTable;

%Export processed data if desired
if strcmpi(ExportData,'on')
    if strcmpi(SensitivityAnalysis,'on')
       writecell(RunParameters ,strjoin({'RunParameters',DataFile},{'_'}));
       writetable(SensTable,DataFile,'WriteVariableNames',true,'WriteRowNames',true);
       writetable(MedianTable,DataFile,'WriteVariableNames',true);
       writetable(UpperQuantileTable,DataFile,'WriteVariableNames',true);
       writetable(LowerQuantileTable,DataFile,'WriteVariableNames',true);
       writetable(MeanTable,DataFile,'WriteVariableNames',true);       
    else
       writecell(RunParameters ,strjoin({'RunParameters',DataFile},{'_'}));
       writetable(MedianTable,strjoin({'MedianResults',DataFile},{'_'}),'WriteVariableNames',true);
       writetable(UpperQuantileTable,strjoin({'UpperQuantileResults',DataFile},{'_'}),'WriteVariableNames',true);
       writetable(LowerQuantileTable,strjoin({'LowerQuantileResults',DataFile},{'_'}),'WriteVariableNames',true);
       writetable(MeanTable,strjoin({'MeanResults',DataFile},{'_'}),'WriteVariableNames',true);
    end
end

%Plot processed data
if strcmpi(PlotResults,'on')
   CVB3ODEPlots(MedianTable,UpperQuantileTable,LowerQuantileTable,MeanTable,'PlotType',plotType,'PopulationSetting',PopulationSetting);
end

%% Outputs run summary in command window

%Displays molecule counts of key species
fprintf('\nPositive Strands = %3.2e molecules, %3.2e nM\n',TotalPosStrands(end)* 6.022*10^23 * CellCytVolume / 10^24,TotalPosStrands(end))
fprintf('Negative Strands = %3.2e molecules, % 3.2e nM\n',TotalNegStrands(end)* 6.022*10^23 * CellCytVolume / 10^24,TotalNegStrands(end))
fprintf('Double Strands = %3.2e molecules, % 3.2e nM\n',TotalDoubleStrandsVRO(end)* 6.022*10^23 * CellCytVolume / 10^24,TotalDoubleStrandsVRO(end))
fprintf('RNA Ratio = %3.2e\n',TotalRNARatio(end));
fprintf('Polyproteins from: Capsid = %3.2e nM, Polymerases = %3.2e nM, Proteases = %3.2e nM, ATPases = %3.2e nM\n',TotalP1(end),TotalP2(end),TotalP3(end),TotalP4(end))
fprintf('Virions = %3.2e molecules, % 3.2e nM\n',Virion(end)* 6.022*10^23 * CellCytVolume / 10^24,Virion(end))
fprintf('Empty Capsid = %3.2e molecules, % 3.2e nM\n',EmptyProvirion(end)* 6.022*10^23 * CellCytVolume / 10^24,EmptyProvirion(end))

end