function [Results,Inputs,InputLabels] = IFNMechanismInteractions(PrimeTime,EC50_RNAdeg,EC50_DetectorDeg,EC50_Protease,EC50_Translate,Output,OutputTime,MOI);
%Function: Assesses the model results of simultaneously altering multiple IFN mechanism parameters
%
%INPUTS:
    %PrimeTime: Time of IFN stimulation. Negative in the case of pre-stimulation.
    %EC50_RNAdeg: Base Value = 5 nM. Set 50% effect at 5 nM (https://www.ncbi.nlm.nih.gov/pubmed/17344297)
    %EC50_DetectorDeg: Base Value = 1 pM. Arbitrary value of 1 pM right now
    %EC50_Protease: Base value = 20 nM. Half-effective concentration of ISG protein for inhibition of 2Apro by ISG15. Currently arbitrary. See: https://www.ncbi.nlm.nih.gov/pubmed/25165091
    %EC50_Translate: Base Value = 10 nM. Half-effective concentration of ISG protein for translation inhibition. Currently arbitrary.
    %Output: 'Plus RNA','Minus RNA','dsRNA','Polyprotein','Virion', or 'Empty Provirion' to select a species to study. 'Other' prompts input for a specific parameter name.
    %OutputTime: Time post-infection we want to study
    %MOI: Multiplicity  of infection
%
%OUTPUTS:
    %Results: 6D tensor of all results
    %Inputs: Matrix of all input values
    %InputLabels: Vector of strings corresponging to the inputs
    
%IMPORTANT NOTE: Input vectors should have a median = base value. If even in length, the base value should be the lower index.

tTotalStart = tic; %Starts a timer

%Prompts the user to enter additional information for custom analysis    
if strcmpi(Output,'Other')
    fprintf(['The options are: uCVB3, Defective uCVB3, uDAF, bDAF, Defective bDAF, uDAF TJ, bDAF TJ, Defective bDAF TJ, uCVB3 TJ, Defective uCVB3 TJ, uCAR, bCAR, Defective bCAR, R p endo, Defective R p endo,\n' ...
    'R p cyt, Defective R p cyt, T c, R n cyt, VRO Avail, R p VRO, R n VRO, Pol3D VRO, R Ip VRO, R In VRO, ATPase2C,\n' ...
    'Pentamer cyt, Pentamer VRO, RNA-Bound Pentamer, P2Filled, P3Filled, P4Filled, P5Filled, P6Filled, P7Filled, P8Filled, P9Filled, P10Filled, P11Filled, Virion,\n' ...
    'P2Empty, P3Empty, P4Empty, P5Empty, P6Empty, P7Empty, P8Empty, P9Empty, P10Empty, P11Empty, Empty Provirion,\n'
    'Viral Ribosomes, Host Ribosomes, Protease, ISG Protein\n'])
    OtherSpecies = input('Enter the name(s) of the species you wish to study as a vector of their names as strings separated by spaces.');   
else
    OtherSpecies = [];
end
    
%Pre-defining quick options for OutputIndex:
Positive_Strands = {'R p endo' 'Defective R p endo' 'R p cyt' 'Defective R p cyt' 'T c' 'R p VRO' 'RNA-Bound Pentamer' 'Virion' 'P2Filled' 'P3Filled' 'P4Filled' 'P5Filled' 'P6Filled' 'P7Filled' 'P8Filled' 'P9Filled' 'P10Filled' 'P11Filled'}; %'R Ip VRO', Previously Didn't Have P__Filled
Negative_Strands = {'R n cyt' 'R n VRO'}; %'R In VRO'
dsRNA = {'R Ip VRO' 'R In VRO'};
Polyproteins = {'Protease'};
Virion = {'Virion'};
EmptyProvirion = {'Empty Provirion'};

%Assigns OutputIndex based on the SensAnalysis input. This is a vector of indices indicating the desired species of study.
if strcmpi(Output,'Plus RNA')
    OutputIndex = Positive_Strands;
elseif strcmpi(Output,'Minus RNA')
    OutputIndex = Negative_Strands;
elseif strcmpi(Output,'dsRNA')
    OutputIndex = dsRNA;
elseif strcmpi(Output,'Polyprotein')
    OutputIndex = Polyproteins;
elseif strcmpi(Output,'Virion')
    OutputIndex = Virion;    
elseif strcmpi(Output,'Empty Provirion')
    OutputIndex = EmptyProvirion;
elseif strcmpi(Output,'Other')
    OutputIndex = OtherSpecies;
else
    fprintf('\nError: Invalid Output Input.\n')
    return
end

%Creates labels vector for the above initial values. If label is changed, ensure it is changed in SensitivityAnalysis and CellBurst as well.
InitValLabels = {'uCVB3' 'Defective uCVB3' 'uDAF' 'bDAF' 'Defective bDAF' 'uDAF TJ' 'bDAF TJ' 'Defective bDAF TJ' 'uCVB3 TJ' 'Defective uCVB3 TJ' 'uCAR' 'bCAR' 'Defective bCAR' 'R p endo' 'Defective R p endo' ...
    'R p cyt' 'Defective R p cyt' 'T c' 'R n cyt' 'R p VRO' 'R n VRO' 'Pol3D VRO' 'R Ip VRO' 'R In VRO' 'ATPase2C' ...
    'Pentamer cyt' 'Pentamer VRO' 'RNA-Bound Pentamer' 'P2Filled' 'P3Filled' 'P4Filled' 'P5Filled' 'P6Filled' 'P7Filled' 'P8Filled' 'P9Filled' 'P10Filled' 'P11Filled' 'Virion' ... 
    'P2Empty' 'P3Empty' 'P4Empty' 'P5Empty' 'P6Empty' 'P7Empty' 'P8Empty' 'P9Empty' 'P10Empty' 'P11Empty' 'Empty Provirion' ...
    'Viral Ribosomes' 'Host Ribosomes' 'Protease' 'ISG Protein'};

%Initializing bookkeeping values, Priming inputs, and results array
Priming = [1 0]; %1 is 'on', 0 is 'off'
Value = 0;
InputSize = length(EC50_Translate);
Results = zeros(length(Priming),length(PrimeTime),length(EC50_RNAdeg),length(EC50_DetectorDeg),length(EC50_Protease),length(EC50_Translate));
Inputs = zeros(6,InputSize);

for a = 1:length(Priming)
    Inputs(1,a) = Priming(a);
    for b = 1:length(PrimeTime)
        Inputs(2,b) = PrimeTime(b);
        for c = 1:length(EC50_RNAdeg)
            Inputs(3,c) = EC50_RNAdeg(c);
            for d = 1:length(EC50_DetectorDeg)
                Inputs(4,d) = EC50_DetectorDeg(d);
                for e = 1:length(EC50_Protease)
                    Inputs(5,e) = EC50_Protease(e);
                    for f = 1:length(EC50_Translate)
                        
                        %get Priming in correct format
                        if Priming(a) == 1
                            PrimingIn='on';
                        else
                            PrimingIn='off';
                        end

                        [Solutions] = CVB3ODEEval(MOI,'MaxTime',OutputTime,'IFNSwitch','on','VirResponse','on','Priming',PrimingIn,'PrimeTime',PrimeTime(b),'EC50_RNAdeg',EC50_RNAdeg(c),'EC50_DetectorDeg',EC50_DetectorDeg(d),'EC50_Protease',EC50_Protease(e),'EC50_Translate',EC50_Translate(f),'plotsOn','off','sensOn','off');
                        
                        LabeledSolutions = table(Solutions(end,:)','RowNames',InitValLabels,'VariableNames',{'Value'}); %Adding labels to the solutions matrix. Only uses end values.
                        
                        for g = 1:length(OutputIndex)
                            Value = Value + LabeledSolutions{OutputIndex(g),:};
                        end
                        
                        Results(a,b,c,d,e,f) = Value;
                        Inputs(6,f) = EC50_Translate(f);
                        Value = 0; 
                                
                    end
                end
            end
        end
    end
end

%Places NaN in extra slots of the inputs vector (Specifically for Priming)
if InputSize > 2
    for i = 3:InputSize
        Inputs(1,i) = NaN;
    end
end

%Vector for the heatmap color bar
Colors = [.3137 .3647 .4549; .3333 .4706 .6235; .3922 .5686 .7098; .4902 .6353 .7686; .6157 .7255 .8314; ...
    .9647 .9647 .9686; .9412 .7059 .5373; .8902 .6157 .4667; .8353 .5333 .4039; .7176 .3843 .3294; .5020 .3373 .3255];

InputLabels = ['Priming' 'Priming Time' 'EC50 RNAdeg' 'EC50 DetectorDeg' 'EC50 Protease' 'EC50 Translate'];

%% Plots

%Determines the median index for the base value of the model
if mod(InputSize,2) == 0 %If input vectors are even in length
    BaseValIndex = InputSize/2;
else %If input vectors are odd in length
    BaseValIndex = InputSize/2 + 0.5;
end

figure
%With IFN priming ON
subplot(5,5,2)
h=heatmap(flipud(reshape(Results(1,:,:,BaseValIndex,BaseValIndex,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(3,:); xlabel('EC50 RNAdeg (nM)')
h.YData=fliplr(Inputs(2,:)); ylabel('Priming Time (hpi)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,3)
h=heatmap(flipud(reshape(Results(1,:,BaseValIndex,:,BaseValIndex,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(4,:).*1000; xlabel('EC50 DetectorDeg (pM)')
h.YData=fliplr(Inputs(2,:)); ylabel('Priming Time (hpi)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,4)
h=heatmap(flipud(reshape(Results(1,:,BaseValIndex,BaseValIndex,:,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(5,:); xlabel('EC50 Protease (nM)')
h.YData=fliplr(Inputs(2,:)); ylabel('Priming Time (hpi)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,5)
h=heatmap(flipud(reshape(Results(1,:,BaseValIndex,BaseValIndex,BaseValIndex,:),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(6,:); xlabel('EC50 Translate')
h.YData=fliplr(Inputs(2,:)); ylabel('Priming Time (hpi)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,8)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,:,:,BaseValIndex,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(4,:).*1000; xlabel('EC50 DetectorDeg (pM)')
h.YData=fliplr(Inputs(3,:)); ylabel('EC50 RNAdeg (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,9)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,:,BaseValIndex,:,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(5,:); xlabel('EC50 Protease (nM)')
h.YData=fliplr(Inputs(3,:)); ylabel('EC50 RNAdeg (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,10)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,:,BaseValIndex,BaseValIndex,:),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(6,:); xlabel('EC50 Translate (nM)')
h.YData=fliplr(Inputs(3,:)); ylabel('EC50 RNAdeg (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,14)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,BaseValIndex,:,:,BaseValIndex),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(5,:); xlabel('EC50 Protease (nM)')
h.YData=fliplr(Inputs(4,:)).*1000; ylabel('EC50 DetectorDeg (pM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,15)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,BaseValIndex,:,BaseValIndex,:),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(6,:); xlabel('EC50 Translate (nM)')
h.YData=fliplr(Inputs(4,:)).*1000; ylabel('EC50 DetectorDeg (pM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,20)
h=heatmap(flipud(reshape(Results(1,BaseValIndex,BaseValIndex,BaseValIndex,:,:),InputSize,InputSize)),'Colormap',Colors);
h.XData=Inputs(6,:); xlabel('EC50 Translate (nM)')
h.YData=fliplr(Inputs(5,:)); ylabel('EC50 Protease (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

%With IFN priming OFF
%Results are transposed so axes can flip (Any time index(XData) < index(Ydata) in Inputs)
subplot(5,5,6)
h=heatmap(flipud(reshape(Results(2,:,:,BaseValIndex,BaseValIndex,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(2,:); xlabel('Priming Time (hpi)')
h.YData=fliplr(Inputs(3,:)); ylabel('EC50 RNAdeg (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,11)
h=heatmap(flipud(reshape(Results(2,:,BaseValIndex,:,BaseValIndex,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(2,:); xlabel('Priming Time (hpi)')
h.YData=fliplr(Inputs(4,:)).*1000; ylabel('EC50 DetectorDeg (pM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,12)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,:,:,BaseValIndex,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(3,:); xlabel('EC50 RNAdeg (nM)')
h.YData=fliplr(Inputs(4,:)).*1000; ylabel('EC50 DetectorDeg (pM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,16)
h=heatmap(flipud(reshape(Results(2,:,BaseValIndex,BaseValIndex,:,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(2,:); xlabel('Priming Time (hpi)')
h.YData=fliplr(Inputs(5,:)); ylabel('EC50 Protease (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,17)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,:,BaseValIndex,:,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(3,:); xlabel('EC50 RNAdeg (nM)')
h.YData=fliplr(Inputs(5,:)); ylabel('EC50 Protease (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,18)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,BaseValIndex,:,:,BaseValIndex),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(4,:).*1000; xlabel('EC50 DetectorDeg (pM)')
h.YData=fliplr(Inputs(5,:)); ylabel('EC50 Protease (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,21)
h=heatmap(flipud(reshape(Results(2,:,BaseValIndex,BaseValIndex,BaseValIndex,:),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(2,:); xlabel('Priming Time (hpi)')
h.YData=fliplr(Inputs(6,:)); ylabel('EC50 Translate (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,22)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,:,BaseValIndex,BaseValIndex,:),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(3,:); xlabel('EC50 RNAdeg (nM)')
h.YData=fliplr(Inputs(6,:)); ylabel('EC50 Translate (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,23)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,BaseValIndex,:,BaseValIndex,:),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(4,:).*1000; xlabel('EC50 DetectorDeg (pM)')
h.YData=fliplr(Inputs(6,:)); ylabel('EC50 Translate (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,24)
h=heatmap(flipud(reshape(Results(2,BaseValIndex,BaseValIndex,BaseValIndex,:,:),InputSize,InputSize)'),'Colormap',Colors);
h.XData=Inputs(5,:); xlabel('EC50 Protease (nM)')
h.YData=fliplr(Inputs(6,:)); ylabel('EC50 Translate (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

%Diagonal of subplot- Comparing with and without IFN priming
subplot(5,5,1)
h=heatmap(reshape(Results(:,:,BaseValIndex,BaseValIndex,BaseValIndex,BaseValIndex),2,InputSize),'Colormap',Colors);
h.YData=Inputs(1,1:2); ylabel('Priming (On/Off)')
h.XData=Inputs(2,:); xlabel('Priming Time (hpi)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,7)
h=heatmap(reshape(Results(:,:,BaseValIndex,BaseValIndex,BaseValIndex,BaseValIndex),2,InputSize),'Colormap',Colors);
h.YData=Inputs(1,1:2); ylabel('Priming (On/Off)')
h.XData=Inputs(3,:); xlabel('EC50 RNAdeg (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,13)
h=heatmap(reshape(Results(:,:,BaseValIndex,BaseValIndex,BaseValIndex,BaseValIndex),2,InputSize),'Colormap',Colors);
h.YData=Inputs(1,1:2); ylabel('Priming (On/Off)')
h.XData=Inputs(4,:).*1000; xlabel('EC50 DetectorDeg (pM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,19)
h=heatmap(reshape(Results(:,:,BaseValIndex,BaseValIndex,BaseValIndex,BaseValIndex),2,InputSize),'Colormap',Colors);
h.YData=Inputs(1,1:2); ylabel('Priming (On/Off)')
h.XData=Inputs(5,:); xlabel('EC50 Protease (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

subplot(5,5,25)
h=heatmap(reshape(Results(:,:,BaseValIndex,BaseValIndex,BaseValIndex,BaseValIndex),2,InputSize),'Colormap',Colors);
h.YData=Inputs(1,1:2); ylabel('Priming (On/Off)')
h.XData=Inputs(6,:); xlabel('EC50 Translate (nM)')
h.ColorLimits=[min(min(min(min(min(min(Results)))))) max(max(max(max(max(max(Results))))))];

Runtime = toc(tTotalStart)/60; %Ends timer, stores in minutes
fprintf('\nThe total required time was %4f minutes\n\n',Runtime)
%sgtitle('') This will be the overall title for the subplot, but this function does not exist in my version of Matlab. Include MOI, time being studied, and output type in the title

%% Notes

%caxis(round([-log2(caxislimit) log2(caxislimit)],0))- Set at every map or have 1 unified color bar

%1. Only check "PrimeTime" changes if Priming is on? Do we want to? Would move it to the last of the 6 and use a conditional break
%2. How to suppress output from each run?

% %Figure out color bars and normalization:
% Color Bars:
% Log or linear scaling?
% Unified color bar for every heatmap or keep them separate?
% Normalization:
% Should we just normalize everything to the output at all base values? If so:
% Should there be a separate normalized value for ?on? and ?off? priming? And what PrimeTime value is the base value?

% Add overall title for subplots that will include the MOI, output species, and time being studied.
% This function exists in a newer version of Matlab but does not exist for mine.

% Find the best way to note on the plot that the lower triangle has priming off and the upper triangle has priming on
