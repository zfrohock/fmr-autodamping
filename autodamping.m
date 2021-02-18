function [dataArray,Hpp,deltaH0,damping,gamma,Ms] = autodamping(filename,frequency,sampleName,varargin)
% Automatically calculates Gilbert damping from BBFMR data. WARNING : This function will close all open figures!
%% INPUTS
% filename = column array of filenames (strings) corresponding to the FMR data to
%            be analyzed. Input in order of increasing frequency
% frequency = column array of frequencies (in GHz) associated with the
%            respective data file referenced by 'filename'
% sampleName = string consisting of the name of your sample. Avoid using
%             special characters or place a / before them to avoid any
%              strange characters appearing in graph titles. 
% varargin = key-value pairs as defined below.
%           key = 'geometry' : values 'IP', 'OOP'
%              What geometry is the sample in. In plane or out of plane.
%           key = 'AutomateGaussianSelection' : values 'True,' 'False'
%              Allow the algorithm to select the data region without user
%              input
%           key = 'FitType' : values 'GUI,' 'Lorentzian,' 'Gaussian'
%              Dictate how fitting takes place. Via a graphical interface
%              with options or simply using a given form.
%           key = 'SuppressFigures' : values 'True,' 'False'
%              Choose to have various figures pop up throughout the fitting
%              process or simply have them saved for future viewing
%           key = 'DeviceType' : values 'BBFMR,' 'VNAFMR'
%              In development. Use only for BBFMR system
%           key = 'AutomateFitting' : values 'True,' 'False'
%              Currently unsused.
%% OUTPUTS
% dataArray = Structure of all data
% Hpp = Peak-to-peak linewidths
% deltaH0 = x intercept of the damping fit
% damping = calculated Gilbert damping coefficient
% gamma = calculated absolute gyromagnetic ratio
% Ms = saturation magnetization
% MATLAB figure files = Saves figures of the damping graph, all Lorentzian
%                       fits, and the gamma graph

%% MAIN FUNCTION BODY
% SYSTEM INFORMATION
set(0,'units','pixels');
Pix_SS = get(0,'screensize');

% Input name-value pair parsing
parser = inputParser;
% Provide defaults
defaultGeometry = 'IP';
defaultAutomateGaussianSelection = 'False';
defaultDeviceType = 'BBFMR';
defaultFitType = 'GUI';
defaultAutomateFitting = 'False';
defaultSuppressFigures = 'False';
% Validate responses
validGeometry = {'IP','OOP'};
validAutomateGaussianSelection = {'False','True'};
validDeviceType = {'BBFMR','VNAFMR'};
validFitType = {'GUI','Lorentzian','Gaussian'};
validAutomateFitting = {'False','True'};
validSuppressFigures = {'False','True'};
% Checking functions
checkGeometry = @(x) any(validatestring(x,validGeometry));
checkAutomateGaussianSelection = @(x) any(validatestring(x,validAutomateGaussianSelection));
checkDeviceType = @(x) any(validatestring(x,validDeviceType));
checkFitType = @(x) any(validatestring(x,validFitType));
checkAutomateFitting = @(x) any(validatestring(x,validAutomateFitting));
checkSuppressFigures = @(x) any(validatestring(x,validSuppressFigures));
% Build parser
addRequired(parser,'filenames')
addRequired(parser,'frequency')
addRequired(parser,'sampleName',@ischar)
addParameter(parser,'geometry',defaultGeometry,checkGeometry);
addParameter(parser,'AutomateGaussianSelection',defaultAutomateGaussianSelection,checkAutomateGaussianSelection);
addParameter(parser,'DeviceType',defaultDeviceType,checkDeviceType);
addParameter(parser,'FitType',defaultFitType,checkFitType);
addParameter(parser,'AutomateFitting',defaultAutomateFitting,checkAutomateFitting);
addParameter(parser,'SuppressFigures',defaultSuppressFigures,checkSuppressFigures);
% Parse input and give everything helpful names
parse(parser,filename,frequency,sampleName,varargin{:});
filename = parser.Results.filenames;
frequency = parser.Results.frequency;
geometry = parser.Results.geometry;
sampleName = parser.Results.sampleName;
deviceType = parser.Results.DeviceType;
fitType = parser.Results.FitType;
automateFitting = parser.Results.AutomateFitting;
suppressFigures = parser.Results.SuppressFigures;
AutomateGaussianSelection = parser.Results.AutomateGaussianSelection;

% Perform different fitting based upon the device. Each has its own
% separate code block
switch deviceType
    case 'BBFMR'
        % DAMPING CALCULATION
        % Get a handle on the shape of the incoming data
        dataDimension = size(filename);
        numFiles = dataDimension(1);
        namelist = [];
        dataArray = [];
        for i = 1:numFiles
            name = strcat('file',num2str(i));
            % Build a list of filenames 
            namelist{i} = name;
            % Load field and data
            [dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data'))] = LoadFMRData(filename{i,1:end});
            % Find beginning and end of Lorentzian, either via kmeans clustering or
            % gui input
            switch AutomateGaussianSelection
                case 'False'
                    [dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField'))] = FindLorentzian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')));
                case 'True'
                    [dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField'))] = AutomaticFindLorentzian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')));
            end
        end
        % Change how fitting takes place based on user input
        switch fitType
            % Walk the user through visual steps to perform fitting. Best
            % used the data is largely unknown. Slow and manual compared to
            % automatic fitting with suppressed figures.
            case 'GUI'
                % Do fitting on the first file to initialize Hpp
                [HppLorentz(1),fitData.(strcat(namelist{1},'fitLorentz')),fitData.(strcat(namelist{1},'fieldLorentz')),fitData.(strcat(namelist{1},'dataLorentz'))] = FitLorentzian(dataArray.(strcat(namelist{1},'field')),dataArray.(strcat(namelist{1},'data')),dataArray.(strcat(namelist{1},'startField')),dataArray.(strcat(namelist{1},'endField')));
                [HppGauss(1),fitData.(strcat(namelist{1},'fitGauss')),fitData.(strcat(namelist{1},'fieldGauss')),fitData.(strcat(namelist{1},'dataGauss'))] = FitGaussian(dataArray.(strcat(namelist{1},'field')),dataArray.(strcat(namelist{1},'data')),dataArray.(strcat(namelist{1},'startField')),dataArray.(strcat(namelist{1},'endField')));
                % Make a plot of Lorentzian and Gaussian fits of the fist data file. Have
                % the user select which fit type is superior
                figure('Renderer', 'painters', 'Position', [(1/15)*Pix_SS(3) .25*Pix_SS(4) (13/15)*Pix_SS(3) .5*Pix_SS(4)])
                subplot(1,2,1)
                % On the left hand side, the default Lorentzian fit
                plot(fitData.(strcat(namelist{1},'fitLorentz')),fitData.(strcat(namelist{1},'fieldLorentz')),fitData.(strcat(namelist{1},'dataLorentz')))
                xlabel('FMR Field (Oe)')
                ylabel('FMR Power Derivative (arb)') 
                title('Lorentzian Fitting of FMR Power Derivative')
                subplot(1,2,2)
                % On the right hand side the optional Gaussian fit
                plot(fitData.(strcat(namelist{1},'fitGauss')),fitData.(strcat(namelist{1},'fieldGauss')),fitData.(strcat(namelist{1},'dataGauss')))
                xlabel('FMR Field (Oe)')
                ylabel('FMR Power Derivative (arb)') 
                title('Gaussian Fitting of FMR Power Derivative')
                % Ask user which fit is better
                answer = questdlg('Are either the Gaussian or Lorentzian Fits Acceptable?',...
                                    'Fit Selection',...
                                    'Lorentzian','Gaussian','Neither','Lorentzian');
                % Go through relevant cases
                switch answer
                    case 'Lorentzian'
                        % Fit all files to Lorentzians
                        close all
                        Hpp(1) = HppLorentz(1);
                        fitData.(strcat(namelist{1},'fit'))= fitData.(strcat(namelist{1},'fitLorentz'));
                        fitData.(strcat(namelist{1},'field')) = fitData.(strcat(namelist{1},'fieldLorentz'));
                        fitData.(strcat(namelist{1},'data')) = fitData.(strcat(namelist{1},'dataLorentz'));
                        % Loop through fitting of all previous files using previous Hpp as a
                        % lowerbound
                        for i = 2:numFiles
                            [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = FitLorentzian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')),Hpp(i-1));
                        end
                    case 'Gaussian'
                        % Fit all files to Gaussians
                        close all
                        Hpp(1) = HppGauss(1);
                        fitData.(strcat(namelist{1},'fit'))= fitData.(strcat(namelist{1},'fitGauss'));
                        fitData.(strcat(namelist{1},'field')) = fitData.(strcat(namelist{1},'fieldGauss'));
                        fitData.(strcat(namelist{1},'data')) = fitData.(strcat(namelist{1},'dataGauss'));
                        % Loop through fitting of all previous files using previous Hpp as a
                        % lowerbound
                        for i = 2:numFiles
                             [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = FitGaussian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')),Hpp(i-1));
                        end
                    case 'Neither'
                        % Go into interactive fit improvement
                        answer_interactive = questdlg('Which fit was better, Lorentzian or Gaussian?',...
                                                        'Fit Selection',...
                                                        'Lorentzian','Gaussian','Cancel','Lorentzian');
                        switch answer_interactive
                            case 'Lorentzian'
                                close all
                                % Manipulate the bounds and tolerances of the lorentzian
                                % fit file
                                for i = 1:numFiles
                                    [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = InteractiveFitLorentzian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')));
                                end
                            case 'Gaussian'
                                close all
                                % Manipulate the bounds and uncertainty of the gaussian fit
                                % file
                                 for i = 1:numFiles
                                    [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = InteractiveFitGaussian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')));
                                end
                            case 'Cancel'
                                close all
                                error('Must select a fit type for interactive improvement')
                        end

                end
            case 'Lorentzian'
                [HppLorentz(1),fitData.(strcat(namelist{1},'fitLorentz')),fitData.(strcat(namelist{1},'fieldLorentz')),fitData.(strcat(namelist{1},'dataLorentz'))] = FitLorentzian(dataArray.(strcat(namelist{1},'field')),dataArray.(strcat(namelist{1},'data')),dataArray.(strcat(namelist{1},'startField')),dataArray.(strcat(namelist{1},'endField')));
                Hpp(1) = HppLorentz(1);
                        fitData.(strcat(namelist{1},'fit'))= fitData.(strcat(namelist{1},'fitLorentz'));
                        fitData.(strcat(namelist{1},'field')) = fitData.(strcat(namelist{1},'fieldLorentz'));
                        fitData.(strcat(namelist{1},'data')) = fitData.(strcat(namelist{1},'dataLorentz'));
                        % Loop through fitting of all previous files using previous Hpp as a
                        % lowerbound
                        for i = 2:numFiles
                            [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = FitLorentzian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')),Hpp(i-1));
                        end
            case 'Gaussian'
                 [HppGauss(1),fitData.(strcat(namelist{1},'fitGauss')),fitData.(strcat(namelist{1},'fieldGauss')),fitData.(strcat(namelist{1},'dataGauss'))] = FitGaussian(dataArray.(strcat(namelist{1},'field')),dataArray.(strcat(namelist{1},'data')),dataArray.(strcat(namelist{1},'startField')),dataArray.(strcat(namelist{1},'endField')));
                 Hpp(1) = HppGauss(1);
                        fitData.(strcat(namelist{1},'fit'))= fitData.(strcat(namelist{1},'fitGauss'));
                        fitData.(strcat(namelist{1},'field')) = fitData.(strcat(namelist{1},'fieldGauss'));
                        fitData.(strcat(namelist{1},'data')) = fitData.(strcat(namelist{1},'dataGauss'));
                        % Loop through fitting of all previous files using previous Hpp as a
                        % lowerbound
                        for i = 2:numFiles
                             [Hpp(i),fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data'))] = FitGaussian(dataArray.(strcat(namelist{i},'field')),dataArray.(strcat(namelist{i},'data')),dataArray.(strcat(namelist{i},'startField')),dataArray.(strcat(namelist{i},'endField')),Hpp(i-1));
                        end
        end
        % Find error bounds for linewidth so we can make error bars on a plot
        for i = 1:numFiles
            tempConfStore = confint(fitData.(strcat(namelist{i},'fit')),0.95);
            HppErrorLowerbound(i) = tempConfStore(1,4);
            HppErrorUpperbound(i) = tempConfStore(2,4);
            HfmrErrorLowerbound(i) = tempConfStore(1,5);
            HfmrErrorUpperbound(i) = tempConfStore(2,5);
            Hpp_Error(i) = abs(HppErrorUpperbound(i) - HppErrorLowerbound(i))/2;
            Hfmr_Error(i) = abs(HfmrErrorUpperbound(i) - HfmrErrorLowerbound(i))/2;
        end
        % Calculate Damping
        [damping,dampingVariation,alphaFitData,deltaH0] = CalculateDamping(Hpp,frequency);
        % Plot the linewidth data points and fit that determines damping
        switch suppressFigures
            case 'False'
                dampingFigure = figure;
            case 'True'
                dampingFigure = figure('visible','off');
        end
        % Create the damping figure
        hold on
        errorbar(frequency,Hpp,Hpp_Error,'o');
        plot(frequency,alphaFitData,'--')   
        hold off
        xlabel('Frequency (GHz)')
        ylabel('Peak to Peak Linewidth (Oe)')
        title(strcat(date,' Damping in ',sampleName));
        txt = ['Damping = ' num2str(damping) ' +/- ' num2str(dampingVariation)];
        xlimits = xlim;
        ylimits = ylim;
        axis([xlimits(1)*.999,xlimits(2)*1.001,ylimits(1)*.999,ylimits(2)*1.001])
        text(xlimits(1)*1.001,ylimits(2)*.999,txt)
        % Save the damping figure
        dampingFigureName = strcat('Damping',strcat(sampleName,'.fig'));
        set(dampingFigure,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
        savefig(dampingFigure,dampingFigureName);
        close all
        % Save all the individual fit figures
        for i = 1:numFiles
            switch suppressFigures
                case 'False'
                    figure;
                case 'True'
                    figure('visible','off')
            end
            plot(fitData.(strcat(namelist{i},'fit')),fitData.(strcat(namelist{i},'field')),fitData.(strcat(namelist{i},'data')))
            xlabel('FMR Field (Oe)')
            ylabel('FMR Power Derivative (arb)') 
            title(strcat(strcat('Lorentzian Fitting of FMR Power Derivative ',num2str(frequency(i))),' GHz'))
            fitFigureName = strcat('LorentzianFit',strcat(strcat(strcat(num2str(frequency(i)),'GHz'),'.fig')));
            set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
            savefig(fitFigureName)
            close all
        end

        % GAMMA CALCULATION
        % Calculate gamma and Ms from Lorentzian fits
        % Convert frequency to megahertz to match standard gamma units
        frequencyKittel = frequency.*1e3;
        for i = 1:numFiles
            tempFitStore = fitData.(strcat(namelist{i},'fit'));
            Hfmr(i) = tempFitStore.Hfmr;
        end
        switch geometry
            case 'IP'
                % Fit the data to the in plane Kittel equation
                [kittelFit,gamma,Ms] = IPKittel(frequencyKittel,Hfmr);
            case 'OOP'
                % Fit the data to the out of plane Kittel equation
                % UNTESTED, USE WITH CARE
                [kittelFit,gamma,Ms] = OOPKittel(frequencyKittel,Hfmr);      
        end
        % Plot the data and the in plane Kittel equation fit
        switch suppressFigures
                case 'False'
                    figure;
                case 'True'
                    figure('visible','off')
        end
        plot(kittelFit,frequencyKittel,Hfmr)
        % Find the 95% confidence interval for uncertainty purposes
        gammaConfidence = confint(kittelFit,0.95);
        gammaUpperbound = gammaConfidence(1,1);
        gammaLowerbound = gammaConfidence(2,1);
        gammaVariation = abs(gammaUpperbound - gammaLowerbound)/2;
        xlabel('Frequency (MHz)')
        ylabel('FMR Field (Oe)')
        title(strcat('Kittel Equation Fit ',geometry))
        xlimits = xlim;
        ylimits = ylim;
        axis([xlimits(1)*.999,xlimits(2)*1.001,ylimits(1)*.999,ylimits(2)*1.001])
        txt = ['Gamma = ' num2str(gamma) ' +/- ' num2str(gammaVariation)];
        text(xlimits(1)*1.001,ylimits(2)*.99,txt)
        % Name and save the figure
        kittelFigureName = strcat('Gamma',strcat(sampleName,'.fig'));
        set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
        savefig(kittelFigureName)
        close all
    case 'VNAFMR'
        
end
end

%% Function that calculates damping from Hpp values
function [alpha,alphaVariation,ycalc,deltaH0] = CalculateDamping(Hpp,frequency)
% Inputs
% Hpp = row array of peak to peak linewidths
% frequency = column array of fmr frequencies at which data was taken
% Calculate damping from given values of peak to peak linewidth and
% frequency
frequency = frequency';
% Make sure there are more than two data points so the standard error
% doesn't blow up
if length(frequency) <= 2
    warning('Damping standard error computation will not function properly. Consider taking more data points')
end
% Build the fit
frequency_fit = [ones(length(frequency),1), frequency'];
fit = frequency_fit\Hpp';
ycalc = frequency_fit*fit;
alpha = fit(2)*2.86/2/1e3;
% Calculate the standard error
stdError = sqrt(1/(length(frequency)-2)*sum((Hpp-ycalc').^2)/sum((frequency-mean(frequency)).^2));
alphaVariation = stdError*2.86/1e3;  
deltaH0 = fit(1);
end     

%% Function that finds the location of the Lorentzian/Gaussian via kmeans clustering or, failing that, user input
function [startField,endField] = FindLorentzian(field,data)
% Inputs
% field = field in Oe over which the data was taken
% data = FMR power derivative data

% SYSTEM INFORMATION
set(0,'units','pixels');
Pix_SS = get(0,'screensize');

% Determine where we think the Lorentzian is via Kmeans clustering
% Make sure that the number of points is divisible by the number of
% clusters
while rem(length(field),3) ~= 0
    field = field(1:end-1);
    data = data(1:end-1);
end
data = double(data);

% Build the features in dataArray
dataArray = zeros(length(field),2);
dataArray3 = abs(movmean(data,9))-mean(abs(movmean(data,9)));
dataArray3(dataArray3 < 0) = 0;
dataArray(:,1) = dataArray3;
dataArray(:,2) = movmean(dataArray(:,2),51);
dataArray(:,3) = abs(movingslope(data,7,1));
dataArray(:,4) = movstd(data,27);
% Perform kmeans clustering
[idx,C] = kmeans(dataArray,3,'Distance','sqeuclidean','MaxIter',10000,'OnlinePhase','On');
% Plot the data with points marking the beginning and end of the proposed
% lorentzian
figure('Renderer', 'painters', 'Position', [(.25)*Pix_SS(3) .1*Pix_SS(4) (.5)*Pix_SS(3) .5*Pix_SS(4)])
hold on
indexFirst = find(idx~=idx(1),1,'first');
indexLast = find(idx~=idx(1),1,'last');
plot(field,data,field(indexFirst),data(indexFirst),'o',field(indexLast),data(indexLast),'o')
xlabel('Field')
ylabel('YChannel Data')
title('Kmeans Clustering Attempt at Finding FMR Signal')
% Ask if we found the Lorentzian
answer = questdlg('Here is the Lorentzian computed by kmeans clustering. Is this correct?');
hold off
close all
% Deal with cases
if strcmp(answer,'Yes') == 1
    % If we found it well, use the points found by kmeans
    startField = field(indexFirst);
    endField = field(indexLast);
elseif strcmp(answer,'No') == 1
    % Ask the user where the Lorentzian is. 
    figure('Renderer', 'painters', 'Position', [(.25)*Pix_SS(3) .1*Pix_SS(4) (.5)*Pix_SS(3) .5*Pix_SS(4)])
    hold on
    plot(field,data)
    xlabel('FMR Field (Oe)')
    ylabel('yChannel Data')
    title('Select the bounds of the region you would like to fit. Left, then right')
    % Gather user input
    [x,y] = ginput(2);
    hold off
    close all
    startField = x(1);
    endField = x(2);
    % Error if the lower bound is above the upper bound
    if startField > endField
        error('Invalid bounds: Select the lowerbound of the fitting region first')
    end
else
    error('Must select Lorentzian bounds.')
end
end

%% Function that finds the location of the Lorentzian/Gaussian when the user specifies that the process is to be automated
function [startField,endField] = AutomaticFindLorentzian(field,data)
% Inputs
% field = field in Oe over which the data was taken
% data = FMR power derivative data

% SYSTEM INFORMATION
set(0,'units','pixels');
Pix_SS = get(0,'screensize');

% Determine where we think the Lorentzian is via Kmeans clustering
% Make sure that the number of points is divisible by the number of
% clusters
while rem(length(field),3) ~= 0
    field = field(1:end-1);
    data = data(1:end-1);
end
% Build the features in dataArray
dataArray = zeros(length(field),2);
dataArray3 = abs(movmean(data,9))-mean(abs(movmean(data,9)));
dataArray3(dataArray3 < 0) = 0;
dataArray(:,1) = dataArray3;
dataArray(:,2) = movmean(dataArray(:,2),51);
dataArray(:,3) = abs(movingslope(data,7,1));
dataArray(:,4) = movstd(data,27);
% Perform kmeans clustering
[idx,C] = kmeans(dataArray,3,'Distance','sqeuclidean','MaxIter',10000,'OnlinePhase','On');
indexFirst = find(idx~=idx(1),1,'first');
indexLast = find(idx~=idx(1),1,'last');
startField = field(indexFirst);
endField = field(indexLast);
end

%% Function which automatically fits a gaussian curve to provided data
function [Hpp,fit_output,fieldFit,dataFit] = FitGaussian(field,data,startField,endField,varargin)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the gaussian starts
% endField = field at which the gaussian ends
% varargin = lowerbound for Hpp

% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the gaussian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in gaussian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end

% STAGE 1 FIT
% Height, b, c all range. Hfmr and Hpp range within bounds.
heightInitial = abs(max(dataMaxMinRange));
bInitial = 0;
cInitial = 0;
HfmrInitial = fieldMaxMinRange(dataZeroIndex);
HppInitial = fieldFit(maxDataInRangeIndex)-fieldFit(minDataInRangeIndex);
bUpperbound = Inf;
bLowerbound = -Inf;
cUpperbound = Inf;
cLowerbound = -Inf;
heightUpperbound = Inf;
heightLowerbound = 0;
HppUpperbound = HppInitial+.1*HppInitial;
if isempty(varargin) == 1
    HppLowerbound = HppInitial-.1*HppInitial;
elseif (HppInitial-.1*HppInitial) > varargin{1}
    HppLowerbound = HppInitial-.1*HppInitial;
else
    HppLowerbound = varargin{1};
    HppUpperbound = varargin{1}+.1*varargin{1};
end
HfmrUpperbound = HfmrInitial+.01*HfmrInitial;
HfmrLowerbound = -HfmrInitial-.01*HfmrInitial;
% Build the fit objects
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[bInitial,cInitial,heightInitial,HppInitial,HfmrInitial],...
                        'Upper',[bUpperbound,cUpperbound,heightUpperbound,HppUpperbound,HfmrUpperbound],...
                        'Lower',[bLowerbound,cLowerbound,heightLowerbound,HppLowerbound,HfmrLowerbound]);
fitType = fittype('GaussianDerivativeFunction(b,c,height,Hpp,Hfmr,x)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'b','c','height','Hpp','Hfmr'});
stage1Fit = fit(fieldFit,dataFit,fitType,fitOptions);

% STAGE 2 FIT
% Hpp and Hfmr range to reduce error bounds
if isempty(varargin) == 1
    HppLowerbound_stage2 = 0;
else
    HppLowerbound_stage2 = varargin{1}-1*varargin{1};
end
fitOptions_stage2 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage1Fit.b,stage1Fit.c,stage1Fit.height,stage1Fit.Hpp,stage1Fit.Hfmr],...
                                'Upper',[stage1Fit.b,stage1Fit.c,stage1Fit.height,Inf,Inf],...
                                'Lower',[stage1Fit.b,stage1Fit.c,stage1Fit.height,HppLowerbound_stage2,-Inf],...
                                'Robust','LAR');
stage2Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage2);

% STAGE 3 FIT
if isempty(varargin) == 1
    HppLowerbound_stage3 = stage2Fit.Hpp - .1*stage2Fit.Hpp;
else
    HppLowerbound_stage3 = varargin{1}-.5*varargin{1};
end
fitOptions_stage3 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage2Fit.b,stage2Fit.c,stage2Fit.height,stage2Fit.Hpp,stage2Fit.Hfmr],...
                                'Upper',[stage2Fit.b,stage2Fit.c,stage2Fit.height,Inf,stage2Fit.Hfmr],...
                                'Lower',[stage2Fit.b,stage2Fit.c,stage2Fit.height,HppLowerbound_stage3,stage2Fit.Hfmr]);
stage3Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage3);
% Build the output variables
Hpp = stage3Fit.Hpp;
fit_output = stage3Fit;
Hfmr = stage3Fit.Hfmr;
end

%% Function that allows parameters in gaussian fitting to be manipulated by GUI
function [Hpp,fit_output,fieldFit,dataFit] = FitGaussianFreeInput(field,data,startField,endField,initial_val)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the lorentzian starts
% endField = field at which the lorentzian ends
% varargin = lowerbound for Hpp

% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the lorentzian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in lorentzian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end

% Stage 1 Fit
% Everything ranges
% Build the fit objects
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[initial_val.bInitial,initial_val.cInitial,initial_val.heightInitial,initial_val.HppInitial,initial_val.HfmrInitial],...
                        'Upper',[initial_val.bUpperbound,initial_val.cUpperbound,initial_val.heightUpperbound,initial_val.HppUpperbound,initial_val.HfmrUpperbound],...
                        'Lower',[initial_val.bLowerbound,initial_val.cLowerbound,initial_val.heightLowerbound,initial_val.HppLowerbound,initial_val.HfmrLowerbound]);
fitType = fittype('GaussianDerivativeFunction(b,c,height,Hpp,Hfmr,x)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'b','c','height','Hpp','Hfmr'});
stage1Fit = fit(fieldFit,dataFit,fitType,fitOptions);

% STAGE 2 FIT
% Hpp and Hfmr range to reduce error bounds
HppLowerbound_stage2 = stage1Fit.Hpp - initial_val.HppTolerance*stage1Fit.Hpp;
fitOptions_stage2 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage1Fit.b,stage1Fit.c,stage1Fit.height,stage1Fit.Hpp,stage1Fit.Hfmr],...
                                'Upper',[stage1Fit.b,stage1Fit.c,stage1Fit.height,Inf,Inf],...
                                'Lower',[stage1Fit.b,stage1Fit.c,stage1Fit.height,HppLowerbound_stage2,-Inf],...
                                'Robust','LAR');
stage2Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage2);

% STAGE 3 FIT
% Only Hpp ranges to further hone in on the correct value
HppLowerbound_stage3 = stage1Fit.Hpp - initial_val.HppTolerance*stage1Fit.Hpp;
fitOptions_stage3 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage2Fit.b,stage2Fit.c,stage2Fit.height,stage2Fit.Hpp,stage2Fit.Hfmr],...
                                'Upper',[stage2Fit.b,stage2Fit.c,stage2Fit.height,Inf,stage2Fit.Hfmr],...
                                'Lower',[stage2Fit.b,stage2Fit.c,stage2Fit.height,HppLowerbound_stage3,stage2Fit.Hfmr]);
stage3Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage3);
% Build the output variables
Hpp = stage3Fit.Hpp;
fit_output = stage3Fit;
Hfmr = stage3Fit.Hfmr;
end

%% Function that automatically fits lorentzians to data
function [Hpp,fit_output,fieldFit,dataFit] = FitLorentzian(field,data,startField,endField,varargin)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the lorentzian starts
% endField = field at which the lorentzian ends
% varargin = lowerbound for Hpp

% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the lorentzian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in lorentzian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end

% Stage 1 Fit
% a,b,c range. Theta,Hpp,Hfmr fixed
% Declare all initial values and bounds. Changing these in stage 1 can make
% a big difference in the quality of fit. If the fit you've found is poor,
% start here
aInitial = 0;
bInitial = 0;
cInitial = 0;
% Initialize Hpp by finding the width between the peaks
HppInitial = fieldFit(maxDataInRangeIndex)-fieldFit(minDataInRangeIndex);
% Initialize Hfmr by finding where the lorentzian passes through zero
HfmrInitial = fieldMaxMinRange(dataZeroIndex);
thetaInitial = 0;
aUpperbound = Inf;
aLowerbound = -Inf;
bUpperbound = Inf;
bLowerbound = -Inf;
cUpperbound = Inf;
cLowerbound = -Inf;
% Declare bounds for Hpp and Hfmr
HppUpperbound = HppInitial+.25*HppInitial;
if isempty(varargin) == 1
    HppLowerbound = HppInitial-.25*HppInitial;
elseif (HppInitial-.25*HppInitial) > varargin{1}
    HppLowerbound = HppInitial-.25*HppInitial;
else
    HppLowerbound = varargin{1};
    HppUpperbound = varargin{1}+.25*varargin{1};
end
HfmrUpperbound = HfmrInitial+.1*HfmrInitial;
HfmrLowerbound = -HfmrInitial-.1*HfmrInitial;
thetaUpperbound = 0;
thetaLowerbound = 0;
% Build the fit objects
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[aInitial,bInitial,cInitial,HppInitial,HfmrInitial,thetaInitial],...
                        'Upper',[aUpperbound,bUpperbound,cUpperbound,HppInitial,HfmrInitial,thetaUpperbound],...
                        'Lower',[aLowerbound,bLowerbound,cLowerbound,HppInitial,HfmrInitial,thetaLowerbound],...
                        'Robust','bisquare');
fitType = fittype('LorentzianDerivativeFunction(A,b,c,Hpp,Hfmr,theta,x)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'A','b','c','Hpp','Hfmr','theta'});               
stage1Fit = fit(fieldFit,dataFit,fitType,fitOptions);

% Stage 2 Fit
% Theta ranges. Else fixed
fitOptions_stage2 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,stage1Fit.theta],...
                                'Upper',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,180],...
                                'Lower',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,-180]);
stage2Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage2);

% Stage 3 Fit
% Hpp,Hfmr range. a,b,c,theta fixed
% HfmrUpperbound = max(fieldMaxMinRange);
% HfmrLowerbound = min(fieldMaxMinRange);
if isempty(varargin) == 1
    HppLowerbound_stage3 = 0;
else
    HppLowerbound_stage3 = varargin{1}-varargin{1}*.25;
end
fitOptions_stage3 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage2Fit.A,stage2Fit.b,stage2Fit.c,stage2Fit.Hpp,stage2Fit.Hfmr,stage2Fit.theta],...
                                'Upper',[stage2Fit.A,stage2Fit.b,stage2Fit.c,Inf,Inf,stage2Fit.theta],...
                                'Lower',[stage2Fit.A,stage2Fit.b,stage2Fit.c,HppLowerbound_stage3,-Inf,stage2Fit.theta]);
stage3Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage3);
% Build the output variables
Hpp = stage3Fit.Hpp;
fit_output = stage3Fit;
Hfmr = stage3Fit.Hfmr;
end

%% Function which fits a lorentzian curve to data based on user input
function [Hpp,fit_output,fieldFit,dataFit] = FitLorentzianFreeInput(field,data,startField,endField,initial_val)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the lorentzian starts
% endField = field at which the lorentzian ends
% varargin = lowerbound for Hpp

% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the lorentzian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in lorentzian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end

% Stage 1 Fit
% a,b,c range. Theta,Hpp,Hfmr fixed
thetaInitial = 0;
thetaUpperbound = 0;
thetaLowerbound = 0;
% Build the fit objects
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[initial_val.aInitial,initial_val.bInitial,initial_val.cInitial,initial_val.HppInitial,initial_val.HfmrInitial,thetaInitial],...
                        'Upper',[initial_val.aUpperbound,initial_val.bUpperbound,initial_val.cUpperbound,initial_val.HppUpperbound,initial_val.HfmrUpperbound,thetaUpperbound],...
                        'Lower',[initial_val.aLowerbound,initial_val.bLowerbound,initial_val.cLowerbound,initial_val.HppLowerbound,initial_val.HfmrLowerbound,thetaLowerbound]);
fitType = fittype('LorentzianDerivativeFunction(A,b,c,Hpp,Hfmr,theta,x)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'A','b','c','Hpp','Hfmr','theta'});
stage1Fit = fit(fieldFit,dataFit,fitType,fitOptions);

% Stage 2 Fit
% Theta ranges. Else fixed
fitOptions_stage2 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,stage1Fit.theta],...
                                'Upper',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,180],...
                                'Lower',[stage1Fit.A,stage1Fit.b,stage1Fit.c,stage1Fit.Hpp,stage1Fit.Hfmr,-180]);
stage2Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage2);

% Stage 3 Fit
% Hpp,Hfmr range. a,b,c,theta fixed
% HfmrUpperbound = max(fieldMaxMinRange);
% HfmrLowerbound = min(fieldMaxMinRange);
fitOptions_stage3 = fitoptions('Method','NonlinearLeastSquares',...
                                'Startpoint',[stage2Fit.A,stage2Fit.b,stage2Fit.c,stage2Fit.Hpp,stage2Fit.Hfmr,stage2Fit.theta],...
                                'Upper',[stage2Fit.A,stage2Fit.b,stage2Fit.c,stage2Fit.Hpp+initial_val.HppTolerance*stage2Fit.Hpp,stage2Fit.Hfmr+initial_val.HfmrTolerance*stage2Fit.Hfmr,stage2Fit.theta],...
                                'Lower',[stage2Fit.A,stage2Fit.b,stage2Fit.c,stage2Fit.Hpp-initial_val.HppTolerance*stage2Fit.Hpp,stage2Fit.Hfmr-initial_val.HfmrTolerance*stage2Fit.Hfmr,stage2Fit.theta]);
stage3Fit = fit(fieldFit,dataFit,fitType,fitOptions_stage3);
% Build the output variables
Hpp = stage3Fit.Hpp;
fit_output = stage3Fit;
Hfmr = stage3Fit.Hfmr;
end

%% Equation used for Gaussian fitting
function [y] = GaussianDerivativeFunction(b,c,height,Hpp,Hfmr,x)
x = x';
y = zeros(length(x));
y = b + c.*x + -1*((Hpp/2)*exp(.5))*(height*(-x+Hfmr).*exp(-((-x+Hfmr).^2)/(2*(Hpp/2)^2)))/((Hpp/2)^2);
y = y';
end

%% Routine defining the interactive gaussian fitting process
function [Hpp,fit_output,fieldFit,dataFit] = InteractiveFitGaussian(field,data,startField,endField)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the lorentzian starts
% endField = field at which the lorentzian ends
% varargin = lowerbound for Hpp

% SYSTEM INFORMATION
set(0,'units','pixels');
Pix_SS = get(0,'screensize');

% Function
% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the lorentzian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in lorentzian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end
% Allow user to alter values for initialization bounds
isSatisfied = 'No';
% Provide universal fit parameters as a starting point
initial_val.HppTolerance = .25;
initial_val.HfmrTolerance = .05;
initial_val.bInitial = 0;
initial_val.cInitial = 0;
initial_val.heightInitial = 0;
initial_val.HppInitial = fieldFit(maxDataInRangeIndex)-fieldFit(minDataInRangeIndex);
initial_val.HfmrInitial = fieldMaxMinRange(dataZeroIndex);
initial_val.bUpperbound = Inf;
initial_val.cUpperbound = Inf;
initial_val.heightUpperbound = Inf;
initial_val.HppUpperbound = initial_val.HppInitial+initial_val.HppTolerance*initial_val.HppInitial;
initial_val.HfmrUpperbound = initial_val.HfmrInitial+initial_val.HfmrTolerance*initial_val.HfmrInitial;
initial_val.bLowerbound = -Inf;
initial_val.cLowerbound = -Inf;
initial_val.heightLowerbound = -Inf;
initial_val.HppLowerbound = initial_val.HppInitial-initial_val.HppTolerance*initial_val.HppInitial;
initial_val.HfmrLowerbound = initial_val.HfmrInitial-initial_val.HfmrTolerance*initial_val.HfmrInitial;
userIteration = 1;
% Let the user change the initial values and bounds until they are
% satisfied with the fit
while strcmp(isSatisfied,'No') == 1
    % Fit based on previously defined attributes
    [Hpp,fit_output,fieldFit,dataFit] = FitGaussianFreeInput(field,data,startField,endField,initial_val);
    figure('Renderer', 'painters', 'Position', [(.61)*Pix_SS(3) .5*Pix_SS(4) (.4)*Pix_SS(3) .4*Pix_SS(4)])
    plot(fit_output,fieldFit,dataFit)
    xlabel('Field (Oe)')
    ylabel('FMR Power Derivative (arb)')
    title(strcat("Test Gaussian Fit Number ",num2str(userIteration)))
    if userIteration == 1
    else
         isSatisfiedAnswer = questdlg('Is the Gaussian fit acceptable?');
         switch isSatisfiedAnswer
             case 'Yes'
                 close all
                 break
             case 'No'
             otherwise
                 close all
                 error('Must provide information about fit quality')
         end
    end
    % Prompt the user for new values
    prompt = {strcat("Enter Initial Value for b: previous value was ",num2str(initial_val.bInitial)),...
              strcat("Enter Lowerbound for b: previous value was ",num2str(initial_val.bLowerbound)),...
              strcat("Enter Upperbound for b: previous value was ",num2str(initial_val.bUpperbound)),...
              strcat("Enter Initial Value for c: previous value was ",num2str(initial_val.cInitial)),...
              strcat("Enter Lowerbound for c: previous value was ",num2str(initial_val.cLowerbound)),...
              strcat("Enter Upperbound for c: previous value was ",num2str(initial_val.cUpperbound)),...
              strcat("Enter Initial Valvue for height: previous value was ",num2str(initial_val.heightInitial)),...
              strcat("Enter Lowerbound for height: previous value was ",num2str(initial_val.heightLowerbound)),...
              strcat("Enter Upperbound for height: previous value was ",num2str(initial_val.heightUpperbound)),...
              strcat("Enter Initial Value for Hpp: previous value was ",num2str(initial_val.HppInitial)),...
              strcat("Enter Lowerbound for Hpp: previous value was ",num2str(initial_val.HppLowerbound)),...
              strcat("Enter Upperbound for Hpp: previous value was ",num2str(initial_val.HppUpperbound)),...
              strcat("Enter Initial Value for Hfmr: previous value was ",num2str(initial_val.HfmrInitial)),...
              strcat("Enter Lowerbound for Hfmr: previous value was ",num2str(initial_val.HfmrLowerbound)),...
              strcat("Enter Upperbound for Hfmr: previous value was ",num2str(initial_val.HfmrUpperbound)),...
              strcat("Enter Tolerance for Hpp: previous value was ",num2str(initial_val.HppTolerance)),...
              strcat("Enter Tolerance for Hfmr: previous value was ",num2str(initial_val.HfmrTolerance))};
    dlgtitle = "Fit Parameters";
    dims = [1 75];
    definput = {num2str(initial_val.bInitial),num2str(initial_val.bLowerbound),num2str(initial_val.bUpperbound),num2str(initial_val.cInitial),num2str(initial_val.cLowerbound),...
                num2str(initial_val.cUpperbound),num2str(initial_val.heightInitial),num2str(initial_val.heightLowerbound),num2str(initial_val.heightUpperbound),num2str(initial_val.HppInitial),...
                num2str(initial_val.HppLowerbound),num2str(initial_val.HppUpperbound),num2str(initial_val.HfmrInitial),num2str(initial_val.HfmrLowerbound),num2str(initial_val.HfmrUpperbound),...
                num2str(initial_val.HppTolerance),num2str(initial_val.HfmrTolerance)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    close all
    % Set all the initial values to those dictated by the user
    initial_val.HppTolerance = str2num(answer{16});
    initial_val.HfmrTolerance = str2num(answer{17});
    initial_val.bInitial = str2num(answer{1});
    initial_val.cInitial = str2num(answer{4});
    initial_val.heightInitial = str2num(answer{7});
    initial_val.HppInitial = str2num(answer{10});
    initial_val.HfmrInitial = str2num(answer{13});
    initial_val.bUpperbound = str2num(answer{3});
    initial_val.cUpperbound = str2num(answer{6});
    initial_val.heightUpperbound = str2num(answer{9});
    initial_val.HppUpperbound = str2num(answer{12});
    initial_val.HfmrUpperbound = str2num(answer{15});
    initial_val.bLowerbound = str2num(answer{2});
    initial_val.cLowerbound = str2num(answer{5});
    initial_val.heightLowerbound = str2num(answer{8});
    initial_val.HppLowerbound = str2num(answer{11});
    initial_val.HfmrLowerbound = str2num(answer{14});
    userIteration = userIteration + 1;
end
end

%% Routine defining the interactive Lorentzian fitting process
% To do: combine this with InteractiveFitGaussian to reduce
% repitition
function [Hpp,fit_output,fieldFit,dataFit] = InteractiveFitLorentzian(field,data,startField,endField)
% Inputs
% field = FMR field sweep
% data = FMR power derivative
% startField = field at which the lorentzian starts
% endField = field at which the lorentzian ends
% varargin = lowerbound for Hpp

% SYSTEM INFORMATION
set(0,'units','pixels');
Pix_SS = get(0,'screensize');

% Function
% Find the point at which data really starts since sometimes it bounces
% around at zero field for a while
realDataStartIndex = find(field-field(1),1);
field = field(realDataStartIndex:end);
data = data(realDataStartIndex:end);
% Find the index of the start and end of the lorentzian
startIndex = interp1(unique(field),1:length(unique(field)),startField,'nearest');
endIndex = interp1(unique(field),1:length(unique(field)),endField,'nearest');
fieldFit = field(startIndex:endIndex);
dataFit = data(startIndex:endIndex);
maxDataInRange = max(dataFit);
minDataInRange = min(dataFit);
% Find location of two peaks in lorentzian derivative
maxDataInRangeIndex = find(dataFit == maxDataInRange,1,'First');
minDataInRangeIndex = find(dataFit == minDataInRange,1,'Last');
dataMaxMinRange = dataFit(minDataInRangeIndex:maxDataInRangeIndex);
fieldMaxMinRange = fieldFit(minDataInRangeIndex:maxDataInRangeIndex);
dataZeroIndex = find(dataMaxMinRange == min(abs(dataMaxMinRange)),1,'First');
% Find the index of the center of the distribution, where it passes closest
% to zero
if isempty(dataZeroIndex) == 1
    dataZeroIndex = find(dataMaxMinRange == -1*min(abs(dataMaxMinRange)),1,'First');
end
% STAGE 1 USER INPUT
% Allow user to alter values for initialization bounds
isSatisfied = 'No';
% Provide universal fit parameters as a starting point
initial_val.HppTolerance = .25;
initial_val.HfmrTolerance = .05;
initial_val.aInitial = 0;
initial_val.bInitial = 0;
initial_val.cInitial = 0;
initial_val.HppInitial = fieldFit(maxDataInRangeIndex)-fieldFit(minDataInRangeIndex);
initial_val.HfmrInitial = fieldMaxMinRange(dataZeroIndex);
initial_val.aUpperbound = Inf;
initial_val.bUpperbound = Inf;
initial_val.cUpperbound = Inf;
initial_val.HppUpperbound = initial_val.HppInitial+initial_val.HppTolerance*initial_val.HppInitial;
initial_val.HfmrUpperbound = initial_val.HfmrInitial+initial_val.HfmrTolerance*initial_val.HfmrInitial;
initial_val.aLowerbound = -Inf;
initial_val.bLowerbound = -Inf;
initial_val.cLowerbound = -Inf;
initial_val.HppLowerbound = initial_val.HppInitial-initial_val.HppTolerance*initial_val.HppInitial;
initial_val.HfmrLowerbound = initial_val.HfmrInitial-initial_val.HfmrTolerance*initial_val.HfmrInitial;
userIteration = 1;
while strcmp(isSatisfied,'No') == 1
    [Hpp,fit_output,fieldFit,dataFit] = FitLorentzianFreeInput(field,data,startField,endField,initial_val);
    figure('Renderer', 'painters', 'Position', [(.61)*Pix_SS(3) .5*Pix_SS(4) (.4)*Pix_SS(3) .4*Pix_SS(4)])
    plot(fit_output,fieldFit,dataFit)
    xlabel('Field (Oe)')
    ylabel('FMR Power Derivative (arb)')
    title(strcat("Test Lorentzian Fit Number ",num2str(userIteration)))
    if userIteration == 1
    else
         isSatisfiedAnswer = questdlg('Is the Lorentzian fit acceptable?');
         switch isSatisfiedAnswer
             case 'Yes'
                 close all
                 break
             case 'No'
             otherwise
                 close all
                 error('Must provide information about fit quality')
         end
    end
    prompt = {strcat("Enter Initial Value for a: previous value was ",num2str(initial_val.aInitial)),...
              strcat("Enter Lowerbound for a: previous value was ",num2str(initial_val.aLowerbound)),...
              strcat("Enter Upperbound for a: previous value was ",num2str(initial_val.aUpperbound)),...
              strcat("Enter Initial Value for b: previous value was ",num2str(initial_val.bInitial)),...
              strcat("Enter Lowerbound for b: previous value was ",num2str(initial_val.bLowerbound)),...
              strcat("Enter Upperbound for b: previous value was ",num2str(initial_val.bUpperbound)),...
              strcat("Enter Initial Valvue for c: previous value was ",num2str(initial_val.cInitial)),...
              strcat("Enter Lowerbound for c: previous value was ",num2str(initial_val.cLowerbound)),...
              strcat("Enter Upperbound for c: previous value was ",num2str(initial_val.cUpperbound)),...
              strcat("Enter Initial Value for Hpp: previous value was ",num2str(initial_val.HppInitial)),...
              strcat("Enter Lowerbound for Hpp: previous value was ",num2str(initial_val.HppLowerbound)),...
              strcat("Enter Upperbound for Hpp: previous value was ",num2str(initial_val.HppUpperbound)),...
              strcat("Enter Initial Value for Hfmr: previous value was ",num2str(initial_val.HfmrInitial)),...
              strcat("Enter Lowerbound for Hfmr: previous value was ",num2str(initial_val.HfmrLowerbound)),...
              strcat("Enter Upperbound for Hfmr: previous value was ",num2str(initial_val.HfmrUpperbound)),...
              strcat("Enter Tolerance for Hpp: previous value was ",num2str(initial_val.HppTolerance)),...
              strcat("Enter Tolerance for Hfmr: previous value was ",num2str(initial_val.HfmrTolerance))};
    dlgtitle = "Fit Parameters";
    dims = [1 75];
    definput = {num2str(initial_val.aInitial),num2str(initial_val.aLowerbound),num2str(initial_val.aUpperbound),num2str(initial_val.bInitial),num2str(initial_val.bLowerbound),...
                num2str(initial_val.bUpperbound),num2str(initial_val.cInitial),num2str(initial_val.cLowerbound),num2str(initial_val.cUpperbound),num2str(initial_val.HppInitial),...
                num2str(initial_val.HppLowerbound),num2str(initial_val.HppUpperbound),num2str(initial_val.HfmrInitial),num2str(initial_val.HfmrLowerbound),num2str(initial_val.HfmrUpperbound),...
                num2str(initial_val.HppTolerance),num2str(initial_val.HfmrTolerance)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    close all
    initial_val.HppTolerance = str2num(answer{16});
    initial_val.HfmrTolerance = str2num(answer{17});
    initial_val.aInitial = str2num(answer{1});
    initial_val.bInitial = str2num(answer{4});
    initial_val.cInitial = str2num(answer{7});
    initial_val.HppInitial = str2num(answer{10});
    initial_val.HfmrInitial = str2num(answer{13});
    initial_val.aUpperbound = str2num(answer{3});
    initial_val.bUpperbound = str2num(answer{6});
    initial_val.cUpperbound = str2num(answer{9});
    initial_val.HppUpperbound = str2num(answer{12});
    initial_val.HfmrUpperbound = str2num(answer{15});
    initial_val.aLowerbound = str2num(answer{2});
    initial_val.bLowerbound = str2num(answer{5});
    initial_val.cLowerbound = str2num(answer{8});
    initial_val.HppLowerbound = str2num(answer{11});
    initial_val.HfmrLowerbound = str2num(answer{14});
    userIteration = userIteration + 1;
end
end

%% Function for Fitting the in plane kittel equation
function [kittelFit,gamma,Ms] = IPKittel(frequency,Hfmr)
% Inputs
% frequency = column array of frequencies data was taken at
% Hfmr = row array of FMR fields at the center of the Lorentzian
% Fit the in plane Kittel equation to provided data
Hfmr = Hfmr';
% Build the fit
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[2.86,1750]);
fitType = fittype('IPKittelEquation(x,gamma,Ms)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'gamma','Ms'});
stage1Fit = fit(frequency,Hfmr,fitType,fitOptions);
% Build output variables
gamma = stage1Fit.gamma;
Ms = stage1Fit.Ms;
kittelFit = stage1Fit;
end

%% The in plane Kittel equation
function [y] = IPKittelEquation(x,gamma,Ms)
% The in plane Kittel equation. Y is FMR field, x is FMR frequency
y = (-Ms + sqrt((Ms).^2 + 4.*x.^2./gamma^2))./2;
end

%% Function that handles data import
function [field,ychannelData] = LoadFMRData(filename)
% Loads FMR data into a structure array
data = tdfread(filename);
ychannelData = data.Ychannel;
field = data.Field0x28G0x29;
end

%% Equation defining the Lorentzian fit
function [y] = LorentzianDerivativeFunction(A,b,c,Hpp,Hfmr,theta,x)
% Function to which we fit the FMR power derivative. 
x = x';
y = zeros(length(x));
y = b + c.*x + (-16).*A.*cos(theta.*pi./180).*Hpp.*sqrt(3).*(x-Hfmr)./pi./(4.*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2).^2 ...
    +   4.*A.*sin(theta.*pi./180).*((-4).*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2)/pi/(4.*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2).^2;
y = y';
end

%% Function that calculated a moving slope. Used for kmeans feature generation
function Dvec = movingslope(vec,supportlength,modelorder,dt)
% Robust calculation of moving slope for a given order and number of points
% around a center value
if (nargin==0)
  help movingslope
  return
end
if ~isvector(vec)
  error('vec must be a row or column vector')
end
n = length(vec);
% supply defaults
if (nargin<4) || isempty(dt)
  dt = 1;
end
if (nargin<3) || isempty(modelorder)
  modelorder = 1;
end
if (nargin<2) || isempty(supportlength)
  supportlength = 3;
end
% check the parameters for problems
if (length(supportlength)~=1) || (supportlength<=1) || (supportlength>n) || (supportlength~=floor(supportlength))
  error('supportlength must be a scalar integer, >= 2, and no more than length(vec)')
end
if (length(modelorder)~=1) || (modelorder<1) || (modelorder>min(10,supportlength-1)) || (modelorder~=floor(modelorder))
  error('modelorder must be a scalar integer, >= 1, and no more than min(10,supportlength-1)')
end
if (length(dt)~=1) || (dt<0)
  error('dt must be a positive scalar numeric variable')
end
% now build the filter coefficients to estimate the slope
if mod(supportlength,2) == 1
  parity = 1; % odd parity
else
  parity = 0;
end
s = (supportlength-parity)/2;
t = ((-s+1-parity):s)';
coef = getcoef(t,supportlength,modelorder);
% Apply the filter to the entire vector
f = filter(-coef,1,vec);
Dvec = zeros(size(vec));
Dvec(s+(1:(n-supportlength+1))) = f(supportlength:end);
% patch each end
vec = vec(:);
for i = 1:s
  % patch the first few points
  t = (1:supportlength)' - i;
  coef = getcoef(t,supportlength,modelorder);
  
  Dvec(i) = coef*vec(1:supportlength);
  
  % patch the end points
  if i<(s + parity)
    t = (1:supportlength)' - supportlength + i - 1;
    coef = getcoef(t,supportlength,modelorder);
    Dvec(n - i + 1) = coef*vec(n + (0:(supportlength-1)) + 1 - supportlength);
  end
end
% scale by the supplied spacing
Dvec = Dvec/dt;
end % mainline end
% =========================================================
% Functin that computes the filter coefficients for moving slope
% calculation
function coef = getcoef(t,supportlength,modelorder)
% Note: bsxfun would have worked here as well, but some people
% might not yet have that release of matlab.
A = repmat(t,1,modelorder+1).^repmat(0:modelorder,supportlength,1);
pinvA = pinv(A);
% we only need the linear term
coef = pinvA(2,:);
end % nested function end

%% Function used for out of plane Kittel equation fitting
function [kittelFit,gamma,Ms] = OOPKittel(frequency,Hfmr)
% Inputs
% frequency = column array of frequencies data was taken at
% Hfmr = row array of FMR fields at the center of the Lorentzian
% Fit the in plane Kittel equation to provided data
Hfmr = Hfmr';
% Build the fit
fitOptions = fitoptions('Method','NonlinearLeastSquares',...
                        'Startpoint',[2.86,1750]);
fitType = fittype('OOPKittelEquation(x,gamma,Ms)',...
                    'independent',{'x'},...
                    'dependent',{'y'},...
                    'coefficients',{'gamma','Ms'});
stage1Fit = fit(frequency,Hfmr,fitType,fitOptions);
% Build output variables
gamma = stage1Fit.gamma;
Ms = stage1Fit.Ms;
kittelFit = stage1Fit;
end

%% The out of plane Kittel equation
function [y] = OOPKittelEquation(x,gamma,Ms)
% The out of plane Kittel equation. X is operational frequency, Y is FMR
% field
y = x./gamma + 4*pi*Ms;
end


%% Incomplete VNA Fitting Pipeline - Complete when I have some VNAFMR data to fit
function [field,real,imaginary] = LoadVNAData(filename)
data = csvread(filename,1,0);
field = data(:,2);
real = data(:,3);
imaginary = data(:,4);
end

function [startField,endField] = FindVNASignal(field,realData,imaginaryData)
set(0,'units','pixels');
Pix_SS = get(0,'screensize');
while rem(length(field),3) ~= 0
    field = field(1:end-1);
    realData = realData(1:end-1);
    imaginaryData = imaginaryData(1:end-1);
end
% Build features
% First off of real data
dataArray3 = abs(movmean(realData,9))-mean(abs(movmean(realData,9)));
dataArray3(dataArray3 < 0) = 0;
dataArray(:,1) = dataArray3;
dataArray(:,2) = movmean(realData,21);
for i = 1:length(realData)
    if dataArray(i,2) < -1*max(dataArray(:,2))*.1 + max(dataArray(:,2))
        dataArray(i,2) = 0;
    elseif dataArray(i,2) > max(dataArray(:,2))*.1 + max(dataArray(:,2))
    end
end
dataArray(:,3) = abs(movingslope(realData,7,1));
dataArray(:,4) = movstd(realData,27);
% Then off of imaginary data
dataArray4 = abs(movmean(imaginaryData,9))-mean(abs(movmean(imaginaryData,9)));
dataArray4(dataArray4 < 0) = 0;
dataArray(:,5) = dataArray4;
dataArray(:,6) = movmean(imaginaryData,21);
for i = 1:length(imaginaryData)
    if dataArray(i,6) < -1*max(dataArray(:,6))*.1 + max(dataArray(:,6))
        dataArray(i,6) = 0;
    elseif dataArray(i,6) > max(dataArray(:,6))*.1 + max(dataArray(:,6))
    end
end
dataArray(:,7) = abs(movingslope(imaginaryData,7,1));
dataArray(:,8) = movstd(imaginaryData,27);
% Perform the clustering
[idx,C] = kmeans(dataArray,3,'Distance','sqeuclidean','MaxIter',10000,'OnlinePhase','On');
% Plot the data with points marking the beginning and end of the proposed
% lorentzian
figure('Renderer', 'painters', 'Position', [(.25)*Pix_SS(3) .1*Pix_SS(4) (.5)*Pix_SS(3) .5*Pix_SS(4)])
hold on
indexFirst = find(idx~=idx(1),1,'first');
indexLast = find(idx~=idx(1),1,'last');
indexLength = indexLast - indexFirst;
% Pad the range a small asymmetric amount
indexFirst = round(indexFirst - indexLength*.35);
indexLast = round(indexLast + indexLast*.15);
yyaxis left
plot(field,realData)
ylabel('Real Part of S21')
yyaxis right
plot(field,imaginaryData)
xlabel('Field')
ylabel('Imaginary Part of S21')
title('Kmeans Clustering Attempt at Finding FMR Signal')
% Ask if we found the signal
answer = questdlg('Here is the Lorentzian computed by kmeans clustering. Is this correct?');
hold off
close all
% Deal with cases
if strcmp(answer,'Yes') == 1
    % If we found it well, use the points found by kmeans
    startField = field(indexFirst);
    endField = field(indexLast);
elseif strcmp(answer,'No') == 1
    % Ask the user where the signal is. 
    figure('Renderer', 'painters', 'Position', [(.25)*Pix_SS(3) .1*Pix_SS(4) (.5)*Pix_SS(3) .5*Pix_SS(4)])
    hold on
    yyaxis left
    plot(field,realData,field(indexFirst),realData(indexFirst),'o',field(indexLast),realData(indexLast),'o')
    ylabel('Real Part of S21')
    yyaxis right
    plot(field,imaginaryData,field(indexFirst),imaginaryData(indexFirst),'o',field(indexLast),imaginaryData(indexLast),'o')
    xlabel('Field')
    ylabel('Imaginary Part of S21')
    title('Select the bounds of the region you would like to fit. Left, then right')
    % Gather user input
    [x,y] = ginput(2);
    hold off
    close all
    startField = x(1);
    endField = x(2);
    % Error if the lower bound is above the upper bound
    if startField > endField
        error('Invalid bounds: Select the lowerbound of the fitting region first')
    end
else
    error('Must select Lorentzian bounds.')
end
end

% Define the equation to which VNAFMR data is to be fit
% Need to go back and re-derive, this doesn't look right at all.
function [y] = RealDataFit(Si,Di,Z,Xr,Xi,Hk,Ms,F,G,Dh,x)
y = Si + Di.*x + Z*( (Dh.*F/G).*(4.*pi.*Ms.*Xr.*(x + Hk - 4.*pi.*Ms)+4.*pi.*Ms.*Xi.*Dh./2)...
    +(4.*pi.*Ms.*Xr.*Dh./2-4.*pi.*Ms.*Xi.*(x + Hk - 4.*pi.*Ms)).*(x.^2+2.*x.*Hk-2.*x.*4.*pi.*Ms+Hk.^2-2.*Hk.*4.*pi.*Ms+(4.*pi.*Ms).^2-(F./G).^2) )./((Xr.^2+Xi.^2).*(((x + Hk - 4.*pi.*Ms).^2-(F./G).^2).^2+(Dh.*F./G).^2));
end















     

