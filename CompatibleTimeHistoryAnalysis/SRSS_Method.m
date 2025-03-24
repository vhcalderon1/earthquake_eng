% =========================================================================
% EARTHQUAKE SPECTRAL ANALYSIS TOOL: SRSS SCALING AND CODE COMPLIANCE CHECK
% =========================================================================
% Author: Victor Calderón (August 2020)
% Updated: Jefferson De la Cuba (February 2025)
% --------------------------------------------------------------------------
% Processes earthquake ground motion spectra to compute SRSS-
% scaled acceleration/displacement spectra and compare against 
% building code reference spectrum (E.031). Produces comparative
% visualizations of averaged earthquake spectra vs. code requirements.
%
% INPUTS:
% - Reference spectrum:    E031_IsolatedBaseSpectrum.txt
% - Earthquake datasets:   Paired *_EW_Spectrum.txt and *_NS_Spectrum.txt files for 7 events       
% - Path configuration:    Datasets located in ../datasets/
% 
% OUTPUTS:
% - SRSS_Acceleration_Comparison.png
% - SRSS_Displacement_Comparison.png 
% - Creates/uses ../outputs/ directory for results
%
% =========================================================================

clear; close all; clc;

%% Directory and File Configuration
projectRoot = fileparts(mfilename('fullpath'));
datasetsPath = fullfile(projectRoot, '..', 'datasets');
outputsPath  = fullfile(projectRoot, '..', 'outputs');

% Create outputs directory if non-existent
if ~exist(outputsPath, 'dir')
    mkdir(outputsPath);
end

% Building Code Reference Spectrum (TXT: period [s], Sa [g])
refSpectrumData = load(fullfile(datasetsPath, 'E031_IsolatedBaseSpectrum.txt'));

% Earthquake dataset list (Paired EW/NS components)
earthquakeList = {...
    'Arequipa2001';...
    'Lima1966';...
    'Lima1974';...
    'Pisco2007';...
    'Concepcion2010';...
    'Curico2010';...
    'Hualane2010'};

% SRSS Component Scaling Factors (EW = 90%, NS = 44%)
factor_EW = 0.90;
factor_NS = 0.44;
%% Load Data, Compute SRSS-Scaled Acceleration Spectra, and Average Them
allScaledSpectra = [];  % Each column corresponds to an event
globalPeriods = [];      % Common period grid for all events

for i = 1:numel(earthquakeList)
    eqName = earthquakeList{i};
    ewFile = [eqName, '_EW_Spectrum.txt'];
    nsFile = [eqName, '_NS_Spectrum.txt'];
    
    % Load paired spectrum data
    [dataEW, dataNS] = loadPair(ewFile, nsFile, datasetsPath);
    
    % Extract period vectors and acceleration values
    periodsEW = dataEW(:,1);
    accelEW   = dataEW(:,2);
    
    periodsNS = dataNS(:,1);
    accelNS   = dataNS(:,2);
    
    % Check if period vectors are the same. If not, interpolate NS onto EW periods.
    if (length(periodsEW) ~= length(periodsNS)) || any(abs(periodsEW - periodsNS) > 1e-6)
        accelNS = interp1(periodsNS, accelNS, periodsEW, 'linear', 'extrap');
        currentPeriods = periodsEW;
    else
        currentPeriods = periodsEW; % (or periodsNS—they are identical)
    end
    
    % Compute SRSS-scaled acceleration spectrum using the vectorized function.
    srssAccel = computeSRSS(accelEW, accelNS, factor_EW, factor_NS);
    
    % For the first event, set the global period grid.
    if isempty(globalPeriods)
        globalPeriods = currentPeriods;
    else
        % If current event's period grid differs from the global one, interpolate.
        if (length(currentPeriods) ~= length(globalPeriods)) || any(abs(currentPeriods - globalPeriods) > 1e-6)
            srssAccel = interp1(currentPeriods, srssAccel, globalPeriods, 'linear', 'extrap');
        end
    end
    
    % Collect the computed SRSS acceleration spectrum
    allScaledSpectra = [allScaledSpectra, srssAccel];
end

% Calculate the average SRSS acceleration spectrum over all events.
averageSRSSAccel = mean(allScaledSpectra, 2);

%% Compute Pseudo-Displacement Spectra from the Acceleration Data
gravity = 9.81;  % Gravity in m/s²
averageSRSSDisp = computeSRSS_Displacement(globalPeriods, averageSRSSAccel, gravity);
refDisp         = computeSRSS_Displacement(refSpectrumData(:,1), refSpectrumData(:,2), gravity);

%% Plot and Save the Acceleration Spectra Comparison

% Create a figure handle for acceleration spectra
figAccel = figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
hold on;
% Plot all individual SRSS-scaled acceleration spectra (light gray lines)
for i = 1:size(allScaledSpectra,2)
    plotSpectrum(globalPeriods, allScaledSpectra(:,i), 'Individual');
end
% Plot the building code reference acceleration spectrum (solid black line)
hRef = plotSpectrum(refSpectrumData(:,1), refSpectrumData(:,2), 'Reference');
% Plot the average SRSS acceleration spectrum (dashed black line)
hAvg = plotSpectrum(globalPeriods, averageSRSSAccel, 'Average');
xlabel('Period (s)', 'FontSize', 14);
ylabel('Spectral Acceleration (Sa/g)', 'FontSize', 14);
legend([hRef, hAvg], {'Building Code E.031', 'Average SRSS Spectrum'},...
    'Location', 'northeast', 'FontSize', 12);
title('SRSS-Scaled Acceleration Spectra');
grid on;
drawnow;  % Ensure the figure is fully rendered

% Export the acceleration spectra figure using the explicit handle
exportgraphics(figAccel, fullfile(outputsPath, 'SRSS_Acceleration_Comparison.png'), 'Resolution', 300);
close(figAccel);  % Now close the figure


%% Plot and Save the Pseudo-Displacement Spectra Comparison

% Create a figure handle for displacement spectra
figDisp = figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
hold on;
% Plot the building code reference displacement spectrum
hRefDisp = plotSpectrum(refSpectrumData(:,1), refDisp, 'Reference');
% Plot the average pseudo-displacement spectrum from SRSS acceleration
hAvgDisp = plotSpectrum(globalPeriods, averageSRSSDisp, 'Average');
xlabel('Period (s)', 'FontSize', 14);
ylabel('Pseudo-Displacement (m)', 'FontSize', 14); % Corrected unit label
legend([hRefDisp, hAvgDisp], {'Reference Disp.', 'Average SRSS Disp.'},...
    'Location', 'northeast', 'FontSize', 12);
title('Pseudo-Displacement Spectra');
grid on;
drawnow;  % Ensure the figure is fully rendered

% Export the displacement spectra figure using the explicit handle
exportgraphics(figDisp, fullfile(outputsPath, 'SRSS_Displacement_Comparison.png'), 'Resolution', 300);
close(figDisp);  % Close the figure when done


%% ------------------- Local Functions -------------------

function [ewData, nsData] = loadPair(ewFile, nsFile, datasetsPath)
% loadPair loads a pair of spectrum files (East-West and North-South) 
% from the datasets folder.
%
% Inputs:
%   ewFile       - Filename for the East-West spectrum data.
%   nsFile       - Filename for the North-South spectrum data.
%   datasetsPath - Path to the datasets folder.
%
% Outputs:
%   ewData       - Data loaded from the East-West file.
%   nsData       - Data loaded from the North-South file.

    ewData = load(fullfile(datasetsPath, ewFile));
    nsData = load(fullfile(datasetsPath, nsFile));
end

function srss = computeSRSS(spectrumEW, spectrumNS, factorEW, factorNS)
% computeSRSS computes the SRSS-scaled acceleration spectrum.
%
% Inputs:
%   spectrumEW - Acceleration data for the East-West direction.
%   spectrumNS - Acceleration data for the North-South direction.
%   factorEW   - Scaling factor for the EW component.
%   factorNS   - Scaling factor for the NS component.
%
% Output:
%   srss       - SRSS-scaled acceleration spectrum.
%
% The calculation is performed vectorially:
%   srss = sqrt( (factorEW * spectrumEW).^2 + (factorNS * spectrumNS).^2 ).

    scaledEW = factorEW .* spectrumEW;
    scaledNS = factorNS .* spectrumNS;
    srss = sqrt(scaledEW.^2 + scaledNS.^2);
end

function displacement = computeSRSS_Displacement(periods, srssAcceleration, gravity)
% computeSRSS_Displacement computes the pseudo-displacement spectrum from
% the SRSS acceleration spectrum.
%
% Inputs:
%   periods          - Vector of vibration periods.
%   srssAcceleration - SRSS-scaled acceleration spectrum.
%   gravity          - Acceleration due to gravity (e.g., 9.81 m/s²).
%
% Output:
%   displacement     - Pseudo-displacement spectrum.
%
% The pseudo-displacement is computed using:
%   displacement = (srssAcceleration * gravity) / ( (2*pi./periods).^2 ).
% A small value check avoids division by zero.

    % Calculate angular frequencies (avoid division by zero for T==0)
    angularFrequencies = 2 * pi ./ periods;
    angularFrequencies(periods == 0) = Inf;
    displacement = srssAcceleration .* gravity ./ (angularFrequencies .^ 2);
end

function plotHandle = plotSpectrum(periods, values, spectrumType)
% plotSpectrum provides a unified way to plot spectra with different styles.
%
% Inputs:
%   periods      - Vector of periods.
%   values       - Spectrum values (acceleration or displacement).
%   spectrumType - String indicating the type of spectrum:
%                  'Reference'  => building code or reference data.
%                  'Average'    => average SRSS spectrum.
%                  otherwise    => individual earthquake spectrum.
%
% Output:
%   plotHandle   - Handle to the plotted graphic object.

    switch spectrumType
        case 'Reference'
            plotHandle = plot(periods, values, 'k', 'LineWidth', 2.5);
        case 'Average'
            plotHandle = plot(periods, values, '--k', 'LineWidth', 2);
        otherwise  % For individual spectra
            plotHandle = plot(periods, values, 'Color', [0.7 0.7 0.7],...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end