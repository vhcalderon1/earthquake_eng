% =========================================================================
% EARTHQUAKE ACCELEROGRAM PROCESSING: NORMATIVE SPECTRUM SCALING AND ANALYSIS
% =========================================================================
% Author: Victor Calderón (January 2019)
% Updated: Jefferson De la Cuba (February 2025)
% -----------------------------------------------------
% Processes raw accelerograms from the 2010 Maule earthquake,
% scales horizontal components to match Peru's E.031 seismic
% design spectrum, and generates analysis outputs.
% 
% INPUTS:
% - Constitucion_2010cms2.txt: Raw triaxial accelerations (cm/s²)
% - Espectro_NormaE031.txt: Normative spectrum [Period (s), Sa (g)]
%
% OUTPUTS:
% - ScaledAccelerograms_Constitucion2010.xlsx: Time histories (cm/s²)
% - Accelerogram_*.png: Unscaled/scaled acceleration plots

clear; close all; clc;

%% ======================== PATH CONFIGURATION ============================
% Configure paths using file traversal
baseDir = fileparts(pwd);
outputDir = fullfile(baseDir, 'outputs');
datasetsDir = fullfile(baseDir, 'datasets');

% Create output directory if non-existent
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ===================== DATA IMPORT AND PROCESSING =======================
% Load earthquake recordings and normative spectrum
rawData = load(fullfile(datasetsDir, 'Constitucion_2010_cm_s2.txt'));
normSpec = load(fullfile(datasetsDir, 'E031_IsolatedBase_Spectrum.txt'));

% Physical constants and signal parameters
gravity = 981;      % cm/s²
dt = 0.005;         % Time step (s)
nsamples = size(rawData,1);
timeVector = (0:nsamples-1)*dt;

% Component extraction and unit conversion
ewAccel = rawData(:,1);     % East-West (cm/s²)
nsAccel = rawData(:,2);     % North-South (cm/s²)
vertAccel = rawData(:,3);   % Vertical (cm/s²)

%% ==================== ACCELEROGRAM VISUALIZATION ========================
% Generate and save unscaled accelerogram plots
plotComponents(timeVector, {ewAccel/gravity, nsAccel/gravity, vertAccel/gravity}, ...
    'Original Accelerograms', 'Constitucion Earthquake 2010, 8.8 Mw', outputDir);

%% ================== SPECTRAL SCALING PROCEDURE ==========================
% Spectral matching parameters
designPGA = normSpec(1,2);  % PGA from normative spectrum

% Calculate scaling factors
originalPGAs = [max(abs(ewAccel/gravity)), max(abs(nsAccel/gravity))];
scaleFactors = designPGA ./ originalPGAs;

% Apply scaling to horizontal components
ewScaled = ewAccel * scaleFactors(1);  % Scaled E-W (cm/s²)
nsScaled = nsAccel * scaleFactors(2);  % Scaled N-S (cm/s²)

%% ================== SCALED ACCELEROGRAM OUTPUTS =========================
% Generate and save scaled accelerogram plots
plotComponents(timeVector, {ewScaled/gravity, nsScaled/gravity}, ...
    'Scaled Accelerograms', 'Normative Spectrum Scaled Components', outputDir);

% Export scaled time histories
outputTable = table(timeVector', ewScaled, nsScaled, ...
    'VariableNames', {'Time_s', 'EW_cms2', 'NS_cms2'});
writetable(outputTable, fullfile(outputDir, 'ScaledAccelerograms.xlsx'));

%% ========================== HELPER FUNCTIONS ============================
function plotComponents(t, components, figName, ~, outputDir)
    % Unified plotting function with auto-saving
    fig = figure('Name', figName, 'Position',[100 100 800 900], 'Visible','off');
    
    for i = 1:length(components)
        subplot(length(components),1,i)
        [peakVal, peakTime] = getAbsolutePeak(components{i}, t(2)-t(1));
        
        plot(t, components{i}, 'k', 'LineWidth', 0.7)
        hold on
        plot(peakTime, peakVal, 'ro', 'MarkerSize', 6, 'LineWidth',1.5)
        plot([t(1) t(end)], [0 0], 'k', 'LineWidth',1.5)
        
        % Axis formatting
        ylim(max(abs(components{i}))*[-1.1 1.1])
        xlim([t(1) t(end)])
        set(gca, 'LineWidth',1.5, 'FontName','Times New Roman')
        
        % Annotation
        text(0.8*t(end), 0.8*max(ylim), ...
            sprintf('Peak: %.2f g', peakVal), ...
            'FontName','Times New Roman', 'FontSize',12)
    end
    
    % Save and close figure
    saveas(fig, fullfile(outputDir, [figName '.png']));
    close(fig)
end

function [absPeak, peakTime] = getAbsolutePeak(signal, dt)
    % Peak detection with temporal localization
    [maxVal, maxIdx] = max(signal);
    [minVal, minIdx] = min(signal);
    
    if abs(maxVal) >= abs(minVal)
        absPeak = maxVal;
        peakTime = (maxIdx-1)*dt;
    else
        absPeak = minVal;
        peakTime = (minIdx-1)*dt;
    end
end