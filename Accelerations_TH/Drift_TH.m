% ==============================================================================================
% BUILDING_STORY_DRIFT_ANALYSIS - Vertical profile visualization of structural drift percentages
% ==============================================================================================
% Author: Victor Calderón (February 2021)
% Updated: Jefferson De la Cuba (February 2025)
% -----------------------------------------------------------
% This script analyzes interstory drift ratios in X/Y directions across building levels, 
% visualizing compliance with seismic performance criteria.
%
% Inputs (from external datasets directory):
%   - Drift_Max_TH_X.txt : Maximum X-direction drift ratios (decimal format)
%   - Drift_Max_TH_Y.txt : Maximum Y-direction drift ratios (decimal format)
%
% Outputs (to external outputs directory):
%   - Building_Drift_Analysis.png : Combined drift visualization figure
%

%% ======================= INITIALIZATION ==============================
clear; close all; clc;

%% ================== INPUT/OUTPUT CONFIGURATION ========================
% Path configuration using relative path traversal
basePath = fileparts(mfilename('fullpath')); % Get current script location
inputFolder = fullfile(basePath, '..', 'datasets');  % External datasets
outputFolder = fullfile(basePath, '..', 'outputs');  % Analysis outputs

% Validate directory structure
if ~exist(inputFolder, 'dir')
    error('Dataset directory not found: Ensure ../datasets exists');
end

% Create outputs directory with existence check
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output directory: %s\n', outputFolder);
end

%% ===================== DATA IMPORT PROCESSING =========================
% Define input file paths
yDriftFile = fullfile(inputFolder, 'Drift_Max_TH_Y.txt');
xDriftFile = fullfile(inputFolder, 'Drift_Max_TH_X.txt');

% Validate input files
if ~isfile(yDriftFile) || ~isfile(xDriftFile)
    error('Missing input files: Verify Drift_Max_TH_X/Y.txt in datasets');
end

% Load and convert drift ratios to percentages
driftYData = [NaN; NaN; load(yDriftFile)] * 100; % Y-direction (%)
driftXData = [NaN; NaN; load(xDriftFile)] * 100; % X-direction (%)

%% ================== BUILDING HEIGHT CALCULATION ========================
levelNames = {'Basement'; 'Ground'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'};
baseElevation = 2.5;    % Ground level height (m)
storyHeight = 5;        % Interstory height (m)

% Cumulative height calculation
storyHeights = (1:8)' * storyHeight + baseElevation;
totalHeightVector = [0; baseElevation; storyHeights]; % Full vertical profile

%% ====================== VISUALIZATION ENGINE ===========================
figure('Color', 'w', 'Units', 'centimeters', 'Position', [0 0 9 10]);
hold on;

% X-direction drift plot
plot(driftXData, totalHeightVector, '-ok', 'LineWidth', 0.9,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8 0.8 0.8],...
    'MarkerSize', 7, 'DisplayName', 'X-Direction');

% Y-direction drift plot
plot(driftYData, totalHeightVector, '-^k', 'LineWidth', 0.9,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8 0.8 0.8],...
    'MarkerSize', 7, 'DisplayName', 'Y-Direction');

%% =================== AXIS FORMATTING & LABELING ========================
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;
ax.LineWidth = 1.5;

% Vertical axis configuration
ax.YTick = totalHeightVector;
ax.YTickLabel = levelNames;
ax.YLim = [0 max(totalHeightVector)];
ylabel('Story Level', 'FontSize', 13)

% Horizontal axis configuration with dynamic scaling
validDriftData = [driftXData(:); driftYData(:)];
validDriftData = validDriftData(isfinite(validDriftData)); % Remove NaNs
maxDrift = max(validDriftData);
ax.XLim = [0 ceil(maxDrift*1.1)]; % 10% buffer above max drift
ax.XTick = linspace(0, ceil(maxDrift*1.1), 6); % 6 ticks across range
xlabel('Drift (%)', 'FontSize', 13)

%% =================== VISUAL ENHANCEMENTS ==============================
grid on;
ax.GridLineStyle = ':';
ax.GridColor = [0.2 0.2 0.2];
ax.GridAlpha = 0.3;
box off;

% Legend configuration
legend('Location', 'northeast', 'FontSize', 10, 'EdgeColor', 'none');

%% ===================== OUTPUT MANAGEMENT ==============================
outputFile = fullfile(outputFolder, 'Building_Drift_Analysis');
print(gcf, outputFile, '-dpng', '-r300'); % 300 DPI PNG export
close(gcf); % Prevent memory accumulation
fprintf('Analysis complete: Output saved to %s\n', outputFile);