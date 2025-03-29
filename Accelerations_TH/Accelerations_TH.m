% =========================================================================
% BUILDING SEISMIC RESPONSE ANALYSIS: VERTICAL ACCELERATION PROFILE VISUALIZATION
% =========================================================================
% Author: Victor Calderón (July 2019)
% Updated: Jefferson De la Cuba (February 2025)
% --------------------------------------------------------------------------
% Analyzes and visualizes peak floor accelerations across building levels during seismic events.
% 
% Inputs:
%   - Accel_Max_TH_X.txt: X-direction acceleration time history data (text file)
%   - Accel_Max_TH_Y.txt: Y-direction acceleration time history data (text file)
% 
% Outputs:
%   - Average_Acceleration_TH.png: Acceleration profile visualization (outputs/)

clear; close all; clc;  % Initialize clean workspace

%% ====================== DATA INGESTION ================================
% Configure platform-independent path handling
base_dir = fileparts(mfilename('fullpath'));  % Get current script location
datasets_folder = fullfile(base_dir, '..', 'datasets');
outputs_folder = fullfile(base_dir, '..', 'outputs');

% Validate directory structure
if ~isfolder(datasets_folder)
    error('Dataset directory not found: %s', datasets_folder);
end

% Ensure output directory existence
if ~isfolder(outputs_folder)
    mkdir(outputs_folder);
    fprintf('Created output directory: %s\n', outputs_folder);
end

% Import acceleration records
x_acceleration = load(fullfile(datasets_folder, 'Accel_Max_TH_X.txt'));
y_acceleration = load(fullfile(datasets_folder, 'Accel_Max_TH_Y.txt'));

%% ====================== VISUALIZATION ENGINE ==========================
% Initialize figure with publication-quality settings
figure_handle = figure('Units', 'centimeters', 'Position', [0 0 18 20],...
    'Color', 'w',...
    'PaperPositionMode', 'auto');

% Generate acceleration profiles
level_labels = {'Ground level'; 'Basement'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'};
num_levels = length(level_labels);

hold on;
plot(x_acceleration, 1:num_levels, '-ok', 'LineWidth', 1.2,...
    'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);
plot(y_acceleration, 1:num_levels, '-^k', 'LineWidth', 1.2,...
    'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 8);

%% ====================== VISUALIZATION FORMATTING ======================
% Configure axis system
ax = gca;
set(ax, 'YTick', 1:num_levels, 'YTickLabel', level_labels,...
    'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1.5,...
    'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.4);

axis([0 1.1 0.5 num_levels+0.5]);
xlabel('Peak Acceleration (g)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Building Level', 'FontSize', 14, 'FontWeight', 'bold');

% Add optimized legend
legend({'X Direction', 'Y Direction'},...
    'Location', 'northeast',...
    'Box', 'off',...
    'FontSize', 12);

%% ====================== OUTPUT MANAGEMENT ============================
% Generate filename with timestamp
output_filename = fullfile(outputs_folder, 'Average_Acceleration_TH.png');

% Save figure in multiple formats
print(figure_handle, output_filename, '-dpng', '-r600');  % High-res PNG
fprintf('Saved visualization: %s\n', output_filename);

% Clean up graphics memory
close(gcf);  % Close current figure after saving