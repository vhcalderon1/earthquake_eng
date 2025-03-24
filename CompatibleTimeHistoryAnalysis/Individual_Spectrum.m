% =================================================
% Seismic Spectrum Analysis and Visualization Tool
% =================================================
% Author: Victor Calderón (August 20120)
% Updated: Jefferson De la Cuba (February 2025)
% -----------------------------------------------
% Processes earthquake acceleration spectra and generates 
% standardized comparison plots against NTP E.031 specifications.
% Performs SRSS scaling of multi-directional components and
% converts acceleration spectra to displacement spectra.
%
% INPUTS:
% - Text files containing spectral data:
%        * EspectroE031_BaseAislada.txt (NTP standard spectrum)
%        * Concepcion2010_g_EW_espectro.txt (EW accelerogram)
%        * Concepcion2010_g_NS_espectro.txt (NS accelerogram)
% - Gravity constant: 981 cm/s²
% - Scaling factors: 0.9 (EW), 0.44 (NS)
%
% OUTPUTS:
% - ScaledAccelerationSpectrum.png
% - DisplacementSpectrum.png
%


%% ======================= INITIALIZATION ================================
clear; close all; clc;

%% ======================= PARAMETERS & PATHS ============================
% Global constants
GRAVITY_CM_PER_S2 = 981;       % Acceleration due to gravity in cm/s²

% Configure cross-platform path handling
INPUT_FOLDER  = fullfile('..', 'datasets');   % Parent directory for inputs
OUTPUT_FOLDER = fullfile('..', 'outputs');    % Parent directory for outputs

% Create outputs directory with existence check
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
    fprintf('Created output directory: %s\n', OUTPUT_FOLDER);
end

%% ===================== LOAD INPUT DATA =================================
% Load standard (NTP E.031) base isolated spectrum 
% (Assume file format: first column = period, second column = S_a/g)
standardSpectrum = load(fullfile(INPUT_FOLDER, 'E031_IsolatedBase_Spectrum.txt'));

% Load Chilean accelerograms (2010 Concepcion earthquake)
% (Assume file format: first column = period, second column = acceleration in g)
ewAcceleration = load(fullfile(INPUT_FOLDER, 'Concepcion2010_g_EW_Spectrum.txt'));
nsAcceleration = load(fullfile(INPUT_FOLDER, 'Concepcion2010_g_NS_Spectrum.txt'));

%% ====================== DATA PROCESSING ================================
% Define scaling factors for the directional components
scaleFactor_EW = 0.9;    % East-West component scaling factor
scaleFactor_NS = 0.44;   % North-South component scaling factor

% Compute the SRSS-scaled acceleration spectrum from the earthquake records.
% (Using the second code's algorithm, but in vectorized form.)
scaledSRSSAcceleration = computeSRSS(...
    ewAcceleration(:,2), nsAcceleration(:,2), scaleFactor_EW, scaleFactor_NS);

%% ==================== ACCELERATION SPECTRUM PLOT =======================
figure;
hold on;

% Plot the standard (NTP E.031) acceleration spectrum in light gray.
hStandard = plot(standardSpectrum(:,1), standardSpectrum(:,2), ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 2);

% Plot the computed SRSS earthquake acceleration spectrum in black.
hEarthquake = plot(ewAcceleration(:,1), scaledSRSSAcceleration, 'k', 'LineWidth', 1.5);

% Configure axes and appearance.
axis([0 5 0 2.0]);
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'FontName', 'Times New Roman');
box off;
xticks(0:1:5);
yticks(0:0.5:2.0);
set(gca, 'YTickLabel', num2str(get(gca, 'ytick')', '%.1f'));
set(gca, 'XTickLabel', []);  % Remove x-axis labels for this plot

% Add axis labels and legend.
ylabel('S_a/g', 'FontSize', 13, 'FontName', 'Times New Roman');
legend([hStandard, hEarthquake], {'NTP E.031 2019', 'Scaled Spectrum'}, ...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Box', 'off');

% Export the acceleration spectrum figure.
outputFilenameAccel = fullfile(OUTPUT_FOLDER, 'ScaledAccelerationSpectrum');
print(outputFilenameAccel, '-dpng', '-r300');
close(gcf);

%% ==================== DISPLACEMENT SPECTRUM CALCULATION ================
% Compute the displacement spectrum for the standard spectrum.
standardDisplacement = computeDisplacementSpectrum(...
    standardSpectrum(:,1), standardSpectrum(:,2), GRAVITY_CM_PER_S2);

% Compute the displacement spectrum for the earthquake (scaled SRSS) data.
earthquakeDisplacement = computeDisplacementSpectrum(...
    ewAcceleration(:,1), scaledSRSSAcceleration, GRAVITY_CM_PER_S2);

%% ==================== DISPLACEMENT SPECTRUM PLOT =======================
figure;
hold on;

% Plot the standard displacement spectrum in light gray.
hStandardDisp = plot(standardSpectrum(:,1), standardDisplacement, ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 2);

% Plot the earthquake displacement spectrum in black.
hEarthquakeDisp = plot(ewAcceleration(:,1), earthquakeDisplacement, ...
    'k', 'LineWidth', 1.5);

% Configure axes and appearance.
axis([0 5 0 60]);
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'FontName', 'Times New Roman');
box off;
xticks(0:1:5);
yticks(0:10:60);
set(gca, 'YTickLabel', num2str(get(gca, 'ytick')', '%.0f'));
set(gca, 'XTickLabel', num2str(get(gca, 'xtick')', '%.0f'));

% Add axis labels and legend.
xlabel('Period T (s)', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('S_d (cm)', 'FontSize', 13, 'FontName', 'Times New Roman');
legend([hStandardDisp, hEarthquakeDisp], {'NTP E.031 2019', 'Scaled Spectrum'}, ...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Location', 'best','Box', 'off');

% Export the displacement spectrum figure.
outputFilenameDisp = fullfile(OUTPUT_FOLDER, 'DisplacementSpectrum');
print(outputFilenameDisp, '-dpng', '-r300');
close(gcf);

%% ==================== LOCAL FUNCTIONS =============================

function SRS_Acceleration = computeSRSS(Spectrum_EW, Spectrum_NS, factor_EW, factor_NS)
    % computeSRSS computes the Square Root of the Sum of Squares (SRSS)
    % of two orthogonal earthquake acceleration spectra.
    %
    % Inputs:
    %   Spectrum_EW - Acceleration spectrum in the East-West direction (vector)
    %   Spectrum_NS - Acceleration spectrum in the North-South direction (vector)
    %   factor_EW   - Scaling factor for the EW spectrum
    %   factor_NS   - Scaling factor for the NS spectrum
    %
    % Output:
    %   SRS_Acceleration - Combined SRSS acceleration spectrum (vector)
    
    % Ensure the inputs are column vectors.
    Spectrum_EW = Spectrum_EW(:);
    Spectrum_NS = Spectrum_NS(:);
    
    % Vectorized computation of SRSS acceleration.
    SRS_Acceleration = sqrt((Spectrum_EW * factor_EW).^2 + (Spectrum_NS * factor_NS).^2);
end

function displacementSpectrum = computeDisplacementSpectrum(periods, accelerations, gravity)
    % computeDisplacementSpectrum converts an acceleration spectrum to a
    % displacement spectrum using the relation:
    %   S_d = (S_a * gravity * T^2) / (4*pi^2)
    %
    % Inputs:
    %   periods       : Vector of vibration periods (T) in seconds.
    %   accelerations : Spectral accelerations (S_a/g) corresponding to each period.
    %   gravity       : Acceleration due to gravity in cm/s².
    %
    % Output:
    %   displacementSpectrum : Computed displacement spectrum in cm.
    
    % Ensure inputs are column vectors.
    periods = periods(:);
    accelerations = accelerations(:);
    
    % Compute displacement spectrum (vectorized).
    displacementSpectrum = (accelerations .* gravity .* periods.^2) / (4 * pi^2);
    
    % Handle potential division-by-zero if any period is zero.
    displacementSpectrum(periods == 0) = 0;
end