# Seismic Analysis Toolkit for Structural Engineering

_A MATLAB-based repository for earthquake acceleration, drift, and spectral analysis_

## References
Calderon, V. H. (2021). “Feasibility of protecting with seismic isolation a limited-ductility-wall social housing.” (In Spanish) Pontifical Catholic University of Peru Digital Repository. Available at: https://tesis.pucp.edu.pe/items/033969e7-0b01-494f-ac01-6d2167bd67d4

Originally created by: Victor Calderon (November 2018 - February 2019 and August 2020 - February 2021)

Updated by: Jefferson de la Cuba (February 2025)

## What is this repository?
This repository contains a suite of MATLAB scripts designed for advanced seismic analysis of structures, including:
Building Acceleration Profile Visualization: Vertical acceleration profiles across building levels in X/Y directions.
Story Drift Analysis: Inter-story drift vs. height plots for earthquake resilience evaluation.
Earthquake Accelerogram Processing: Scaling raw accelerograms to normative design spectra (e.g., Peru’s E.031).
SRSS Spectrum Analysis: Combining orthogonal earthquake components using the Square Root of Sum of Squares (SRSS) method.
Multi-Event Spectrum Averaging: Comparing averaged earthquake spectra with building code requirements.
The tools are optimized for academic research, structural design validation, and earthquake engineering workflows.

## How does it work?
### Workflow Overview
#### Data Ingestion
- Load acceleration/drift data or normative spectra from standardized .txt files.
- Validate input formats (e.g., period vs. spectral acceleration pairs).

#### Core Computations
- Compute SRSS-combined spectra, pseudo-displacement, and peak accelerations.
- Scale accelerograms to match code-prescribed spectra (e.g., PGA alignment).
- Interpolate mismatched period grids for multi-event comparisons.

#### Visualization
- Generate publication-ready plots (Times New Roman fonts, consistent styling).
- Overlay earthquake records with normative spectra for compliance checks.

#### Output Export
- Save figures in .eps, .tiff, and .png formats.
- Export scaled data to .xlsx or .csv for further analysis.

## Why use this repository?
### Academic
Research Reproducibility: Preprocessing, scaling, and visualization workflows follow peer-reviewed methodologies.
Efficiency: Vectorized computations and automated interpolation ensure rapid analysis of large datasets.

### Key Advantages
Rigorous Spectral Scaling: Align recorded earthquake data with normative spectra for design-level comparisons.
Multi-Event Averaging: Statistically robust SRSS results across historical earthquakes (e.g., 2001 Arequipa, 2010 Concepción).
Transparent Workflow: Functions include error handling, dynamic axis scaling, and unit consistency checks.
