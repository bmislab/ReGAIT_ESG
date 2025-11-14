# ReGAIT
[![GitHub issues](https://img.shields.io/github/issues/sccn/bmsilab?color=%23fa251e&logo=GitHub)](https://github.com/bmislab/ReGAIT_ESG/issues)
![Twitter Follow](https://img.shields.io/twitter/follow/BMISLab?style=social)

# What is BMIsLAB?
The Brain-Machine Interface Systems Lab (BMIsLAB) is a team of researchers led by José María Azorín. Our work focuses on human-machine interaction through brain control in order to improve human capabilities in neural rehabilitation.
This repository is an open source signal processing environment for electrospinography (ESG) signals running on Matlab. This folder contains original Matlab functions from the BMISLAB that have been adapted to the context of the project ReGAIT.

# Signal Processing and Validation Toolbox
This MATLAB toolbox is developed to analyze database files with a focus on signal preprocessing and technical validation.

# Overview
The main function, process_signal, orchestrates the execution of preprocessing and validation procedures through a modular set of subfunctions:

# Preprocessing Functions
- ATS_filter_signals: Performs ECG denoising. Proper operation of this function requires the inclusion of Petersen’s cardiac artifact removal toolbox, available at: https://github.com/ime-luebeck/ecg-removal/releases/tag/1.1

- compute_noisy_outliers: Identifies outlier samples and detects noisy electrodes.

# Technical Validation Functions
- compute_metric_table: Generates a table containing key quality metrics including the ratio of null and noisy electrodes, the proportion of outlier samples, mean range and outlier amplitude.

- compute_correlation_time: Computes a time-based correlation matrix among electrodes and quaternion signals.

# Auxiliary Functions
- comb_filter, read_func, state_filter, state_filter_sample, and filter_setup: Perform bandpass filtering and apply multiple notch filters using a state-variable filter design.

# Optional Flags
To support flexible execution, the toolbox provides the following optional flags:

- flag_ATS: Runs preprocessing and applies ATS denoising.

- flag_outlier: Identifies noisy electrodes and outlier samples.

- flag_filtered_data: Replaces noisy and outlier data with NaNs in the matrix.

- flag_correlation_matrix: Computes inter-electrode correlation.

- flag_correlation_quaternions: Computes correlation between electrodes and quaternion data.

- flag_save: Saves the updated data_structure.

# Run

    process_signal.m   
    
to read the provided data sets, apply ECG-removal algorithms (or reload precomputed results) and execute a validation of the data.

# Requirements
- MATLAB R2018b or newer
- Petersen's ECG Artifact Removal Toolbox (see link above)

# License
This software is distributed under the MIT License. See below for full terms.

# MIT License
Copyright Copyright 2025, Brain-Machine Interfaces Systems Lab, 
Universidad Miguel Hernández de Elche
Author: Desirée I. Gracia, Eduardo Iáñez, Mario Ortiz, Jose M. Azorin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# In publications, please reference:

# Funding
This publication is part of grant PID2021-124111OB-C31, funded by MICIU/AEI/10.13039/501100011033 and by ERDF, EU. This research has been also supported by grant CIACIF/2022/108 funded by ``Consellería de Educación, Universidades y Empleo (Generalitat Valenciana)" and the European Social Fund, and grant PRE2022-103336 funded by MICIU/AEI/10.13039/501100011033 and by ESF+.




