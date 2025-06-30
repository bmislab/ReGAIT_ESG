% This MATLAB script is designed to process all files within the database
% focusing on preprocessing and technical validation tasks.

% The main function, 'process_signal', integrates two subfunctions
% responsible for different aspects of signal preprocessing:
%   - 'ATS_filter_signals': Performs ECG denoising. For correct operation,
%     Petersen’s cardiac artifact removal toolbox must be included.
%   - 'compute_noisy_outliers': Identifies outlier samples and detects noisy electrodes.

% Technical validation is performed using the following functions:
%   - 'compute_metric_table': Generates a table containing key signal quality metrics.
%   - 'compute_correlation_time': Computes correlation matrices over time.

% Auxiliary functions to perform bandpass filtering and apply multiple notch filters using a state-variable filter design.
%   - 'comb_filter', 'read_func', 'state_filter', 'state_filter_sample' and 'filter_setup'.

% The script includes several optional flags to allow flexible execution of
% specific processing steps:
%   - 'flag_ATS': Executes the preprocessing pipeline and applies ATS denoising.
%   - 'flag_outlier': Detects outlier samples and noisy electrodes.
%   - 'flag_filtered_data': Updates the data matrix by assigning NaN to
%     values identified as outliers or noisy.
%   - 'flag_correlation_matrix': Computes correlations among electrode signals.
%   - 'flag_correlation_quaternions': Computes correlations between electrode
%     signals and quaternion data.
%   - 'flag_save': Saves the updated 'data_structure' variable.

% Note: In order to use 'ATS_filter_signals', it is necessary to include the
% cardiac artifact removal toolbox from:
% https://github.com/ime-luebeck/ecg-removal/releases/tag/1.1

% Copyright 2025, Brain-Machine Interfaces Systems Lab, 
% Universidad Miguel Hernández de Elche
% Desirée I. Gracia, Eduardo Iáñez, Mario Ortiz, Jose M. Azorin

% This software is provided under the MIT License:
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.

% Define Trial Names for Each Experiment %
% Each experiment contains a set of trial identifiers, which are strings referencing data files.
experiment(1).name_trial = {'S01_20231114', 'S02_20231115', 'S03_20231204', 'S04_20231205', 'S05_20231218', 'S06_20231219', 'S11_20250220'};
experiment(2).name_trial = {'S05_20240507', 'S03_20240507', 'S01_20240523', 'S02_20240611', 'S07_20240627', 'S07_20241015', 'S10_20250211'};
experiment(3).name_trial = {'S08_20250210', 'S09_20250210'};

% Initialize Result Structures %
% Pre-allocate result matrices with NaNs to hold metrics and correlations
results.metrics.brachial = NaN(5,7,3,24);
results.metrics.lumbar = NaN(5,7,3,24);
results.metrics.linear = NaN(5,7,3,24);
results.correlation.test1 = NaN(68,68,7,24);
results.correlation.test2 = NaN(68,68,7,24);
results.correlation.test3 = NaN(68,68,7,24);

% Define Processing Flags %
% These boolean flags control which parts of the pipeline are executed
flag_ATS = true;                     % Apply Adaptive Template Subtraction (ECG removal)
flag_outiler = true;                 % Detect and mask outliers
flag_filtered_data = true;           % Generate filtered data matrix
flag_correlation_matrix = true;      % Compute electrode correlation metrics
flag_correlation_quaternions = true; % Include quaternions in correlation
flag_save = true;                    % Save updated data structure to disk

% Main Loop: Process Trials for each subject of each gorup of experiments %
for test = 1:3
    for subj = 1:numel(experiment(test).name_trial)
        trial_base = experiment(test).name_trial{subj};

        for trial = 1:24
            % Build trial file name
            trial_file = sprintf('%s_trial_%02d.mat', trial_base, trial);

            if exist(trial_file, 'file')  % Only proceed if the file exists
                load(trial_file);         % Load the data structure

                % Truncate quaternion data to match brachial matrix length
                data_estructure.quaternions.data(:, size(data_estructure.brachial_matrix.filtered_data, 2) + 1:end) = [];

                % Preprocess and compute brachial matrix metrics
                data_estructure.brachial_matrix = process_signal(data_estructure.brachial_matrix, ...
                    data_estructure.experiment_data, flag_ATS, flag_outiler, flag_filtered_data);

                results.metrics.brachial(:, subj, test, trial) = ...
                    compute_metric_table(data_estructure.brachial_matrix, ...
                    data_estructure.experiment_data.sample_frequency, ...
                    flag_correlation_matrix);

                switch test
                    case 1
                        % Correlate brachial matrix and quaternion signals
                        results.correlation.test1(:, :, subj, trial) = ...
                            compute_correlation_time(data_estructure.brachial_matrix.filtered_data, ...
                            [], ...
                            data_estructure.quaternions.data, ...
                            data_estructure.experiment_data.sample_frequency, ...
                            flag_correlation_quaternions);

                    case 2
                        % Preprocess and compute lumbar matrix metrics and correlation
                        data_estructure.lumbar_matrix = ...
                            process_signal(data_estructure.lumbar_matrix, ...
                            data_estructure.experiment_data, ...
                            flag_ATS, flag_outiler, flag_filtered_data);

                        results.metrics.lumbar(:, subj, test, trial) = ...
                            compute_metric_table(data_estructure.lumbar_matrix, ...
                            data_estructure.experiment_data.sample_frequency, ...
                            flag_correlation_matrix);

                        results.correlation.test2(:, :, subj, trial) = ...
                            compute_correlation_time(data_estructure.brachial_matrix.filtered_data, ...
                            data_estructure.lumbar_matrix.filtered_data, ...
                            data_estructure.quaternions.data, ...
                            data_estructure.experiment_data.sample_frequency, ...
                            flag_correlation_quaternions);

                    case 3
                        % Preprocess and compute linear array metrics and correlation
                        data_estructure.linear_array = ...
                            process_signal(data_estructure.linear_array, ...
                            data_estructure.experiment_data, ...
                            flag_ATS, flag_outiler, flag_filtered_data);

                        results.metrics.linear(:, subj, test, trial) = ...
                            compute_metric_table(data_estructure.linear_array, ...
                            data_estructure.experiment_data.sample_frequency, ...
                            flag_correlation_matrix);

                        results.correlation.test3(:, :, subj, trial) = ...
                            compute_correlation_time(data_estructure.brachial_matrix.filtered_data, ...
                            data_estructure.linear_array.filtered_data, ...
                            data_estructure.quaternions.data, ...
                            data_estructure.experiment_data.sample_frequency, ...
                            flag_correlation_quaternions);
                end

                % Save updated data if required
                if flag_save
                    save(trial_file, 'data_estructure');
                end
            end
        end
    end
end

% Helper Functions %

function electrode_data = process_signal(electrode_data, experiment_data, flag_ATS, flag_outiler, flag_filtered_data)
    % Applies signal processing steps to an electrode matrix.
    %
    % Inputs:
    % - electrode_data: Struct containing raw_data and metadata
    % - experiment_data: Struct with time vector and sample frequency
    % - flag_ATS: Enable/disable  ECG artifact removal (Adaptive Template Subtraction)
    % - flag_outiler: Enable/disable  noisy outlier detection
    % - flag_filtered_data: Enable/disable  storing filtered result
    %
    % Output:
    % - electrode_data: Struct updated with processed fields
    
    if flag_ATS
        % Apply ECG artifact removal via adaptive template subtraction
        electrode_data = ATS_filter_signals(electrode_data, experiment_data.time_vector, experiment_data.sample_frequency);
    end
    
    % Detect and mask noisy electrodes and outlier segments
    electrode_data = compute_noisy_outliers(electrode_data, experiment_data.sample_frequency, flag_outiler, flag_filtered_data);
end


function electrode_data = ATS_filter_signals(electrode_data, time_vector, fm)
    % Applies ECG artifact removal using Adaptive Template Subtraction
    % (Petersen's Cardiac Artifact Removal Toolbox needed)
    %
    % Inputs:
    % - electrode_data: Struct containing raw_data and null_electrodes
    % - time_vector: Time vector for signal alignment
    % - fm: Sampling frequency
    %
    % Output:
    % - electrode_data: Struct with field ATS_data (artifact-reduced signal)
    
    data = electrode_data.raw_data;
    valid_electrodes = find(electrode_data.null_electrodes == 0)';  % Only process valid electrodes
    
    % Comb filtering for powerline and baseline removal
    rawData = cell(1, size(data, 1));
    for el = valid_electrodes
        rawData{1, el} = comb_filter(data(el, :), fm, 10, 500)';
    end
    rawData{1, size(data, 1)+1} = time_vector';  % Append time vector
    
    % R-Peak Detection
    algoRPeak = @(sig, time) [sig, peak_detection(sig, fm, time)'];
    % Variables: Matrix data - function to detect R peak - Electrodes to
    % filter (no nulls) - Auxiliar channel (none) - Auxiliar recording (time vector)
    dataWithRPeaks = filter_signals(rawData, algoRPeak, valid_electrodes, 0, size(data, 1)+1);
    
    % High-Pass Filter (20 Hz)
    hp20 =  @(sig, rpeaks) [butter_filt_stabilized(sig, 20, fm, 'high', true, 6), rpeaks];
    % Variables: Matrix data - function to high filter 20 Hz - Electrodes to filter (no nulls) -
    % Auxiliar channel (time R peak detections) - Auxiliar recording (none)
    dataWithRPeaksHP20 = filter_signals(dataWithRPeaks, hp20, valid_electrodes, 2);
    
    % ECG-Artifact Removal - Adaptative Template Subtraction (ATS)
    ats = @(signal, rpeaks) adaptive_template_subtraction(signal, rpeaks, fm);
    % Variables: Matrix data - ATS algorithm - Electrodes to filter (no nulls) -
    % Auxiliar channel (time R peak detections) - Auxiliar recording (none)
    ATS_data = filter_signals(dataWithRPeaksHP20, ats, valid_electrodes, 2);
    
    % Store result in structured format
    electrode_data.ATS_data = NaN(size(data));
    for el = valid_electrodes
        electrode_data.ATS_data(el,:) = ATS_data{el};
    end
end

function electrode_data = compute_noisy_outliers(electrode_data, fm, flag_outiler, flag_filtered_data)
    % Identifies outlier segments and noisy electrodes based on amplitude
    % statistics and then replaces them with NaNs.
    %
    % Inputs:
    % - electrode_data: Struct with ATS-filtered data
    % - fm: Sampling frequency (used to determine segment size)
    % - flag_outiler: Enable/disable  computation and mask outliers
    % - flag_filtered_data: Enable/disable  update .filtered_data with NaNs
    %
    % Output:
    % - electrode_data: Updated with outlier mask, noisy electrode info and optionally filtered data
    
    if flag_outiler
        % Parameters
        seg_size = fm; % Segment size in samples
        n_segments = floor(size(electrode_data.ATS_data, 2) / seg_size); % Total number of segments
    
        % Reshape electrode data into 3D segments: [electrode, time, segment]
        data_segments = reshape(electrode_data.ATS_data(:, 1:n_segments * seg_size), size(electrode_data.ATS_data, 1), seg_size, []);
    
        % Compute median and MAD across time for each electrode
        median_segments = median(data_segments, 2, 'omitnan');
        mad_segments = mad(data_segments, 1, 2);
    
        % Z-scores per sample within segments
        z_scores = (data_segments - median_segments) ./ (1.4826 * mad_segments);
    
        % Identify and mask outliers
        outlier_mask = abs(z_scores) > 3;
        data_segments(outlier_mask) = NaN;
    
        % Range (amplitude) for each segment, per electrode
        segment_ranges = squeeze(range(data_segments, 2));
        mean_range = mean(segment_ranges, 2, 'omitnan') / 2;
    
        % Identify noisy electrodes via MAD-based z-score
        z_scores_el = (mean_range - median(mean_range, 'omitnan')) ./ (1.4826 * mad(mean_range, 1));
        noisy_electrodes = abs(z_scores_el) > 3;
        data_segments(noisy_electrodes,:,:) = NaN;
    
        % Store results
        electrode_data.outlier_mask = reshape(outlier_mask, size(electrode_data.ATS_data, 1), []);
        electrode_data.noisy_electrodes=noisy_electrodes;
    
        if flag_filtered_data
            electrode_data.filtered_data = reshape(data_segments, size(electrode_data.ATS_data, 1), []);
        end
    
    elseif flag_filtered_data
        % Apply existing outlier mask (if available)
        electrode_data.filtered_data = electrode_data.ATS_data;
        electrode_data.filtered_data(electrode_data.outlier_mask)=NaN;
        electrode_data.filtered_data(electrode_data.noisy_electrodes,:) = NaN;
    end
end

function [metric_table] = compute_metric_table(electrode_data, fm, flag_correlation_matrix)
    % Computes quality metrics for an electrode matrix.
    %
    % Inputs:
    % - electrode_data: Struct with filtered_data and outlier metadata
    % - fm: Sampling frequency (used to determine segment size)
    % - flag_correlation_matrix: Enable/disable metric computation
    %
    % Output:
    % - metric_table: [null_electrodes; mean_amp; noisy_electrodes; outlier_ratio; outlier_amp]
    
    if flag_correlation_matrix
        % Parameters
        seg_size = fm;
        n_segments = floor(size(electrode_data.filtered_data, 2) / seg_size);
    
        % Reshape to 3D: [electrode, time, segment]
        data_segments = reshape(electrode_data.filtered_data(:, 1:n_segments * seg_size), size(electrode_data.filtered_data, 1), seg_size, []);
    
        % Compute average amplitude range
        segment_ranges = squeeze(range(data_segments, 2));
        mean_range = mean(mean(segment_ranges, 2, 'omitnan') / 2, 'omitnan');
    
        % Percentage null, noisy, outlier stats
        null_electrodes = mean(electrode_data.null_electrodes);
        noisy_electrodes = mean(electrode_data.noisy_electrodes);
        outlier_values = electrode_data.ATS_data(electrode_data.outlier_mask);
        outlier_amplitude = mean(abs(outlier_values), 'omitnan');
        outlier_ratio = sum(electrode_data.outlier_mask(:)) / sum(electrode_data.filtered_data(:) ~= 0);
    
        % Format result
        metric_table = [null_electrodes; mean_range; noisy_electrodes; outlier_ratio; outlier_amplitude];
    end
end


function correlation = compute_correlation_time(electrode_data1, electrode_data2, quaternionsdata, fm, flag_correlation_quaternions)
    % Computes correlation matrix between electrode matrix and quaternion data.
    %
    % Inputs:
    % - electrode_data1: Electrode Matrix [electrode x time]
    % - electrode_data2: Electrode Matrix [electrode x time], if present
    % - quaternionsdata: Matrix [4 x time]
    % - fm: Sampling frequency (used to determine window size)
    % - flag_correlation_quaternions: Enable/disable correlation computation
    %
    % Output:
    % - correlation: Mean correlation matrix across windows
    
    if flag_correlation_quaternions
        % Truncate quaternion data to match electrode matrix length
        quaternionsdata(:, size(electrode_data1, 2) + 1:end) = [];
    
        % Combine all signals into one matrix (samples x features)
        combined_data = [electrode_data1', electrode_data2', quaternionsdata'];
    
        % Initialize output: [features x features x time_windows]
        corr_matrix = NaN(size(combined_data, 2), size(combined_data, 2), floor(size(combined_data, 1)/fm) - 1);
    
        % Compute correlation over 1-second sliding windows
        for time = 1:fm:(size(combined_data, 1) - fm)
            corr_matrix(:, :, (time - 1) / fm + 1) = corrcoef(combined_data(time:time + fm - 1, :), 'rows', 'pairwise');
        end
    
        % Average correlation matrices over time
        correlation = mean(corr_matrix, 3, 'omitnan');
    end
end

function y_filtered = comb_filter(data, fs, fc1, fc2)
    % Applies a comb filter (high-pass, low-pass and notch filters) 
    %
    % Inputs:
    % - data: Signal matrix [channels x time]
    % - fs: Sampling frequency (Hz)
    % - fc1: High-pass filter cutoff frequency (Hz)
    % - fc2: Low-pass filter cutoff frequency (Hz)
    %
    % Output:
    % - y_filtered: Filtered signal
    
    % Parameters
    shift = 0.5;       % Shift step in seconds
    epoch_size = 1;    % Epoch length in seconds
    [n_ch, n_t] = size(data);
    y = detrend(data')';  % Detrend each channel
    shift_length = shift * fs; % Calculate number of samples per shift

    % Setup high-pass filter (2nd order)
    [A_h, B_h, C_h, D_h] = filter_setup(fs, fc1, 2, 'high');
    Xnn_h = repmat({zeros(size(A_h,1),1)}, 1, n_ch);

    % Setup low-pass filter (2nd order)
    [A_l, B_l, C_l, D_l] = filter_setup(fs, fc2, 2, 'low');
    Xnn_l = repmat({zeros(size(A_l,1),1)}, 1, n_ch);

    % Setup notch filters for 50 Hz harmonics up to 500 Hz (3rd order)
    notch_freqs = 50:50:500;
    notch_filters = cell(length(notch_freqs), 1);
    Xnn_notch = cell(length(notch_freqs), n_ch);
    for i = 1:length(notch_freqs)
        [A, B, C, D] = filter_setup(fs, [notch_freqs(i)-2, notch_freqs(i)+2], 3, 'stop');
        notch_filters{i} = {A, B, C, D};
        Xnn_notch(i, :) = repmat({zeros(size(A,1),1)}, 1, n_ch);
    end

    % Apply all filters to signal
    y_filtered = read_func(y, shift, n_t, fs, shift_length, ...
        A_h, B_h, C_h, D_h, Xnn_h, ...
        A_l, B_l, C_l, D_l, Xnn_l, ...
        notch_filters, Xnn_notch);
end

function y_filtered = read_func(y, shift, n_t, fs, shift_length, ...
    A_h, B_h, C_h, D_h, Xnn_h, ...
    A_l, B_l, C_l, D_l, Xnn_l, ...
    notch_filters, Xnn_notch)
    % Segments input signal and applies high-pass, low-pass and notch filters.
    %
    % Inputs:
    % - y: Signal matrix [channels x time]
    % - shift: Shift step in seconds
    % - n_t: Total number of time samples
    % - fs: Sampling frequency
    % - shift_length: Number of samples in each shift segment
    % - A_*, B_*, C_*, D_*: State-space matrices for filters
    % - Xnn_*: Filter states for each channel
    % - notch_filters: Cell array of state-space matrices for notch filters
    % - Xnn_notch: State vectors for notch filters
    %
    % Output:
    % - y_filtered: Filtered signal
    
    % Parameters
    y_filtered = zeros(size(y));
    i_1 = 1:fs*shift:(n_t - fs*shift + 1);
    i_2 = fs*shift:fs*shift:n_t;
    n_shifts = length(i_1);

    for k = 1:n_shifts
        input = y(:, i_1(k):i_2(k));

        % Apply high-pass filter
        [output, Xnn_h] = state_filter(input, A_h, B_h, C_h, D_h, Xnn_h, shift_length);
        
        % Apply low-pass filter
        [output, Xnn_l] = state_filter(output, A_l, B_l, C_l, D_l, Xnn_l, shift_length);

        % Apply each notch filter
        for i = 1:length(notch_filters)
            [A, B, C, D] = deal(notch_filters{i}{:});
            [output, Xnn_notch(i, :)] = state_filter(output, A, B, C, D, Xnn_notch(i, :), shift_length);
        end
        
        % Store filtered segment
        y_filtered(:, i_1(k):i_2(k)) = output;
    end
end

function [output, Xnn] = state_filter(input, A, B, C, D, Xnn, shift_length)
    % Applies a state-space filter over multiple time steps
    %
    % Inputs:
    % - input: Input signal [channels x time]
    % - A, B, C, D: State-space matrices
    % - Xnn: Current state vectors for each channel
    % - shift_length: Number of time steps to filter
    %
    % Output:
    % - output: Filtered signal [channels x time]
    % - Xnn: Updated state vectors

    % Parameters
    [n_ch, ~] = size(input);
    output = zeros(n_ch, shift_length);
    
    for m = 1:shift_length
        [buffer, Xnn] = state_filter_sample(input(:, m), A, B, C, D, Xnn);
        output(:, m) = buffer;
    end
end

function [output, Xnn] = state_filter_sample(input, A, B, C, D, Xnn)
    % Applies one step of a state-space filter to input.
    %
    % Inputs:
    % - input: Current sample [channels x 1]
    % - A, B, C, D: State-space matrices
    % - Xnn: Current state vectors for each channel
    %
    % Output:
    % - output: Filtered output [channels x 1]
    % - Xnn: Updated state vectors

    % Parameters
    n_ch = size(input, 1);
    output = zeros(n_ch, 1);

    for k = 1:n_ch
        % State-space update
        Xn1 = A * Xnn{k} + B * input(k);
        output(k) = C * Xnn{k} + D * input(k);
        Xnn{k} = Xn1;
    end
end

function [A, B, C, D] = filter_setup(fs, f_cut, filt_order, filter_type)
    % Designs Butterworth filter as a state-space system
    %
    % Inputs:
    % - fs: Sampling frequency (Hz)
    % - f_cut: Cutoff frequency or range (Hz) (can be scalar or 2-element vector)
    % - filt_order: Order of the Butterworth filter
    % - filter_type: Type of filter ('low', 'high', 'stop', etc.)
    %
    % Outputs:
    % - A, B, C, D: State-space matrices of the filter
    
    Wn = f_cut / (fs / 2);  % Normalize cutoff frequency
    [A, B, C, D] = butter(filt_order, Wn, filter_type);
end
