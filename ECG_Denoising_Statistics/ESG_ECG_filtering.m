% This MATLAB script is designed to process all files within the database,
% focusing on evaluating ECG filtering algorithms performance:
% TS, ATS, HP200, HP200+ATS, % SWT, EKS2 and EMD.

% The main function, 'process_matrix', handles the signal processing.
% For correct operation of the ECG denoising algorithms, Petersen's
% Cardiac Artifact Removal Toolbox must be included.

% The subfunctions 'calculate_parameters' and 'build_results_table'
% obtain the performance metrics used to compare the performance of each
% algorithm and format them for later analysis, respectively.

% Auxiliary functions perform bandpass filtering and apply multiple notch
% filters using a state-variable filter design:
%   - 'comb_filter', 'read_func', 'state_filter', 'state_filter_sample' and 'filter_setup'.

% Note: In order to use 'process_matrix', it is necessary to include the
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

% Define types of processed data that will be generated.
typedata = {'rawData','dataTS','dataATS','dataHP200', 'dataATS_HP200Last','dataEKS2','dataSWT','dataEMD'};

% Initialize structure to store results.
results = table([], [], [], [], [], [], ...
    'VariableNames', {'Algorithm', 'Metric', 'Subject', 'Matrix', 'Repetition', 'Value'});


% Main Loop: Iterate through experiments, subjects and selected trials
for test = 1:2
    for subj = 1:numel(experiment(test).name_trial)
        trial_base = experiment(test).name_trial{subj};

        for trial = [1, 24] % Only first and last trials (Stage 1 trials)
            % Build trial file name
            trial_file = sprintf('%s_trial_%02d.mat', trial_base, trial);
            load(trial_file);         % Load the data structure

            % Process Brachial Matrix
            param = process_matrix(data_estructure.brachial_matrix, data_estructure.experiment_data, typedata, trial_file, 1);

            % Append to results table
            newRows = build_results_table(param, subj, test, trial);
            results = [results; newRows];


            if test == 2
                % Process Lumbar Matrix
                param = process_matrix(data_estructure.lumbar_matrix, data_estructure.experiment_data, typedata, trial_file, 2);

                % Append to results table
                newRows = build_results_table(param, subj, test, trial);
                results = [results; newRows];

            end

        end
    end
end

% Save results table
save('results_summary.mat','results');


% Helper Functions %

function param = process_matrix(electrode_data, experiment_data, typedata, ruta, matrix)
% Applies filtering, artifact removal and parameter computation to an
% electrode matrix. Saves processed results to disk.
%
% Inputs:
%   - electrode_data   : Structure with raw_data and metadata
%   - experiment_data               : Data about experiment as Sampling
%   frequency (Hz) and time vector (s)
%   - typedata         : Cell array of data type names
%   - ruta             : File path reference for saving
%   - matrix           : Identifier (1=Brachial, 2=Lumbar)
%
% Outputs:
%   - param            : Computed signal quality parameters

data = electrode_data.raw_data;
valid_electrodes = find(electrode_data.null_electrodes == 0)';  % Only process valid electrodes
fm = experiment_data.fm;


% Comb filtering for powerline and baseline removal
rawData = cell(1, size(data, 1));
for el = valid_electrodes
    rawData{1, el} = comb_filter(data(el, :), fm, 10, 500)';
end
rawData{1, size(data, 1)+1} = experiment_data.time_vector';  % Append time vector

% R-Peak Detection
algoRPeak = @(sig, time) [sig, peak_detection(sig, fm, time)'];
dataWithRPeaks = filter_signals(rawData, algoRPeak, electrodos_val, 0, size(data, 1)+1);

% High-Pass Filtering
hp20 =  @(sig, rpeaks) [butter_filt_stabilized(sig, 20, fm, 'high', true, 6), rpeaks];
dataWithRPeaksHP20 = filter_signals(dataWithRPeaks, hp20, electrodos_val, 2);

% ECG-Artifact Removal

% Template Subtraction (TS)
dataTS = filter_signals(dataWithRPeaksHP20, @template_subtraction, electrodos_val, 2);

% Adaptive Template Subtraction (ATS)
ats = @(signal, rpeaks) adaptive_template_subtraction(signal, rpeaks, fm);
dataATS = filter_signals(dataWithRPeaksHP20, ats, electrodos_val, 2);

% High-pass filter with 200 Hz cutoff frequency (HP200)
dataHP200 = filter_signals(dataWithRPeaks, ...
    @(signal) butter_filt_stabilized(signal, 200, fm, 'high', true, 3), electrodos_val);

% Adaptive Template Subtraction and subsequent high-pass filtering (ATS + HP200)
dataATS_HP200Last = filter_signals(dataATS, ...
    @(signal) butter_filt_stabilized(signal, 200, fm, 'high', true, 3), electrodos_val);

% Stationary Wavelet Transform (SWT)
swt = @(signal, rPeaks) swtden(signal, rPeaks, fm, 'h', 3, 'db2', 4.5);
dataSWT = filter_signals(dataWithRPeaks, swt, electrodos_val, 2);

% Second order Extended Kalman Smoother (EKS2)
ekf2 = @(signal, rpeaks) kalman_filter_2(signal, rpeaks, fm);
[~, dataEKS2] = filter_signals(dataWithRPeaks, ekf2, electrodos_val, 2);

% Empirical Mode Decomposition (EMD)
emd = @(signal) EMD(signal, 60, 1);
dataEMD = filter_signals(dataWithRPeaksHP20, emd, electrodos_val);

% Save Processed Data
data_filt = struct();
for i = 1:numel(typedata)
    fieldname = typedata{i};
    data_filt.(fieldname) = eval(fieldname);
end


if matrix==1
    save("Brachial_"+ruta(1:end-4), 'data_filt');
elseif matrix==2
    save("Lumbar_"+ruta(1:end-4), 'data_filt');
end
% Compute Signal Quality Parameters
param = NaN(7, length(typedata), size(data,1));
for i = 1:length(typedata)
    for el = electrodos_val
        % Compute parameters for each electrode and data type
        eval(sprintf(['param(:,i,el) = calculate_parameters(rawData{el}(15*fm+1:end), ' ...
            '%s{el}(15*fm+1:end), fm);'], typedata{i}));
    end
end

% Average across electrodes
param = mean(param, 3, 'omitnan');
end

function param = calculate_parameters(rawData, cleanData, fm)
% Computes quantitative quality metrics comparing raw vs. cleaned signals.
%
% Inputs:
%   - rawData          : Matrix with raw data
%   - cleanData        : Matrix with ECG denoised data
%   - fm               : Sampling frequency
%
% Outputs:
%   - param            : Computed signal quality parameters

%%% Revisar si necesario %%%
% Ensure column vectors
y0 = rawData(:);
y1 = cleanData(:);

% Relative Error (RE)
N = length(y0);
xdft0 = fft(y0);
xdft0 = xdft0(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft0).^2;
psdx(2:end-1) = 2 * psdx(2:end-1);
py0 = pow2db(psdx);

xdft1 = fft(y1);
xdft1 = xdft1(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft1).^2;
psdx(2:end-1) = 2 * psdx(2:end-1);
py1 = pow2db(psdx);

param(1) = sum((py0 - py1).^2) / sum(py0.^2);

% Signal-to-Noise Ratio (SNR)
param(2) = abs(10 * log10(std(y0) / std(y0 - y1)));

% Correlation Coefficient (CC)
param(3) = 100 * (sum(y0 .* y1) / sqrt(sum(y0.^2) * sum(y1.^2)));

% Spectral Distortion (SD) [0,20] and [20, 200] Hz
n_data=length(y1);
n_half_data=floor(n_data/2);
f1 = fm*(0:n_half_data)/n_data;
f_20 = find(f1 >= 20, 1, 'first');
f_200 = find(f1 <= 200, 1, 'last');
A0=abs(xdft0)/n_half_data;
A1=abs(xdft1)/n_half_data;
param(4) = trapz(A1(1:f_20)) / trapz(A0(1:f_20)) * 100;
param(5) = trapz(A1(f_20:f_200)) / trapz(A0(f_20:f_200)) * 100;

% Kurtosis Ratio (KR2)
[f, x] = ecdf(y1);
icdf_y1 = arrayfun(@(p) mean(x(max(find(f <= p, 1, 'last')): ...
    min(find(f >= p, 1, 'first')))), ...
    [0.975, 0.025, 0.75, 0.25]);
param(6) = (icdf_y1(1) - icdf_y1(2)) / (icdf_y1(3) - icdf_y1(4)) - 2.91;

% Variation of Kurtosis Ratio ($\Delta$KR2)
[f, x] = ecdf(y0);
icdf_y0 = arrayfun(@(p) mean(x(max(find(f <= p, 1, 'last')): ...
    min(find(f >= p, 1, 'first')))), ...
    [0.975, 0.025, 0.75, 0.25]);
KR20 = (icdf_y1(1) - icdf_y1(2)) / (icdf_y1(3) - icdf_y1(4)) - 2.91;
param(7) = ( param(6)-KR20)/KR20;
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

function T = build_results_table(param, subj, test, trial)
% Converts 3D parameter matrix into a results table.
%
% Inputs:
%   - param : 3D array [metrics x algorithms x electrodes]
%   - subj  : Subject index
%   - test  : Matrix index (1=Brachial, 2=Lumbar)
%   - trial : Repetition index (trial number)
%
% Output:
%   - T : Table with columns Algorithm, Metric, Subject, Matrix, Repetition, Value

    [metricCount, algorithmCount, elCount] = size(param);
    [Metric, Algorithm, Electrode] = ndgrid(1:metricCount, 1:algorithmCount, 1:elCount);

    T = table( ...
        Algorithm(:), ...
        Metric(:), ...
        repmat(subj, numel(Algorithm), 1), ...
        repmat(test, numel(Algorithm), 1), ...
        repmat(trial, numel(Algorithm), 1), ...
        param(:), ...
        'VariableNames', {'Algorithm','Metric','Subject','Matrix','Repetition','Value'} );
end
