% Author: Omer Hakan Yaran
% New Jersey Institute of Technology, 2021
% Advisor: Moshe Kam

%% [1] Obtaining and Classifying the Data
% Window size: how many samples in a window
% Step size: distance between start of each iterations

clc

% Loading .mat file 
load('processed_sensors_data - sub2.mat');
load('sub2_labeled_v2.mat')

% sampling frequency
fs = 100;

% capturing workpace data
workspaceVar = who;
Var = eval(workspaceVar{1});

% categorizing the data
ecg_data = Var(:,1);
g_ppg_data = Var(:,9);


%% [2] First Generation PPG HR Measurement
% Calculates the heart rate using PPG

% HR 42-300
g_ppg_filtered = bandpass(g_ppg_data, [0.5 10], fs, Steepness = 0.99);

% calculate time wrt frequency
time = linspace(0,length(ecg_data)/fs, length(ecg_data));


%% [3] Create Windows

% Window size
ws = 800;

% Step size
step_size = 100;

ecg_length = length(ecg_data);

% Max step
%max_step = ceil(ecg_length/ws);
max_step = floor((ecg_length-ws)/step_size);

% Preallocations for loops
ecg_hr = zeros(1, max_step-1);
scatter_data = zeros(max_step-1, 3);

% START OF FOR LOOP
for i = 1:max_step-1

    % Current step, see the variable max_step for the limit
    step = i;
    
    % Find Windows for ECG
    %st=(step-1)*ws+1;
    st = (step-1)*step_size+1;
    
    ecg_win = ecg_data(st:st+ws);
    ppg_win = g_ppg_filtered(st:st+ws);
   
    
    % Find windows for labels
    labels_wini = (st < labels(:,1))&(labels(:,1) < st+ws);
    labels_win = labels(labels_wini,:);
    
    
    % Average HR calculator ===============================================
    
    % This part is also essential for plotting
    inxi = labels_win(:,1);
    inx = inxi-st;
    val = labels_win(:,2);

    % Calculate mean of difference between each peaks in a window
    hr = 1/mean(diff(inx))*fs;
    ecg_hr(i) = hr*60;
    
    
    % calculate FFTs for a window =========================================

     [ecg_fft_abs, fax_ecg, ~ ] = FFT_ZeroPadded(ecg_win, 100, fs, 0);
     [ppg_fft_abs, fax_ppg, ~ ] = FFT_ZeroPadded(ppg_win, 100, fs, 0);

     ppg_fft_abs = abs(ppg_fft_abs);
%     ecg_fft_abs = abs(FFT_ZeroPadded(ecg_win, 100, 100, 0));
%     ppg_fft_abs = abs(FFT_ZeroPadded(ppg_win, 100, 100, 0));


    % Find frequency axises ===== =========================================

    % PPG
     ppg_fft_abs_len = length(ppg_fft_abs);
%     fax_ppg = fs*(0:(ppg_fft_abs_len/2))/ppg_fft_abs_len;
     fax_ppg_max = floor(ppg_fft_abs_len/2)+1;
% 
%     % ECG
     ecg_fft_abs_len = length(ecg_fft_abs);
%     fax_ecg = fs*(0:(ecg_fft_abs_len/2))/ecg_fft_abs_len;
     fax_ecg_max = floor(ecg_fft_abs_len/2)+1;

    
    % HR Estimation Algorithm =============================================
    
    % Define a space for searching HR
    est_ppg_mid = ecg_hr(i)/60;
    est_ppg_low = est_ppg_mid - 0.05;
    est_ppg_high = est_ppg_mid + 0.05;
    
    % First column is frequencies and the second column is related FFT vals
     indexed_ppg(:,2) = ppg_fft_abs(1:fax_ppg_max);
     indexed_ppg(:,1) = fax_ppg(1:fax_ppg_max)';

    % Find PPG signals and corresponding frequencies within the bonds
    esti = (est_ppg_low < indexed_ppg(:,1))&(indexed_ppg(:,1) < est_ppg_high);
    est_ppg_fft = indexed_ppg(esti,:);

    [~,x] = max(est_ppg_fft(:,2));
    est_freq = est_ppg_fft(x,1);
    est_val = est_ppg_fft(x,2);

    
    % Data gathering for scatter plot =====================================
    % ECG HR / PPG HR / Magnitute of PPG signal

    ppg_hr = est_freq*60;
    
    scatter_data(i,1) = ecg_hr(i);
    scatter_data(i,2) = ppg_hr;
    scatter_data(i,3) = est_val;

    % Plot each iteration =================================================

    t=tiledlayout(2,2);
    title(t, ['Window Size: ', num2str(ws)], ...
        ['Window: ', num2str(step), '/', num2str(max_step)])
    
    nexttile
    hold on
    plot(ecg_win)
    plot(inx, val,'xr')
    xlabel('Time (10ms)')
    ylabel('Voltage');
    hold off
    xlim([0 ws])
    ylim([-0.0007 0.0012])
    title('ECG');
    
    nexttile
    plot(ecg_hr)
    title('ECG HR', num2str(ecg_hr(i)))
    xlabel('Window num')
    ylabel('HR');
    xlim([1 max_step])
    ylim([60 200])

    nexttile
    plot(fax_ppg(1:fax_ppg_max), ppg_fft_abs(1:fax_ppg_max),'g-');
    title('PPG HR', num2str(est_freq*60))
    xlabel('Hz')
    ylabel('Magnitude');
    title('PPG FFT');
    xlim([0 10])

    hold on
    plot(est_freq, est_val,'*r');
    xlim([0 10])
    xline(est_ppg_low)
    xline(est_ppg_high)
    hold off

%     nexttile
%     stem(est_ppg_mid,'g-');
%     xlabel('Hz')
%     ylabel('Magnitude');
%     title('EST');
    
%      waitforbuttonpress
    pause(0.1);
end
% END OF FOR LOOP

%% Scatter plot ===========================================================

nexttile
scatter(scatter_data(1:end,1), scatter_data(1:end,3))
title('Heart rate PPG magnitute correlation')
xlabel('Heart rate')
ylabel('Magnitude of PPG');
