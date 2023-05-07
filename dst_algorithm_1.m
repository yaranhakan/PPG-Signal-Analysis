% Author O. Hakan Yaran
% Replicating Masimo's DST blood oximetry algorithm
% 2021-2022 NJIT

% This program IS NOT applicable for real-world experiments
% Only works on synthetic data because it does not take into account the
% venous blood which gives the strongest output power to the adaptive
% filter

%%  Import data from .txt

clc
clear
close all

fid = fopen('Amb_Light_PPG.txt','rt');
%   fid = fopen('Clean_PPG.txt','rt');
%   open .txt file 'Clean_PPG.txt with rt argument (read text)

ppg_data = textscan(fid, '%f %f', 'HeaderLines', 1);
% get 2 floating number,

red = ppg_data (1:end,1);
% read red from row 1 to end, column 1

infrared = ppg_data (1:end,2);
% read infrared from row 1 to end, column 2

fclose(fid);
% close 'Clean_PPG.txt'

red_mat = cell2mat(red);
infrared_mat = cell2mat(infrared);

%%  Plot the PPG Signals
% figure
% t=tiledlayout(2,2);
% 
% nexttile
% plot(red_mat,'r-');
% xlabel('Index')
% ylabel('Magnitude');
% title('Red (Unfiltered)');
% 
% nexttile
% plot(infrared_mat,'b-');
% xlabel('Index')
% ylabel('Magnitude');
% title('Infrared (Unfiltered)');

%% Filtering and Normalization
% normalization:
% LPF (0.4 Hz) -> DC
% HPF (0.4 Hz) -> AC
% page 38 figure 3.13

% Sampling frequency
fs = 300;

% calculate time wrt frequency
% time = linspace(0,length(red_mat)/fs, length(red_mat));

% HR 42-300
red_ac = bandpass(red_mat, [0.5 10], fs, Steepness = 0.99);
infrared_ac = bandpass(infrared_mat, [0.5 10], fs, Steepness = 0.99);

% Filter dc
red_dc = lowpass(red_mat,0.5,fs, Steepness = 0.99);
infrared_dc = lowpass(infrared_mat,0.4,fs, Steepness = 0.99);

% Normalization
red_n = red_ac./red_dc;
infrared_n = infrared_ac./infrared_dc;



%% DST algorithm set up

% R Value Calculation
% r = 3.2   30%
% r = 0.4   100%
r = rms(red_n)./rms(infrared_n);

% Estimated SPO2 value using conventional RD/IR method (percentage)
spo2_est = 110-25*r;
spo2_est_max = spo2_est + 2;
spo2_est_min = spo2_est - 2;

% Adaptive Filter step size
% (0.1 MATLAB default)
mu = 0.1;

% Create an object for LMS filter
lms = dsp.LMSFilter(16, "StepSize", mu, 'Method', 'Normalized LMS');

%% Recursive Part
% Try RLS filter in the next version

d = (3.2-0.4)/fs;
i = 0;

% Allocate memory for the array
lms_pow = zeros(1,301);
% rls_pow = zeros(1,301);


for rn = 3.2:-d:0.4

    i = i + 1;
    % counter

    rs = rn.*infrared_n - red_n;
    % calculate reference signal for the adaptive filter


    [lms_signal, err, wts] = lms(rs, infrared_n);
    % lms filtering
    

    lms_pow(1, i) = sum(lms_signal.^2);
    % calculate power of the lms filtered signal
end

rn = 3.2:-d:0.4;
%% Find the SpO2 percentage

% find index of maximum value in lms_pow
[max_val, max_idx] = max(lms_pow);

% find value at same index in rn
value_at_max_idx = rn(max_idx);

% Create an array of percentages from 30% to 100% with the same length as lms_pow
percentages = linspace(30, 100, fs+1);

spo2_dst = percentages(max_idx);

fprintf('SpO2 Found using the conventional method: %.2f\nSpo2 Found using DST algorithm by MASIMO: %.2f\n', spo2_est, spo2_dst);

%%

% scale = 31.25*(0.4:d/10:3.2);
% scale(1) = [];

figure
t=tiledlayout(2,2);

nexttile
plot(red_mat,'b-');
xlabel('Index')
ylabel('Magnitude');
title('Signal with noise');

nexttile
plot(red_n,'r-');
xlabel('Index')
ylabel('Magnitude');
title('Desired Signal');


nexttile
plot(percentages, lms_pow, 'r-');
xlabel('Percentage (%)')
ylabel('Magnitude');
title('Filtered Output');


nexttile
hold on
plot(spo2_est, 'b*')
plot(spo2_dst, 'r*')
xlabel('Percentage (%)')
ylabel('Magnitude');
ylim([95 100])
title('SpO2 obtained from conventional method (blue) and DST (red)')
hold off




