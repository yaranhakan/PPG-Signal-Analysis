% Author: Omer Hakan Yaran
% New Jersey Institute of Technology, 2021
% Advisor: Moshe Kam

%%   PART I: import data from .txt

clc
fid = fopen('Amb_Light_PPG.txt','rt');
% open .txt file 'Clean_PPG.txt with rt argument (read text)

ppg_data = textscan(fid, '%f %f', 'HeaderLines', 1);
% get 2 floating number

red = ppg_data (1:end,1);
% read red from row 1 to end, column 1

infrared = ppg_data (1:end,2);
% read infrared from row 1 to end, column 2

fclose(fid);
% close 'Clean_PPG.txt'

red_mat = cell2mat(red);
infrared_mat = cell2mat(infrared);



%% Part II: Filtering and Normalization
% normalization:
% LPF (0.4 Hz) -> DC
% HPF (0.4 Hz) -> AC
% page 38 figure 3.13

fs = 300;

time = linspace(0,length(red_mat)/fs, length(red_mat));
% calculate time wrt frequency

red_ac = bandpass(red_mat, [0.5 10], fs, Steepness = 0.99);
infrared_ac = bandpass(infrared_mat, [0.5 10], fs, Steepness = 0.99);
% HR 42-300

red_dc = lowpass(red_mat,0.5,fs, Steepness = 0.99);
infrared_dc = lowpass(infrared_mat,0.4,fs, Steepness = 0.99);
% filter dc

red_n = red_ac./red_dc;
infrared_n = infrared_ac./infrared_dc;
% normalization


%%   PART III: Fast Fourier Transform

padding_val = 100;
padding_p = nextpow2(padding_val);
padding_pts = 2^padding_p;
% padding_pts is the closest power of 2 higher than padding_val

lpad = padding_pts*length(red_n);
% padding value (10k)

fax_hz = 0:fs/lpad:fs/2;
% frequency axis

hr_scale = fax_hz*60;

red_n_fft_abs2 = abs(FFT_ZeroPadded(red_n, padding_pts, 300, 0));
infrared_n_fft_abs2 = abs(FFT_ZeroPadded(infrared_n, padding_pts, 300, 0));
% calculate FFTs

%%   PART IV: Calculate HR
% HR(BPM) = (index number - 1) * (fs) / (sample number) * 60

[red_hr_f, index] = max(red_n_fft_abs2);
% get max and index

hr_bpm = hr_scale(index);

% calculate HR BPM

disp(['HR: ' num2str(round(hr_bpm,0)),' bpm']);
% display HR BPM

%%   PART V: Calculate SpO2
% red is used to measure desaturated hemoglobin
% infrared is used to measure oxygenated blood
% R = (ACred - DCred) / (ACinfrared - DCinfrared)

r = rms(red_n)./rms(infrared_n);
% R value calculation
% peak-to-peak value division

spo2 = 110-25*r;
% calculate SpO2

disp(['SpO2: ', num2str(round(spo2,1)), ' %']);
% display SpO2


%% PART VI: Plot the PPG Signals

figure
t=tiledlayout(3,2);

    title(t, ['Heart Rate: ', num2str(round(hr_bpm,0)), 'bpm'], ...
        ['Blood Oxygen Level: ', num2str(round(spo2,1)), '%'])
% Plot noisy signal
nexttile
plot(red_mat,'r-');
xlabel('Index')
ylabel('Magnitude');
title('Red (Unfiltered)');

nexttile
plot(infrared_mat,'b-');
xlabel('Index')
ylabel('Magnitude');
title('Infrared (Unfiltered)');

% Plot filtered signal
nexttile
plot(red_ac,'r');
xlabel('Index')
ylabel('Magnitude');
title('Red (Band-pass Filtered)');

nexttile
plot(infrared_ac,'b');
xlabel('Index')
ylabel('Magnitude');
title('Infrared ( Band-pass Filtered)');

% Plot FFT
nexttile
plot(fax_hz, red_n_fft_abs2,'r-');
xlabel('Frequency (Hz)')
ylabel('Magnitude');
title('Red FFT');
xlim([0,5]);

nexttile
plot(fax_hz, infrared_n_fft_abs2,'b-');
xlabel('Frequency (Hz)')
ylabel('Magnitude');
title('Infrared FFT');
xlim([0,5]);
