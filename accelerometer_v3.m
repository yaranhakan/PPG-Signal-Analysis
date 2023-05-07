% Author: Omer Hakan Yaran
% New Jersey Institute of Technology, 2021
% Advisor: Moshe Kam

%% [1] Obtaining the Data
% ECG, X acc, Y acc, Z acc, X gyro, Y gyro, Z gyro, IR PPG, G PPG, B PPG

clc

% Load '.mat' experimental data file 
load('processed_sensors_data - sub2.mat');

% sampling frequency
fs = 100;


% max/min acc freq
acc_freq_max = 5;
acc_freq_min = 0.5;


% capturing workpace data
workspaceVar = who;
Var = eval(workspaceVar{1});


% categorizing the data
x_acc = Var(:,2);
y_acc = Var(:,3);
z_acc = Var(:,3);

%% Start the for loop here

% Window size
ws = 1000;

% Step size
step_size = 100;

acc_length = length(x_acc);

% Max step
max_step = floor((acc_length-ws)/step_size);

% preallocation
peak_data = zeros(max_step-1,4);
max_pk_x = zeros(max_step-1, 1);
max_pk_y = zeros(max_step-1, 1);
max_pk_z = zeros(max_step-1, 1);


for i = 50:max_step-1

    % Current step, see the variable max_step for the limit
    step = i;

    st = (step-1)*step_size+1;

    x_win = abs(x_acc(st:st+ws));
    y_win = abs(y_acc(st:st+ws));
    z_win = abs(z_acc(st:st+ws));

%     % find max values
%     max_x = max(x_acc);
%     max_y = max(y_acc);
%     max_z = max(z_acc);
    
    
    % make it zero if its less than 50% of the max
%     x_zero = abs(x_acc);
%     y_zero = abs(y_acc);
%     z_zero = abs(z_acc);
% 
%     max_x_zero = max(x_zero);
%     max_y_zero = max(y_zero);
%     max_z_zero = max(z_zero);
% 
%     x_zero(x_zero < max_x_zero/2) = 0;
%     y_zero(y_zero < max_y_zero/2) = 0;
%     z_zero(y_zero < max_y_zero/2) = 0;
    
    % Calculate FFTs
    % [x_acc_fft, fax_x, ~ ] = FFT_ZeroPadded(x_acc, 100, fs, 0);
    % x_acc_fft_abs = abs(x_acc_fft);
    
    [x_win_fft, fax_x_win, ~ ] = FFT_ZeroPadded(x_win, 100, fs, 0);
    x_win_fft = abs(x_win_fft);
    fax_x_win_t = transpose(fax_x_win);
    
    [y_win_fft, fax_y_win, ~ ] = FFT_ZeroPadded(y_win, 100, fs, 0);
    y_win_fft = abs(y_win_fft);
    fax_y_win_t = transpose(fax_y_win);
    
    [z_win_fft, fax_z_win, ~ ] = FFT_ZeroPadded(z_win, 100, fs, 0);
    z_win_fft = abs(z_win_fft);
    fax_z_win_t = transpose(fax_z_win);
    
    % save accelerometer FFT data wrt frequency
    % x_fft_data(:,1) = fax_x;
    % x_fft_data(:,2) = x_acc_fft_abs;
    
    x_win_fft_data = [fax_x_win_t(:), x_win_fft(:)];
    x_win_fft_data = sortrows(x_win_fft_data);
    
    y_win_fft_data = [fax_y_win_t(:), y_win_fft(:)];
    y_win_fft_data = sortrows(y_win_fft_data);
    
    z_win_fft_data = [fax_z_win_t(:), z_win_fft(:)];
    z_win_fft_data = sortrows(z_win_fft_data);
    
    
    % filtering
    
    x_win_fft_data(:,2) = bandpass(x_win_fft_data(:,2), ...
        [acc_freq_min acc_freq_max], fs);
    y_win_fft_data(:,2) = bandpass(abs(y_win_fft_data(:,2)), ...
        [acc_freq_min acc_freq_max], fs);
    z_win_fft_data(:,2) = bandpass(abs(z_win_fft_data(:,2)), ...
        [acc_freq_min acc_freq_max], fs);

    x_win_fft_data(:,2) = abs(x_win_fft_data(:,2));
    y_win_fft_data(:,2) = abs(y_win_fft_data(:,2));
    z_win_fft_data(:,2) = abs(z_win_fft_data(:,2));

    % find peaks

    idx = x_win_fft_data(:,1) > acc_freq_min;

    [pk_x, lk_x] = findpeaks(x_win_fft_data(idx,2),x_win_fft_data(idx,1), ...
        'MinPeakProminence', 0.02);

    idx = y_win_fft_data(:,1) > acc_freq_min;

    [pk_y, lk_y] = findpeaks(y_win_fft_data(idx,2),y_win_fft_data(idx,1), ...
        'MinPeakProminence',0.04);

    idx = z_win_fft_data(:,1) > acc_freq_min;

    [pk_z, lk_z] = findpeaks(z_win_fft_data(idx,2),z_win_fft_data(idx,1), ...
        'MinPeakProminence',0.04);


    % Find fundamental frequencies

    % x
    if isempty(pk_x)
        max_pk_x(:,1) = 0;

    else
        max_pk_x(:,1) = max(pk_x);
    end

    % y
    if isempty(pk_y)
        max_pk_y(:,1) = 0;

    else
        max_pk_y(:,1) = max(pk_y);
    end

    % z
    if isempty(pk_z)
        max_pk_z(:,1) = 0;

    else
        max_pk_z(:,1) = max(pk_z);
    end


    % save peak data

    peak_data(i,1) = max_pk_x(i,1);
    peak_data(i,2) = max_pk_y(i,1);
    peak_data(i,3) = max_pk_z(i,1);

    % find fundamental freq

    if peak_data(i,1) == 0
        ff_x = 'Not periodic';
    else
        ff_x = num2str(peak_data(i,1));
    end

    if peak_data(i,2) == 0
        ff_y = 'Not periodic';
    else
        ff_y = num2str(peak_data(i,2));
    end

    if peak_data(i,3) == 0
        ff_z = 'Not periodic';
    else
        ff_z = num2str(peak_data(i,3));
    end

    %% plotting

    t = tiledlayout(3,2);

    title(t, ['Window Size: ', num2str(ws)], ...
        ['Window: ', num2str(step), '/', num2str(max_step)])
    
    % x
    nexttile
    plot(x_win_fft_data(:,1), x_win_fft_data(:,2),'b');
    xlabel('Hz')
    ylabel('Magnitude');
    title(['X FFT ', ff_x]);
    xlim([0 acc_freq_max])

    hold on
    plot(lk_x, pk_x, '*')
    hold off

    
    nexttile
    plot(x_win, 'b');
    xlabel('Time')
    ylabel('Magnitude');
    title('X Accelerometer');
    
    % y
    nexttile
    plot(y_win_fft_data(:,1), y_win_fft_data(:,2), 'g');
    xlabel('Hz')
    ylabel('Magnitude');
    title(['Y FFT ', ff_y]);
    xlim([0 acc_freq_max])

    hold on
    plot(lk_y, pk_y, '*')
    hold off
    
    nexttile
    plot(y_win, 'g');
    xlabel('Time')
    ylabel('Magnitude');
    title('Y Accelerometer');
    
    % z
    nexttile
    plot(z_win_fft_data(:,1), z_win_fft_data(:,2), 'r');
    xlabel('Hz')
    ylabel('Magnitude');
    title(['Z FFT ', ff_z]);
    xlim([0 acc_freq_max])

    hold on
    plot(lk_z, pk_z, '*')
    hold off
    
    nexttile
    plot(z_win, 'r');
    xlabel('Time')
    ylabel('Magnitude');
    title('Z Accelerometer');

%     waitforbuttonpress
    pause(0.025);
end


