%% Load the data
% ECG HR / PPG HR / Magnitute of PPG signal

fs = 100;

clc

% Loading .mat file 
load('scatter_data_sub1.mat');
sub1_data = scatter_data;
ecg_hr_1 = sub1_data(1:end,1);
ppg_hr_1 = sub1_data(1:end,2);
ppg_mag_1 = sub1_data(1:end,3);

load('scatter_data_sub2.mat');
sub2_data = scatter_data;
ecg_hr_2 = sub2_data(1:end,1);
ppg_hr_2 = sub2_data(1:end,2);
ppg_mag_2 = sub2_data(1:end,3);

load('scatter_data_sub3.mat');
sub3_data = scatter_data;
ecg_hr_3 = sub3_data(1:end,1);
ppg_hr_3 = sub3_data(1:end,2);
ppg_mag_3 = sub3_data(1:end,3);

load('scatter_data_sub4.mat');
sub4_data = scatter_data;
ecg_hr_4 = sub4_data(1:end,1);
ppg_hr_4 = sub4_data(1:end,2);
ppg_mag_4 = sub4_data(1:end,3);

% Normalize PPG Data ======================================================

ppg_mag_1_norm = normalize(ppg_mag_1, 'range', [-1 1]);
ppg_mag_2_norm = normalize(ppg_mag_2, 'range', [-1 1]);
ppg_mag_3_norm = normalize(ppg_mag_3, 'range', [-1 1]);
ppg_mag_4_norm = normalize(ppg_mag_4, 'range', [-1 1]);


% Merge the Data ==========================================================

ecg_hr_meta = [ecg_hr_1; ecg_hr_2; ecg_hr_3; ecg_hr_4];
ppg_mag_meta = [ppg_mag_1_norm; ppg_mag_2_norm; ppg_mag_3_norm; ...
    ppg_mag_4_norm];

data_meta(:, 1) = ecg_hr_meta;
data_meta(:, 2) = ppg_mag_meta;

data_meta = sortrows(data_meta);


% Calculate Statistics ====================================================

% y = px + S

% Sub 1
mdl1 = fitlm(ecg_hr_1, ppg_mag_1);
r2o1 = mdl1.Rsquared.Ordinary;
[p1,S1] = polyfit(ecg_hr_1, ppg_mag_1, 1);
% [rho1, pval1] = corr(ecg_hr_1, ppg_mag_1, 'Type', 'Pearson');

% Sub 2
mdl2 = fitlm(ecg_hr_2, ppg_mag_2);
r2o2 = mdl2.Rsquared.Ordinary;
[p2,S2] = polyfit(ecg_hr_2, ppg_mag_2, 1);
% [rho2, pval2] = corr(ecg_hr_2, ppg_mag_2, 'Type', 'Pearson');

% Sub 3
mdl3 = fitlm(ecg_hr_3, ppg_mag_3);
r2o3 = mdl3.Rsquared.Ordinary;
[p3,S3] = polyfit(ecg_hr_3, ppg_mag_3, 1);
% [rho3, pval3] = corr(ecg_hr_1, ppg_mag_1, 'Type', 'Pearson');

% Sub 4
mdl4 = fitlm(ecg_hr_4, ppg_mag_4);
r2o4 = mdl4.Rsquared.Ordinary;
[p4,S4] = polyfit(ecg_hr_4, ppg_mag_4, 1);
% [rho4, pval4] = corr(ecg_hr_4, ppg_mag_4, 'Type', 'Pearson');

% META
mdl_meta = fitlm(data_meta(1:end, 1), data_meta(1:end, 2));
r2o_meta = mdl_meta.Rsquared.Ordinary;
[p_meta, S_meta] = polyfit(data_meta(1:end, 1), data_meta(1:end, 2), 1);

% plotResiduals(mdl_meta)

% Plotting ================================================================

t=tiledlayout(2,4);

nexttile
plot(ecg_hr_1, ppg_mag_1, '.');
h1 = lsline;
h1.Color = 'r';
legend('Data','Fit')
title('Sub 1')
subtitle({['y = ', num2str(p1(1)), 'x + ', num2str(p1(2))], ...
    ['R Squared = ', num2str(r2o1)]})
xlabel('Heart Rate')
ylabel('PPG Magnitute');

nexttile
plot(ecg_hr_2, ppg_mag_2, '.')
h1 = lsline;
h1.Color = 'r';
legend('Data','Fit')
title('Sub 2')
subtitle({['y = ', num2str(p2(1)), 'x + ', num2str(p2(2))], ...
    ['R Squared = ', num2str(r2o2)]})
xlabel('Heart Rate')
ylabel('PPG Magnitute');

nexttile
plot(ecg_hr_3, ppg_mag_3, '.')
h1 = lsline;
h1.Color = 'r';
legend('Data','Fit')
title('Sub 3')
subtitle({['y = ', num2str(p3(1)), 'x + ', num2str(p3(2))], ...
    ['R Squared = ', num2str(r2o3)]})
xlabel('Heart Rate')
ylabel('PPG Magnitute');

nexttile
plot(ecg_hr_4, ppg_mag_4, '.')
h1 = lsline;
h1.Color = 'r';
legend('Data','Fit')
title('Sub 4')
subtitle({['y = ', num2str(p4(1)), 'x + ', num2str(p4(2))], ...
    ['R Squared = ', num2str(r2o4)]})
xlabel('Heart Rate')
ylabel('PPG Magnitute');

nexttile(t,5,[1 4])
plot(data_meta(1:end,1),data_meta(1:end,2),'.')
h1 = lsline;
h1.Color = 'r';
legend('Data','Fit')
title('Meta Data')
subtitle({['y = ', num2str(p_meta(1)), 'x + ', num2str(p_meta(2))], ...
    ['R Squared = ', num2str(r2o_meta)]})
xlabel('Heart Rate')
ylabel('Normalized PPG Magnitute');








