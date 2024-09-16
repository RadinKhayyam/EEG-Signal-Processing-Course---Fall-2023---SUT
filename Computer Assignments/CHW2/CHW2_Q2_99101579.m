%% CHW2 - EEG Signal Processing - Radin Khayyam - 99101579
%% Q2-Part a
clc; clear;

load('Ex2/Ex2.mat')
fs = 250;

plotEEG(X_org);
title('Original signal','Fontsize',14,'Interpreter','latex');

plotEEG(X_noise_2)
title('Noise 2 signal','Fontsize',14,'Interpreter','latex')

plotEEG(X_noise_4)
title('Noise 4 signal','Fontsize',14,'Interpreter','latex')

P_X = sum(sum(X_org.^2));
P_N_2 = sum(sum(X_noise_2.^2));
P_N_4 = sum(sum(X_noise_4.^2));

SNR_10 = -10;
SNR_20 = -20;

sigma_1_N2 = sqrt((P_X/P_N_2)*10^(-1*(SNR_10)/10));
sigma_2_N2 = sqrt((P_X/P_N_2)*10^(-1*(SNR_20)/10));

sigma_1_N4 = sqrt((P_X/P_N_4)*10^(-1*(SNR_10)/10));
sigma_2_N4 = sqrt((P_X/P_N_4)*10^(-1*(SNR_20)/10));

X_SNR_10_N2 = X_org + sigma_1_N2.*(X_noise_2);
X_SNR_20_N2 = X_org + sigma_2_N2.*(X_noise_2);

X_SNR_10_N4 = X_org + sigma_1_N4.*(X_noise_4);
X_SNR_20_N4 = X_org + sigma_2_N4.*(X_noise_4);

plotEEG(X_SNR_10_N2)
title('-10 db noisy signal with Noise 2','Fontsize',14,'Interpreter','latex');


plotEEG(X_SNR_20_N2)
title('-20 db noisy signal with Noise 2','Fontsize',14,'Interpreter','latex');

plotEEG(X_SNR_10_N4)
title('-10 db noisy signal with Noise 4','Fontsize',14,'Interpreter','latex');


plotEEG(X_SNR_20_N4)
title('-20 db noisy signal with Noise 4','Fontsize',14,'Interpreter','latex');

%% Q2-Part b
clc;
% PCA
[U_10_N2,score_10_N2,~] = pca(X_SNR_10_N2.');
[U_20_N2,score_20_N2,~] = pca(X_SNR_20_N2.');
[U_10_N4,score_10_N4,~] = pca(X_SNR_10_N4.');
[U_20_N4,score_20_N4,~] = pca(X_SNR_20_N4.');

sources_10_N2 = score_10_N2.';
sources_20_N2 = score_20_N2.';
sources_10_N4 = score_10_N4.';
sources_20_N4 = score_20_N4.';

plotEEG(sources_10_N2)
title('Sources for SNR = -10dB, Noise 2, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(sources_20_N2)
title('Sources for SNR = -20dB, Noise 2, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(sources_10_N4)
title('Sources for SNR = -10dB, Noise 4, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(sources_20_N4)
title('Sources for SNR = -20dB, Noise 4, PCA','Fontsize',14,'Interpreter','latex');

[F_10_N2,W_10_N2,~]=COM2R(X_SNR_10_N2,32);
ICA_10_N2_components = W_10_N2*X_SNR_10_N2;

[F_20_N2,W_20_N2,~]=COM2R(X_SNR_20_N2,32);
ICA_20_N2_components = W_20_N2*X_SNR_20_N2;

[F_10_N4,W_10_N4,~]=COM2R(X_SNR_10_N4,32);
ICA_10_N4_components = W_10_N4*X_SNR_10_N4;

[F_20_N4,W_20_N4,~]=COM2R(X_SNR_20_N4,32);
ICA_20_N4_components = W_20_N4*X_SNR_20_N4;

plotEEG(ICA_10_N2_components);
title('Sources for SNR = -10dB, Noise 2, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(ICA_20_N2_components);
title('Sources for SNR = -20dB, Noise 2, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(ICA_10_N4_components);
title('Sources for SNR = -10dB, Noise 4, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(ICA_20_N4_components);
title('Sources for SNR = -20dB, Noise 4, ICA','Fontsize',14,'Interpreter','latex');

%% Q2-Part c
clc; 

desired_PCA_10_N2 = 1:8;
desired_PCA_20_N2 = 1:8;
desired_PCA_10_N4 = 1:8;
desired_PCA_20_N4 = 1:8;

desired_ICA_10_N2 = [6,12,15,25];
desired_ICA_20_N2 = [15,19,21];
desired_ICA_10_N4 = [2,3,8,10,27,28];
desired_ICA_20_N4  = [7,8,16,25,28,30];

%% Q2-Part d
clc;

X_denoised_PCA_10_N2 = inv(U_10_N2.') * [sources_10_N2(desired_PCA_10_N2,:) ; zeros(32-length(desired_PCA_10_N2),10240)];
X_denoised_PCA_20_N2 = inv(U_20_N2.') * [sources_20_N2(desired_PCA_20_N2,:) ; zeros(32-length(desired_PCA_20_N2),10240)];
X_denoised_PCA_10_N4 = inv(U_10_N4.') * [sources_10_N4(desired_PCA_10_N4,:) ; zeros(32-length(desired_PCA_10_N4),10240)];
X_denoised_PCA_20_N4 = inv(U_20_N4.') * [sources_20_N4(desired_PCA_20_N4,:) ; zeros(32-length(desired_PCA_20_N4),10240)];

plotEEG(X_denoised_PCA_10_N2);
title('denoised signal for SNR = -10dB, Noise 2, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_PCA_20_N2);
title('denoised signal for SNR = -20dB, Noise 2, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_PCA_10_N4);
title('denoised signal for SNR = -10dB, Noise 4, PCA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_PCA_20_N4);
title('denoised signal for SNR = -20dB, Noise 4, PCA','Fontsize',14,'Interpreter','latex');


X_denoised_ICA_10_N2 =F_10_N2(:,desired_ICA_10_N2)*ICA_10_N2_components(desired_ICA_10_N2,:);
X_denoised_ICA_20_N2 = F_20_N2(:,desired_ICA_20_N2)*ICA_20_N2_components(desired_ICA_20_N2,:);
X_denoised_ICA_10_N4 = F_10_N4(:,desired_ICA_10_N4)*ICA_10_N4_components(desired_ICA_10_N4,:);
X_denoised_ICA_20_N4 = F_20_N4(:,desired_ICA_20_N4)*ICA_20_N4_components(desired_ICA_20_N4,:);

plotEEG(X_denoised_ICA_10_N2);
title('denoised signal for SNR = -10dB, Noise 2, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_ICA_20_N2);
title('denoised signal for SNR = -20dB, Noise 2, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_ICA_10_N4);
title('denoised signal for SNR = -10dB, Noise 4, ICA','Fontsize',14,'Interpreter','latex');
plotEEG(X_denoised_ICA_20_N4);
title('denoised signal for SNR = -20dB, Noise 4, ICA','Fontsize',14,'Interpreter','latex');

%% Q2 - Part e
clc;

fs = 200;
signal_length = length(X_org)/fs;
t = 1/fs:1/fs:signal_length;

figure;

subplot(5,2,1);
plot(t,X_org(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Original signal Ch13, ICA','Fontsize',14,'Interpreter','latex');


subplot(5,2,3);
plot(t,X_SNR_10_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -10dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,4);
plot(t,X_denoised_ICA_10_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -10dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,5);
plot(t,X_SNR_20_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -20dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,6);
plot(t,X_denoised_ICA_20_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -20dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,7);
plot(t,X_SNR_10_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -10dB - Noise 4, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,8);
plot(t,X_denoised_ICA_10_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -10dB - Noise 4, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,9);
plot(t,X_SNR_20_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -20dB - Noise 4 - ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,10);
plot(t,X_denoised_ICA_20_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -20dB - Noise 4 - ICA','Fontsize',14,'Interpreter','latex');

figure;

subplot(5,2,1);
plot(t,X_org(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Original signal Ch24, ICA','Fontsize',14,'Interpreter','latex');


subplot(5,2,3);
plot(t,X_SNR_10_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -10dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,4);
plot(t,X_denoised_ICA_10_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -10dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,5);
plot(t,X_SNR_20_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -20dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,6);
plot(t,X_denoised_ICA_20_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -20dB - Noise 2, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,7);
plot(t,X_SNR_10_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -10dB - Noise 4, ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,8);
plot(t,X_denoised_ICA_10_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -10dB - Noise 4, ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,9);
plot(t,X_SNR_20_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -20dB - Noise 4 - ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,10);
plot(t,X_denoised_ICA_20_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -20dB - Noise 4 - ICA','Fontsize',14,'Interpreter','latex');

figure;

subplot(5,2,1);
plot(t,X_org(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Original signal Ch13, PCA','Fontsize',14,'Interpreter','latex');


subplot(5,2,3);
plot(t,X_SNR_10_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -10dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,4);
plot(t,X_denoised_PCA_10_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -10dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,5);
plot(t,X_SNR_20_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -20dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,6);
plot(t,X_denoised_PCA_20_N2(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -20dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,7);
plot(t,X_SNR_10_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -10dB - Noise 4, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,8);
plot(t,X_denoised_PCA_10_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -10dB - Noise 4, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,9);
plot(t,X_SNR_20_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch13 - SNR = -20dB - Noise 4 - PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,10);
plot(t,X_denoised_PCA_20_N4(13,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch13 - SNR = -20dB - Noise 4 - PCA','Fontsize',14,'Interpreter','latex');

figure;

subplot(5,2,1);
plot(t,X_org(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Original signal Ch24, PCA','Fontsize',14,'Interpreter','latex');


subplot(5,2,3);
plot(t,X_SNR_10_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -10dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,4);
plot(t,X_denoised_PCA_10_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -10dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,5);
plot(t,X_SNR_20_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -20dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,6);
plot(t,X_denoised_PCA_20_N2(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -20dB - Noise 2, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,7);
plot(t,X_SNR_10_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -10dB - Noise 4, PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,8);
plot(t,X_denoised_PCA_10_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -10dB - Noise 4, PCA','Fontsize',14,'Interpreter','latex');

subplot(5,2,9);
plot(t,X_SNR_20_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Noisy signal Ch24 - SNR = -20dB - Noise 4 - PCA','Fontsize',14,'Interpreter','latex');
subplot(5,2,10);
plot(t,X_denoised_PCA_20_N4(24,:));
xlim([0 signal_length]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('Denosied signal Ch24 - SNR = -20dB - Noise 4 - PCA','Fontsize',14,'Interpreter','latex');



%% Q2 - Part f
clc;

RRMSE_ICA_10_N2 = RRMSE(X_org,X_denoised_ICA_10_N2);
RRMSE_ICA_20_N2 = RRMSE(X_org,X_denoised_ICA_20_N2);
RRMSE_ICA_10_N4 = RRMSE(X_org,X_denoised_ICA_10_N4);
RRMSE_ICA_20_N4 = RRMSE(X_org,X_denoised_ICA_20_N4);

RRMSE_PCA_10_N2 = RRMSE(X_org,X_denoised_PCA_10_N2);
RRMSE_PCA_20_N2 = RRMSE(X_org,X_denoised_PCA_20_N2);
RRMSE_PCA_10_N4 = RRMSE(X_org,X_denoised_PCA_10_N4);
RRMSE_PCA_20_N4 = RRMSE(X_org,X_denoised_PCA_20_N4);

disp('RRMSE of ICA for Noise 2 and SNR -10 dB: ');
disp(RRMSE_ICA_10_N2);
disp('RRMSE of ICA for Noise 2 and SNR -20 dB: ');
disp(RRMSE_ICA_20_N2);
disp('RRMSE of ICA for Noise 4 and SNR -10 dB: ');
disp(RRMSE_ICA_10_N4);
disp('RRMSE of ICA for Noise 4 and SNR -20 dB: ');
disp(RRMSE_ICA_20_N4);

disp('RRMSE of PCA for Noise 2 and SNR -10 dB: ');
disp(RRMSE_PCA_10_N2);
disp('RRMSE of PCA for Noise 2 and SNR -20 dB: ');
disp(RRMSE_PCA_20_N2);
disp('RRMSE of PCA for Noise 4 and SNR -10 dB: ');
disp(RRMSE_PCA_10_N4);
disp('RRMSE of PCA for Noise 4 and SNR -20 dB: ');
disp(RRMSE_PCA_20_N4);

%% Functions
function plotEEG(X) 
load('Ex2/Electrodes.mat') ;
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X,offset,feq,ElecName);
end

function result = RRMSE(X_org,X_den)

temp_num = sum(sum((X_org-X_den).^2,2));
temp_den = sum(sum((X_org.^2),2));

result = sqrt(temp_num)/sqrt(temp_den);

end


