% Radin Khayyam - 99101579
% EEG Signal Processing - CHW_1

%% Q1 - Part A
clc; close all; clear;

fs = 1000;
t = 0:1/fs:2;
f0 = 100;
beta = 100;
x = chirp(t,100,1,f0+beta,'quadratic');
pspectrum(x,1e3,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',99,'Leakage',0.85);
%% Q1 - Part B
clc; close all;
L = 128;
f = (-L/2:L/2-1)*(fs/L);

% Rectangular window
rect_window = rectwin(L);
fft_rect_window = fftshift(fft(rect_window))/L;

% Triangular window
triang_window = triang(L);
fft_triang_window = fftshift(fft(triang_window))/L;

% Gaussian window
gauss_window = gausswin(L, 2.5);
fft_gauss_window = fftshift(fft(gauss_window))/L;

% Hamming window
hamming_window = hamming(L);
fft_hamming_window = fftshift(fft(hamming_window))/L;

% Plot time and frequency domain representations
figure;
subplot(2, 2, 1);
plot(rect_window);
xlim([0 128]);
title('Rectangular Window', 'Interpreter','latex');
subplot(2, 2, 2);
plot(triang_window);
xlim([0 128]);
title('Triangular Window', 'Interpreter','latex');
subplot(2, 2, 3);
plot(gauss_window);
xlim([0 128]);
title('Gaussian Window', 'Interpreter','latex');
subplot(2, 2, 4);
plot(hamming_window);
title('Hamming Window', 'Interpreter','latex');
xlim([0 128]);

figure;
subplot(2, 2, 1);
plot(f,abs(fft_rect_window));
xlim([-100 100]);
title('Frequency Domain - Rectangular Window', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
subplot(2, 2, 2);
plot(f,abs(fft_triang_window));
xlim([-100 100]);
xlabel('frequency(Hz)', 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
title('Frequency Domain - Triangular Window', 'Interpreter','latex');
subplot(2, 2, 3);
plot(f,abs(fft_gauss_window));
xlim([-100 100]);
xlabel('frequency(Hz)', 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
title('Frequency Domain - Gaussian Window', 'Interpreter','latex');
subplot(2, 2, 4);
plot(f,abs(fft_hamming_window));
xlim([-100 100]);
xlabel('frequency(Hz)', 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
title('Frequency Domain - Hamming Window', 'Interpreter','latex');

%% Q1 - Part C
clc; close all;
Noverlap = 0;
nfft = L;

figure;
colormap default

subplot(2, 2, 1);
spectrogram(x, rect_window, Noverlap, nfft, fs, 'yaxis');
title('STFT Using Rectangular Window','FontSize', 14, 'Interpreter','latex');

subplot(2, 2, 2);
spectrogram(x, triang_window, Noverlap, nfft, fs, 'yaxis');
title('STFT Using Triangular Window','FontSize', 14, 'Interpreter','latex');

subplot(2, 2, 3);
spectrogram(x, gauss_window, Noverlap, nfft, fs, 'yaxis');
title('STFT Using Gaussian Window','FontSize', 14, 'Interpreter','latex');

subplot(2, 2, 4);
spectrogram(x, hamming_window, Noverlap, nfft, fs, 'yaxis');
title('STFT Using Hamming Window','FontSize', 14, 'Interpreter','latex');



%% Q1 - Part D
clc; close all;
N_overlap = [0, 64, 127];
L = 128;
nfft = 128;
counter = 0;
for Noverlap = N_overlap            
    hamming_window = hamming(L);
    counter = counter +1;
    subplot(2,2,counter);
    spectrogram(x, hamming_window, Noverlap, nfft, fs, 'yaxis');
    title(['N overlap = ' num2str(Noverlap) ],'FontSize', 14,'Interpreter','latex');
end

%% Q1 - Part E
clc; close all;
window_lengths = [32, 128, 512];
counter = 0;
for L = window_lengths
    nfft = L;               
    Noverlap = L - 1; 
    hamming_window = hamming(L);
    counter = counter +1;
    subplot(2,2,counter);
    spectrogram(x, hamming_window, Noverlap, nfft, fs, 'yaxis');
    title(['L = ' num2str(L) ],'FontSize', 14,'Interpreter','latex');
end
%% Q1 - Part F
clc; close all;
L = 128;
Noverlap = L/2;
hamming_window = hamming(L);
nfft_coef = [1, 2, 4];
counter = 0;
for coef = nfft_coef
    counter = counter +1;
    nfft = coef * L;
    subplot(2,2,counter);
    spectrogram(x, hamming_window, Noverlap, nfft, fs, 'yaxis');
    title(['nfft = ' num2str(nfft) ],'FontSize', 14,'Interpreter','latex');
end

%% Q1 - Part G
clc; close all;
nfft = 128;
window_size = 128;     % Size of the analysis window
overlap = 128;        % Overlap between windows

num_windows = floor(length(x)/window_size);
result = zeros(num_windows,nfft);
for i = 1:num_windows
    selected_part = x(window_size*(i-1)+1:window_size*i);
    fft_selected_part = fft(selected_part,nfft);
    result(i,:) = fft_selected_part;
end
result_2 = spectrogram(x,rectwin(window_size),0,nfft,fs);


%% Q2 - Loading the data
clc;clear; close all;
data = load('NewEEGSignal.mat');
signal = data.NewEEGSignal;
fs = 256;
L = 512;
t= 0:1/fs:(L-1)/fs;


%% Q2 - Part A
clc; close all;
f = (-L/2:L/2-1)*(fs/L);
fft_signal = abs(fftshift(fft(signal)))/L;

subplot(4,2,3);
plot(t,signal,'Linewidth',1);
title('Signal in time domain','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('Time(s)', 'Interpreter','latex');

subplot(4,2,5)
plot(f,fft_signal,'Linewidth',1);
xlim([0 64])
title('DFT of Signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

subplot(4,2,[4,6]);
spectrogram(signal,hamming(128),64,128,fs,'yaxis');
title('Spectrogram of Signal','FontSize', 14, 'Interpreter','latex');

%% Q2 - Part B
clc; close all;

filtered_signal = lowpass(signal,64,fs);
downsampled_signal = downsample(filtered_signal,2);

new_L=length(downsampled_signal);
new_fs = 128;
new_t = 0:1/new_fs:(new_L-1)/new_fs;
new_f = (-new_L/2:new_L/2-1)*(new_fs/new_L);
fft_downsampled_signal = abs(fftshift(fft(downsampled_signal)))/new_L;

subplot(4,2,1);
plot(new_t,downsampled_signal,'Linewidth',1);
title('down sampled signal in time domain','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('Time(s)', 'Interpreter','latex');

subplot(4,2,3);
plot(t,signal,'Linewidth',1);
title('Original signal in time domain','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('Time(s)', 'Interpreter','latex');

subplot(4,2,5)
plot(new_f,fft_downsampled_signal,'Linewidth',1);
xlim([0 64])
title('DFT of down sampled signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

subplot(4,2,7)
plot(f,fft_signal,'Linewidth',1);
xlim([0 64])
title('DFT of original signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

subplot(4,2,[2,4]);
spectrogram(downsampled_signal,hamming(64),32,64,new_fs,'yaxis');
title('Spectrogram of down sampled signal','FontSize', 14, 'Interpreter','latex');

subplot(4,2,[6,8]);
spectrogram(signal,hamming(128),64,128,fs,'yaxis');
title('Spectrogram of original signal','FontSize', 14, 'Interpreter','latex');

%% Q2 - Part C
clc; close all;
fs = 128;

signal_n2 = downsampled_signal(1:new_L/2);
L = length(signal_n2);
fft_signal_n2 = abs(fftshift(fft(signal_n2)))/L;
f = (-L/2:L/2-1)*(fs/L);

subplot(3,1,1)
plot(f,fft_signal_n2,'Linewidth',1);
xlim([0 64])
title('DFT of first N/2 samples','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

signal_n4 = downsampled_signal(1:new_L/4);
L = length(signal_n4);
fft_signal_n4 = abs(fftshift(fft(signal_n4)))/L;
f = (-L/2:L/2-1)*(fs/L);

subplot(3,1,2)
plot(f,fft_signal_n4,'Linewidth',1);
xlim([0 64])
title('DFT of first N/4 samples','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

signal_n8 = downsampled_signal(1:new_L/8);
L = length(signal_n8);
fft_signal_n8 = abs(fftshift(fft(signal_n8)))/L;
f = (-L/2:L/2-1)*(fs/L);

subplot(3,1,3)
plot(f,fft_signal_n8,'Linewidth',1);
xlim([0 64])
title('DFT of first N/8 samples','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

% I don't understand what should i do exactly so, i will calculate N/2, N/4
% and N/8 point DFT of the down sampled signal also.

L = length(downsampled_signal);


fft_downsampled_signal_n2 = abs(fftshift(fft(downsampled_signal, L/2)))/(L/2);
f = (-L/4:L/4-1)*(fs/(L/2));
figure;
subplot(3,1,1)
plot(f,fft_downsampled_signal_n2,'Linewidth',1);
xlim([0 64])
title('N/2 point DFT of down sampled signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

fft_downsampled_signal_n4 = abs(fftshift(fft(downsampled_signal, L/4)))/(L/4);
f = (-L/8:L/8-1)*(fs/(L/4));
subplot(3,1,2)
plot(f,fft_downsampled_signal_n4,'Linewidth',1);
xlim([0 64])
title('N/4 point DFT of down sampled signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

fft_downsampled_signal_n8 = abs(fftshift(fft(downsampled_signal, L/8)))/(L/8);
f = (-L/16:L/16-1)*(fs/(L/8));
subplot(3,1,3)
plot(f,fft_downsampled_signal_n8,'Linewidth',1);
xlim([0 64])
title('N/8 point DFT of down sampled signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

%% Q2 - Part D
clc; close all;
L = length(downsampled_signal);
fs = 128;
signal_n2_zrpad = [signal_n2, zeros(1,L/2)];
signal_n4_zrpad = [signal_n4, zeros(1,3*L/4)];
signal_n8_zrpad = [signal_n8, zeros(1,7*L/8)];
f = (-L/2:L/2-1)*(fs/L);

fft_signal_n2_zrpad = abs(fftshift(fft(signal_n2_zrpad, L)))/L;
fft_signal_n4_zrpad = abs(fftshift(fft(signal_n4_zrpad, L)))/L;
fft_signal_n8_zrpad = abs(fftshift(fft(signal_n8_zrpad, L)))/L;

figure;
subplot(3,1,1)
plot(f,fft_signal_n2_zrpad,'Linewidth',1);
xlim([0 64])
title('DFT of zero paded N/2 point signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

subplot(3,1,2)
plot(f,fft_signal_n4_zrpad,'Linewidth',1);
xlim([0 64])
title('DFT of zero paded N/4 point signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

subplot(3,1,3)
plot(f,fft_signal_n8_zrpad,'Linewidth',1);
xlim([0 64])
title('DFT of zero paded N/8 point signal','FontSize', 14, 'Interpreter','latex');
ylabel('magnitude', 'Interpreter','latex');
xlabel('frequency(Hz)', 'Interpreter','latex');

%% Q3 - Loading the Data
clc; close all; clear;
data = load('NewEEGSignal.mat');
signal = data.NewEEGSignal;
fs = 256;
L = 512;
t= 0:1/fs:(L-1)/fs;
f = (-L/2:L/2-1)*(fs/L);

%% Q3 - All Parts
clc; close all;

% Using correlation
fs = 256;
signal_acorr = xcorr(signal);
L = length(signal_acorr);  
f = (0:L/2)*(fs/L);
m = fft(signal_acorr);
m = m(1:round(L/2));
fft_signal_acorr = abs(m)/L;

figure;
subplot(5,1,1);
plot(f,fft_signal_acorr,'Linewidth',1);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');
title('PSD using auto correlation','FontSize',14,'Interpreter','latex');

% Using Circular Correlation
signal_acorr = ifft(fft(signal).*conj(fft(signal)));
L = length(signal_acorr);  
f = (0:L/2-1)*(fs/L);
m = fft(signal_acorr);
m = m(1:round(L/2));
fft_signal_acorr = abs(m)/L;

subplot(5,1,2);
plot(f,fft_signal_acorr,'Linewidth',1);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');
title('PSD using auto correlation','FontSize',14, 'Interpreter','latex');

% Using Periodogram
subplot(5,1,3);
[pxx,f] = periodogram(signal,[],[],fs);
plot(f,pxx,'Linewidth',1);
title('PSD using periodogram','FontSize',14, 'Interpreter','latex');
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');

% Using Pwelch
[pxx,f] = pwelch(signal,[],[],[],fs);
subplot(5,1,4);
plot(f,pxx,'Linewidth',1);
title('PSD using pwelch','FontSize',14, 'Interpreter','latex');
d

% Using Pwelch with rect window
L = length(signal);
[pxx,f] = pwelch(signal,rectwin(L),[],[],fs);
subplot(5,1,5);
plot(f,pxx,'Linewidth',1);
title('PSD using pwelch with rect window','FontSize',14, 'Interpreter','latex');
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');
