%% EEG Signal Processing - Radin Khayyam - 99101579
%% CHW4
%% Q1
clc;clear; close all;

ERP_EEG = load('ERP_EEG.mat');
signal = ERP_EEG.ERP_EEG;
signal = transpose(signal);

%% Q1 - Part a

fs = 240;
t = 0:1/fs:1-1/fs;
counter = 1;
figure;
for N=100:100:2500
    signal_mean = mean(signal(1:N,:));
    subplot(5,5,counter);
    plot(t,signal_mean);
    xlabel('time')
    grid on;
    hold on;
    title(['mean of ',num2str(N),' samples'],'Interpreter','latex');
    counter = counter + 1;
end

%% Q1 - Part b

clc; close all;

max_values = zeros(1,2550);
mean_values = zeros(2550,240);
for N=1:2550
    mean_values(N,:) = mean(signal(1:N,:));
    signal_max = max(mean_values(N,:));
    max_values(N) = signal_max;
end
plot(max_values,'Linewidth',1);
xlim([1 2550])
title('Max value of first N samples mean','Interpreter','latex');
xlabel('N','Interpreter','latex');
ylabel('max value','Interpreter','latex');

%% Q1 - Part c

clc; close all;

rms_values = zeros(1,2549);
for i = 1:2549
    rms_values(i) = rms(mean_values(i+1)-mean_values(i));
end
plot(rms_values,'Linewidth',1);
xlabel('i');
grid on;
title('rms value between i and i+1 samples','Interpreter','latex');

%% Q1 - Part d

clc; close all;

N0 = 900;
mean_1 = mean(signal(1:N0,:));
mean_2 = mean(signal(1:2550,:));
mean_3 = mean(signal(1:(N0/3),:));
mean_4 = mean(signal(randi(2550,1,N0),:));
mean_5 = mean(signal(randi(2550,1,N0/3),:));

figure();
subplot(2,3,1);
plot(mean_1,'Linewidth',1)
xlim([0 240])
title('N = 900','Interpreter','latex');

subplot(2,3,2);
plot(mean_2,'Linewidth',1)
xlim([0 240])
title('N = 2550','Interpreter','latex');

subplot(2,3,3);
plot(mean_3,'Linewidth',1)
xlim([0 240])
title('N = 300','Interpreter','latex');

subplot(2,3,4);
plot(mean_4,'Linewidth',1)
xlim([0 240])
title('N = 900 random samples','Interpreter','latex');

subplot(2,3,5);
plot(mean_5,'Linewidth',1)
xlim([0 240])
title('N = 300 random samples','Interpreter','latex');

%% Q2
clc; clear; close all;
load("SSVEP_EEG.mat");
fs = 250;

%% Q2 - Part a - 1
clc;
signal = zeros(size(SSVEP_Signal,1),size(SSVEP_Signal,2));
for i = 1:6
    signal(i,:) = bandpass(SSVEP_Signal(i,:),[1 40],fs);
end
%% Q2 - Part a - 2
clc;
Evoked_signals = cell(15,1);
for i = 1:15
    Evoked_signals{i} = signal(:,Event_samples(i):(Event_samples(i)+fs*5-1));
end
%% Q2 - Part a - 3
clc;
figure;
for i = 1:15
    subplot(5,3,i);

    [pxx1,f] = pwelch(Evoked_signals{i}(1,:),[],[],[],fs);
    [pxx2,f] = pwelch(Evoked_signals{i}(2,:),[],[],[],fs);
    [pxx3,f] = pwelch(Evoked_signals{i}(3,:),[],[],[],fs);
    [pxx4,f] = pwelch(Evoked_signals{i}(4,:),[],[],[],fs);
    [pxx5,f] = pwelch(Evoked_signals{i}(5,:),[],[],[],fs);
    [pxx6,f] = pwelch(Evoked_signals{i}(6,:),[],[],[],fs);
    plot(f,pxx1,'Linewidth',1);
    hold on;
    plot(f,pxx2,'Linewidth',1);
    hold on;
    plot(f,pxx3,'Linewidth',1);
    hold on;
    plot(f,pxx4,'Linewidth',1);
    hold on;
    plot(f,pxx5,'Linewidth',1);
    hold on;
    plot(f,pxx6,'Linewidth',1);
    title(['Trial ',num2str(i)],'Interpreter','latex','FontSize',12);
    legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6');
    xlabel('frequency (Hz)','Interpreter','latex','FontSize',8);
end

%% Q2 - Part a - 4
clc; close all; 
figure;
for i = 1:6
    subplot(6,1,i);
    [pxx,f] = pwelch(Evoked_signals{5}(i,:),[],[],[],fs);
    plot(f,pxx,'Linewidth',1);
    title(['Trial 5 - Chanel ',num2str(i)],'Interpreter','latex','FontSize',8);
    xlabel('frequency (Hz)','Interpreter','latex')
end
%% Q2 - Part a - 5
clc; close all; 
figure;
for i = 1:15
    subplot(5,3,i);
    [pxx1,f] = pwelch(Evoked_signals{i}(1,:),[],[],[],fs);
    [pxx2,f] = pwelch(Evoked_signals{i}(2,:),[],[],[],fs);
    [pxx3,f] = pwelch(Evoked_signals{i}(3,:),[],[],[],fs);
    [pxx4,f] = pwelch(Evoked_signals{i}(4,:),[],[],[],fs);
    [pxx5,f] = pwelch(Evoked_signals{i}(5,:),[],[],[],fs);
    [pxx6,f] = pwelch(Evoked_signals{i}(6,:),[],[],[],fs);
    mean_pxx = (pxx1 + pxx2 + pxx3 + pxx4 + pxx5 + pxx6)/6;
    plot(f,mean_pxx,'LineWidth',1);
    title(['Trial ',num2str(i),' - Mean of channels'],'Interpreter','latex','FontSize',8);
end

%% Q2 - Part b - 1
clc; close all; clear;
load("SSVEP_EEG.mat");
fs = 250;
Evoked_signals = cell(15,1);
for i = 1:15
    Evoked_signals{i} = SSVEP_Signal(:,Event_samples(i):(Event_samples(i)+fs*5-1));
end
%% Q2 - Part b - 2
clc;
t = 0 : 1/fs : 5- 1/fs;
f1 = 6.5;
f2 = 7.35;
f3 = 8.3;
f4 = 9.6;
f5 = 11.61;
Y_f1 = [sin(2*pi*f1*t);cos(2*pi*f1*t);sin(2*pi*2*f1*t);cos(2*pi*2*f1*t);sin(2*pi*3*f1*t);cos(2*pi*3*f1*t);sin(2*pi*4*f1*t);cos(2*pi*4*f1*t);sin(2*pi*5*f1*t);cos(2*pi*5*f1*t);sin(2*pi*6*f1*t);cos(2*pi*6*f1*t)];
Y_f2 = [sin(2*pi*f2*t);cos(2*pi*f2*t);sin(2*pi*2*f2*t);cos(2*pi*2*f2*t);sin(2*pi*3*f2*t);cos(2*pi*3*f2*t);sin(2*pi*4*f2*t);cos(2*pi*4*f2*t);sin(2*pi*5*f2*t);cos(2*pi*5*f2*t)];
Y_f3 = [sin(2*pi*f3*t);cos(2*pi*f3*t);sin(2*pi*2*f3*t);cos(2*pi*2*f3*t);sin(2*pi*3*f3*t);cos(2*pi*3*f3*t);sin(2*pi*4*f3*t);cos(2*pi*4*f3*t)];
Y_f4 = [sin(2*pi*f4*t);cos(2*pi*f4*t);sin(2*pi*2*f4*t);cos(2*pi*2*f4*t);sin(2*pi*3*f4*t);cos(2*pi*3*f4*t);sin(2*pi*4*f4*t);cos(2*pi*4*f4*t)];
Y_f5 = [sin(2*pi*f5*t);cos(2*pi*f5*t);sin(2*pi*2*f5*t);cos(2*pi*2*f5*t);sin(2*pi*3*f5*t);cos(2*pi*3*f5*t)];

for i=1:15
    [~,~,r1] = canoncorr(Evoked_signals{i}',Y_f1');
    [~,~,r2] = canoncorr(Evoked_signals{i}',Y_f2');
    [~,~,r3] = canoncorr(Evoked_signals{i}',Y_f3');
    [~,~,r4] = canoncorr(Evoked_signals{i}',Y_f4');
    [~,~,r5] = canoncorr(Evoked_signals{i}',Y_f5');
    r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
    switch find(r==max(r))
        case 1
            disp(['Stimulation frequency of trial number ',num2str(i),' is 6.5Hz']);
        case 2
            disp(['Stimulation frequency of trial number ',num2str(i),' is 7.35Hz']);
        case 3
            disp(['Stimulation frequency of trial number ',num2str(i),' is 8.3Hz']);
        case 4
            disp(['Stimulation frequency of trial number ',num2str(i),' is 9.6Hz']);
        case 5
            disp(['Stimulation frequency of trial number ',num2str(i),' is 11.61Hz']);

    end
end

%% Q2 - Part b - 3
clc;
five_channels_signal = cell(15,1);
four_channels_signal = cell(15,1);
three_channels_signal = cell(15,1);
two_channels_signal = cell(15,1);
one_channel_signal = cell(15,1);

mean_accuracy = zeros(1,5);
for k=1:10
    rp_channels = randperm(6); % random permutation of channels
    for i=1:15
        five_channels_signal{i} = Evoked_signals{i}(rp_channels(1:5),:);
        four_channels_signal{i} = Evoked_signals{i}(rp_channels(1:4),:);
        three_channels_signal{i} = Evoked_signals{i}(rp_channels(1:3),:);
        two_channels_signal{i} = Evoked_signals{i}(rp_channels(1:2),:);
        one_channel_signal{i} = Evoked_signals{i}(rp_channels(1),:);
    end
    
    Ground_truth = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];
    
    selected_frequencies = cell(5,1);
    for i=1:15
        [~,~,r1] = canoncorr(one_channel_signal{i}',Y_f1');
        [~,~,r2] = canoncorr(one_channel_signal{i}',Y_f2');
        [~,~,r3] = canoncorr(one_channel_signal{i}',Y_f3');
        [~,~,r4] = canoncorr(one_channel_signal{i}',Y_f4');
        [~,~,r5] = canoncorr(one_channel_signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{1}(i) = find(r == max(r));
       
        [~,~,r1] = canoncorr(two_channels_signal{i}',Y_f1');
        [~,~,r2] = canoncorr(two_channels_signal{i}',Y_f2');
        [~,~,r3] = canoncorr(two_channels_signal{i}',Y_f3');
        [~,~,r4] = canoncorr(two_channels_signal{i}',Y_f4');
        [~,~,r5] = canoncorr(two_channels_signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{2}(i) = find(r == max(r));
    
        [~,~,r1] = canoncorr(three_channels_signal{i}',Y_f1');
        [~,~,r2] = canoncorr(three_channels_signal{i}',Y_f2');
        [~,~,r3] = canoncorr(three_channels_signal{i}',Y_f3');
        [~,~,r4] = canoncorr(three_channels_signal{i}',Y_f4');
        [~,~,r5] = canoncorr(three_channels_signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{3}(i) = find(r == max(r));
    
        [~,~,r1] = canoncorr(four_channels_signal{i}',Y_f1');
        [~,~,r2] = canoncorr(four_channels_signal{i}',Y_f2');
        [~,~,r3] = canoncorr(four_channels_signal{i}',Y_f3');
        [~,~,r4] = canoncorr(four_channels_signal{i}',Y_f4');
        [~,~,r5] = canoncorr(four_channels_signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{4}(i) = find(r == max(r));
    
        [~,~,r1] = canoncorr(five_channels_signal{i}',Y_f1');
        [~,~,r2] = canoncorr(five_channels_signal{i}',Y_f2');
        [~,~,r3] = canoncorr(five_channels_signal{i}',Y_f3');
        [~,~,r4] = canoncorr(five_channels_signal{i}',Y_f4');
        [~,~,r5] = canoncorr(five_channels_signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{5}(i) = find(r == max(r));
        
    end
    
    accuracies = zeros(1,5);
    for i = 1:5
        accuracies(i) = sum(selected_frequencies{i} == Ground_truth)/15*100;
    end
    
    mean_accuracy = mean_accuracy + accuracies;
end

mean_accuracy = mean_accuracy/10;
disp(['mean accuracy for one channel signal = ',num2str(mean_accuracy(1))]);
disp(['mean accuracy for two channels signal = ',num2str(mean_accuracy(2))]);
disp(['mean accuracy for three channels signal = ',num2str(mean_accuracy(3))]);
disp(['mean accuracy for four channels signal = ',num2str(mean_accuracy(4))]);
disp(['mean accuracy for five channels signal = ',num2str(mean_accuracy(5))]);

%% Q2 - Part b - 4
clc;
signal = cell(15,1);
selected_frequencies = cell(9,1);
for k=0.5:0.5:4.5
    t = 0:1/fs:k-1/fs;
    Y_f1 = [sin(2*pi*f1*t);cos(2*pi*f1*t);sin(2*pi*2*f1*t);cos(2*pi*2*f1*t);sin(2*pi*3*f1*t);cos(2*pi*3*f1*t);sin(2*pi*4*f1*t);cos(2*pi*4*f1*t);sin(2*pi*5*f1*t);cos(2*pi*5*f1*t);sin(2*pi*6*f1*t);cos(2*pi*6*f1*t)];
    Y_f2 = [sin(2*pi*f2*t);cos(2*pi*f2*t);sin(2*pi*2*f2*t);cos(2*pi*2*f2*t);sin(2*pi*3*f2*t);cos(2*pi*3*f2*t);sin(2*pi*4*f2*t);cos(2*pi*4*f2*t);sin(2*pi*5*f2*t);cos(2*pi*5*f2*t)];
    Y_f3 = [sin(2*pi*f3*t);cos(2*pi*f3*t);sin(2*pi*2*f3*t);cos(2*pi*2*f3*t);sin(2*pi*3*f3*t);cos(2*pi*3*f3*t);sin(2*pi*4*f3*t);cos(2*pi*4*f3*t)];
    Y_f4 = [sin(2*pi*f4*t);cos(2*pi*f4*t);sin(2*pi*2*f4*t);cos(2*pi*2*f4*t);sin(2*pi*3*f4*t);cos(2*pi*3*f4*t);sin(2*pi*4*f4*t);cos(2*pi*4*f4*t)];
    Y_f5 = [sin(2*pi*f5*t);cos(2*pi*f5*t);sin(2*pi*2*f5*t);cos(2*pi*2*f5*t);sin(2*pi*3*f5*t);cos(2*pi*3*f5*t)];

    for i=1:15
        signal{i}=Evoked_signals{i}(:,1:k*fs);
        [~,~,r1] = canoncorr(signal{i}',Y_f1');
        [~,~,r2] = canoncorr(signal{i}',Y_f2');
        [~,~,r3] = canoncorr(signal{i}',Y_f3');
        [~,~,r4] = canoncorr(signal{i}',Y_f4');
        [~,~,r5] = canoncorr(signal{i}',Y_f5');
        r = [r1(1),r2(1),r3(1),r4(1),r5(1)];
        selected_frequencies{k/0.5}(i) = find(r == max(r));
    end
end

Ground_truth = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];
accuracies = zeros(1,9);
for i=1:9
    accuracies(i)=sum(selected_frequencies{i}==Ground_truth)/15*100;
end

disp(['accuracy for 0.5 seconds signal = ',num2str(accuracies(1))]);
disp(['accuracy for 1 seconds = ',num2str(accuracies(2))]);
disp(['accuracy for 1.5 seconds signal = ',num2str(accuracies(3))]);
disp(['accuracy for 2 seconds signal = ',num2str(accuracies(4))]);
disp(['accuracy for 2.5 seconds signal = ',num2str(accuracies(5))]);
disp(['accuracy for 3 seconds signal = ',num2str(accuracies(6))]);
disp(['accuracy for 3.5 seconds signal = ',num2str(accuracies(7))]);
disp(['accuracy for 4 seconds signal = ',num2str(accuracies(8))]);
disp(['accuracy for 4.5 seconds signal = ',num2str(accuracies(9))]);

%% Q3 
clear; clc; close all;
load('Ex3.mat');
fs = 256;

%% Q3 - Part a
clc;

Traindata = zeros(size(TrainData,1),size(TrainData,2),size(TrainData,3));
for i = 1:size(TrainData,3)
    Traindata(:,:,i) = TrainData(:,:,i)- mean(TrainData(:,:,i),2);
end

W = getCSPFilters(Traindata, TrainLabel, 2); % find CSP filters

filtered_data = applyCSPFilters(Traindata, W); % apply csp filters on data

 
trials_number = [65 75]; % select some trials
part_of_trials = zeros(2,(trials_number(2)-trials_number(1))*256);
for i=trials_number(1):trials_number(2)-1
    part_of_trials(1,((i-trials_number(1))*256+1):(i-trials_number(1)+1)*256) = filtered_data(1,:,i); 
    part_of_trials(2,((i-trials_number(1))*256+1):(i-trials_number(1)+1)*256) = filtered_data(2,:,i); 
end

t = trials_number(1):1/fs:trials_number(2)-1/fs;
subplot(2,1,1);
plot(t,part_of_trials(1,:));
grid on;
ylim([-15 15]);
xlabel('Time(s)','interpreter','latex');
title(['Labels: ',num2str(TrainLabel(trials_number(1):trials_number(2)))],'interpreter','latex');
subplot(2,1,2);
plot(t,part_of_trials(2,:));
ylim([-15 15]);
xlabel('Time(s)','interpreter','latex');
grid on;
%% Q3 - Part b
clc; close all;
% import data for AllElectrodes.mat file
[X,Y,labels]=importElectrodes();
% based on data set description
selected_electrodes = [37,7,5,38,40,42,10,47,45,15,13,48,50,52,18,32,55,23,22,21,20,31,57,58,59,60,26,63,27,64]; 
X = X(selected_electrodes);
Y = Y(selected_electrodes);
labels = labels(selected_electrodes);

% plot the spatial filters
figure;
subplot(1,2,1);
plottopomap(X,Y,labels,abs(W(:,1)));
title('Filter number 1','FontSize',14, 'Interpreter','latex');
subplot(1,2,2);
plottopomap(X,Y,labels,abs(W(:,2)));
title('Filter number 2','FontSize',14, 'Interpreter','latex');

%% Q3 - Part c
clc; close all;


num_folds = 4;
Traindata = Traindata(:,:,1:164);
TrainLabel = TrainLabel(1:164);
average_accuracy = zeros(15, 1);

csp_filters = getCSPFilters(Traindata, TrainLabel, 30);

for n = 1:15    
    total_accuracy = 0;
    W = csp_filters(:,[1:n,30-n+1:30]);
    features = zeros(2*n,size(Traindata,3));
    for i=1:size(Traindata,3)
        filtered_data = W'*Traindata(:,:,i);
        features(:,i) = var(transpose(filtered_data));
    end
    for fold = 1:num_folds

        size_fold = size(Traindata,3)/num_folds;
        val_idx = (fold-1)*size_fold+1:fold*size_fold;
        train_idx = setdiff(1:size(Traindata,3),val_idx);
        val_data = features(:,val_idx);
        val_label = TrainLabel(val_idx);
        train_data = features(:,train_idx);
        train_label = TrainLabel(train_idx);

        SVMModel = fitcsvm(train_data',train_label');
        predicted_label = predict(SVMModel,val_data');

        accuracy = sum(predicted_label == val_label')/(size(val_data,2));
        total_accuracy = total_accuracy + accuracy;
    end
    average_accuracy(n,1) = total_accuracy/num_folds;      
end

[~, optimal_num_filters] = max(average_accuracy);
% Display results
disp(['Optimal number of CSP filters: ' num2str(optimal_num_filters)]);
disp(['Average accuracy with optimal filters: ' num2str(average_accuracy(optimal_num_filters))]);

% Plot results
figure;
plot(1:15, average_accuracy, 'o-','LineWidth',1);
title('Cross-Validation Accuracy vs. Number of CSP Filters','Interpreter','latex','Fontsize',14);
xlabel('Number of CSP Filters','Interpreter','latex');
ylabel('Average Accuracy','Interpreter','latex');
%% Q3 - Part d
clc;
Testdata = zeros(size(TestData,1),size(TestData,2),size(TestData,3));
for i = 1:size(TestData,3)
    Testdata(:,:,i) = TestData(:,:,i)- mean(TestData(:,:,i),2);
end
n = optimal_num_filters;
W_optimal = W(:,[1:n,30-n+1:30]);
features_test = zeros(2*n,size(Testdata,3));
for i=1:size(Testdata,3)
    filtered_data = W_optimal'*Testdata(:,:,i);
    features_test(:,i) = var(transpose(filtered_data));
end
features_train = zeros(2*n,size(Traindata,3));
for i=1:size(Traindata,3)
    filtered_data = W_optimal'*Traindata(:,:,i);
    features_train(:,i) = var(transpose(filtered_data));
end
SVMModel = fitcsvm(features_train',TrainLabel);
predicted_labels = predict(SVMModel,features_test');
save("TestLabel.mat","predicted_labels");

%% Functions
function plottopomap(elocsX,elocsY,elabels,data)

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX)-.2,max(elocsX)+.25,interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);
% Interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsY,elocsX,data);
topodata = interpFunction(gridX,gridY);

% plot map
contourf(interpY,interpX,topodata);
hold on
scatter(elocsY,elocsX,10,'ro','filled');
for i=1:length(elocsX)
    text(elocsY(i),elocsX(i),elabels(i))
end
set(gca,'xtick',[])
set(gca,'ytick',[])
end
function [X,Y,labels] = importElectrodes()
    X = [80.7840137690914;68.6910763510323;76.1527667684846;59.9127302448179;57.5510633930990;54.0378881132512;49.8713779489202;
        26.2075035325936;28.7628234353576;30.9552849531915;32.4361838987395;2.11480422795274e-15;3.86812533613566e-15;
        4.94950482941819e-15;5.17649253748256e-15;-26.2075035325936;-28.7628234353576;-30.9552849531915;-32.4361838987395;
        -59.9127302448179;-57.5510633930990;-54.0378881132511;-49.8713779489202;-44.4841053285474;-68.6910763510323;-76.1527667684845;
        -80.7840137690914;-77.6332510774614;-84.9812336134463;-79.0255388591416;-60.7384809484625;-32.9278836352560;84.9812336134463;
        80.7840137690914;68.7208994216315;76.1527667684846;79.0255388591416;60.7384809484625;59.8744127660118;57.5840261068105;
        54.0263340465386;49.9265268118817;26.2075035325936;28.7628234353576;30.9552849531915;32.4361838987395;32.9278836352560;
        5.20474889637625e-15;2.11920249382479e-15;3.86788221025119e-15;4.94950482941819e-15;5.17649253748256e-15;-26.2847718099742;
        -28.7628234353576;-30.9552849531915;-32.4361838987395;-59.8744127660117;-57.5840261068105;-54.0263340465386;-49.9265268118817;
        -44.4841053285474;-68.7208994216315;-76.1527667684845;-80.7840137690914];
    Y = [26.1330144040702;49.7094313148880;31.4827967984807;26.0420933899754;48.2004273175388;63.0582218645482;68.4233350269540;80.4100143118706;
        76.2473645099531;59.2749781760892;32.3513771312283;34.5373740318457;63.1712807125907;80.8315480490248;84.5385386396573;80.4100143118706;
        76.2473645099531;59.2749781760892;32.3513771312283;26.0420933899754;48.2004273175389;63.0582218645482;68.4233350269539;59.7082740762203;
        49.7094313148880;31.4827967984807;26.1330144040702;-9.50733124394159e-15;-1.04071995732300e-14;-9.67783732147425e-15;-7.43831862786072e-15;
        -4.03250272966127e-15;0;-26.1330144040702;-49.6689040281160;-31.4827967984807;0;0;-26.0254380421476;-48.1425964684523;-63.0447391225751;
        -68.3835902976096;-80.4100143118706;-76.2473645099531;-59.2749781760892;-32.3513771312283;0;0;-34.6092031645412;-63.1673101655785;-80.8315480490248;
        -84.5385386396573;-80.3851021196760;-76.2473645099531;-59.2749781760892;-32.3513771312283;-26.0254380421476;-48.1425964684523;-63.0447391225751;-68.3835902976096;
        -59.7082740762203;-49.6689040281160;-31.4827967984807;-26.1330144040702];

    labels = {'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';
        'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';
        'P8';'P10';'PO8';'PO4';'O2'};
end

function [train_data, train_labels, val_data, val_labels] = splitDataForCrossValidation(data, labels, num_folds, fold)
    % Split data into training and validation sets for cross-validation
    fold_size = size(data, 3) / num_folds;
    
    val_start_idx = round((fold - 1) * fold_size) + 1;
    val_end_idx = round(fold * fold_size);
    
    val_data = data(:,:,val_start_idx:val_end_idx);
    val_labels = labels(val_start_idx:val_end_idx);
    
    train_data = cat(3, data(:,:,1:val_start_idx-1), data(:,:,val_end_idx+1:end));
    train_labels = cat(2, labels(1:val_start_idx-1), labels(val_end_idx+1:end));
end
function W = getCSPFilters(TrainData, TrainLabel, num_filters)
    C0 = zeros(30,30);
    C1 = zeros(30,30);
    counter = 0;
    num_data = size(TrainData,3);
    for i =1:num_data
        if (TrainLabel(i) == 0)
            C0 = C0 + TrainData(:,:,i)*TrainData(:,:,i)';
            counter = counter+1;
        else
            C1 = C1 + TrainData(:,:,i)*TrainData(:,:,i)';
        end
    end
    C0 = C0/counter;
    C1 = C1/(num_data-counter);
    
    [V,D] = eig(C0,C1);
    [~, ind] = sort(diag(D),'descend');
    V = V(:, ind);
    W = zeros(30,num_filters);
    for i=1:num_filters/2
        W(:,i) = V(:,i)/norm(V(:,i));
        W(:,num_filters-(i-1))=V(:,30-(i-1))/norm(V(:,30-(i-1)));
    end
    
end

function filtered_data = applyCSPFilters(data, filters)
    filtered_data = zeros(size(filters,2),size(data,2),size(data,3));

    for i = 1:size(data, 3)
        filtered_data(:,:,i) = filters' * data(:,:,i);
    end
end



