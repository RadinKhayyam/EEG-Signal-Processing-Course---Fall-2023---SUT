%% CHW2 - EEG Signal Processing - Radin Khayyam - 99101579
%% Q1 - Part a
clc; close all; clear;

signal = load('Ex1.mat').EEG_Sig;
fs = 200;
signal_length = length(signal)/fs;
t = 1/fs:1/fs:signal_length;
% part a
figure;
for i=1:3
    subplot(3,1,i);
    plot(t,signal(i,:));
    xlim([0 signal_length]);
    xlabel('Time(s)')
    ylabel('Amplitude')
    title(["original signal - Channel ",num2str(i)],'FontSize',14,'Interpreter','latex');
end

signal = signal - mean(signal,2); % change the mean to 0

%% Q1 - Part b
clc; close all;
figure();
scatter3(signal(1,:),signal(2,:),signal(3,:),'.');
title('Original Signal Scatter Plot','FontSize',14,'Interpreter','latex');
xlabel('channel 1');
ylabel('channel 2');
zlabel('channel 3');

%% Q1 - Part c
clc; close all;
covariance_matrix = cov(signal.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
whitened_signal = D * signal;

scatter3(signal(1,:),signal(2,:),signal(3,:),'.');
hold on;
for i=1:3
    plot3([0,U(1,i)*10],[0,U(2,i)*10],[0,U(2,i)*10],'LineWidth',2);
    hold on;
end
legend('data',['lambda1 = ',num2str(Lambda(1,1))],['lambda2 = ',num2str(Lambda(2,2))],['lambda3 = ',num2str(Lambda(3,3))])
title('Original Data and Principal Components','FontSize',14,'Interpreter','latex');

figure;
for i=1:3
    subplot(3,1,i);
    plot(t,whitened_signal(i,:));
    xlim([0 signal_length]);
    title(["Channel ",num2str(i)],'FontSize',14,'Interpreter','latex');
    xlabel('Time(s)')
    ylabel('Amplitude')
end

figure;
scatter3(whitened_signal(1,:),whitened_signal(2,:),whitened_signal(3,:),'.');
title('Whitened Signal Scatter Plot','FontSize',14,'Interpreter','latex');
new_covariance_matrix = cov(whitened_signal.');
disp(new_covariance_matrix);
%% Q1 - Part d
clc; close all;
[U,~,Lambda] = pca(signal.');
Lambda = Lambda.';
D = diag(Lambda.^(-1/2)) * U.';
whitened_signal = D * signal;

scatter3(signal(1,:),signal(2,:),signal(3,:),'.');
hold on;
for i=1:3
    plot3([0,U(1,i)*10],[0,U(2,i)*10],[0,U(2,i)*10],'LineWidth',2);
    hold on;
end
legend('data',['lambda1 = ',num2str(Lambda(1))],['lambda2 = ',num2str(Lambda(2))],['lambda3 = ',num2str(Lambda(3))])
title('Original Data and Principal Components','FontSize',14,'Interpreter','latex');

figure;
for i=1:3
    subplot(3,1,i);
    plot(t,whitened_signal(i,:));
    xlim([0 signal_length]);
    title(["Channel ",num2str(i)],'FontSize',14,'Interpreter','latex');
end

figure;
scatter3(whitened_signal(1,:),whitened_signal(2,:),whitened_signal(3,:),'.');
title('Whitened Signal Scatter Plot','FontSize',14,'Interpreter','latex');
new_covariance_matrix = cov(whitened_signal.');
%% Q1 - Part e
clc; close all;
[U,S,V] = svd(signal,0);

Lambda = diag(S.^2)/length(signal);
disp(Lambda);
disp(U);

D = diag(Lambda.^(-1/2)) * U.';
whitened_signal = D * signal;

[U_new,S_new,V_new] = svd(whitened_signal,0);

Lambda_new = diag(S_new.^2)/length(signal);
disp(Lambda_new);
disp(U_new);

