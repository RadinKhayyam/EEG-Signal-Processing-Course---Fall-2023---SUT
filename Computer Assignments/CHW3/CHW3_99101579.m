%% CHW3 - EEG signal processing - 99101579
%% Q1
clc; clear; close all;
load('Q1.mat');
fs = 100;
X = X_org - mean(X_org,2);
L = length(X_org);
num_channels = size(X_org,1);
%% Q1 - Part a - GEVD
clc;

% GEVD
taw=400;
V = GEVD_periodic(X,taw);
s_1_GEVD = V(:,1).'*X;
y_new = [s_1_GEVD;zeros(num_channels-1,L)];
x_1_hat_GEVD = inv(V.')*y_new;
disp('GEVD finished');

%DSS

% whitening
covariance_matrix = cov(X.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
z = D * X;
% algorithm
w_p_0 = randn(1,num_channels)';
w_p_new = DSS_periodic(w_p_0,z,taw);
w_p_previous = w_p_new;
flag=0;
counter = 1;
disp(['DSS step :',num2str(counter)]);
counter = counter+1;
while flag == 0
    w_p_new = DSS_periodic(w_p_previous,z,taw);
    if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
        flag=1;
    end
    clc;
    disp('GEVD finished');
    disp(['DSS step :',num2str(counter)]);
    w_p_previous = w_p_new;
    counter = counter + 1;
end
w_p_final = w_p_new;

s_1_DSS = w_p_final.'*z;
x_1_hat_DSS = inv(D)*w_p_final*s_1_DSS;

t = 1/fs: 1/fs:L/fs;
figure;
subplot(2,1,1);
plot(t,s_1_GEVD);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_1}(t),GEVD$','Fontsize',14,'Interpreter','latex');

subplot(2,1,2);
plot(t,s_1_DSS);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_1}(t),DSS$','Fontsize',14,'Interpreter','latex');

plotEEG(x_1_hat_GEVD);
title('$\hat{x_1}(t),GEVD$','Fontsize',14,'Interpreter','latex');

plotEEG(x_1_hat_DSS);
title('$\hat{x_1}(t),DSS$','Fontsize',14,'Interpreter','latex');

disp (['RRMSE of X1 estimation using GEVD= ',num2str(RRMSE(x_1_hat_GEVD,X1))]);
disp (['RRMSE of X1 Estimation using DSS= ',num2str(RRMSE(x_1_hat_DSS,X1))]);
%% Q1 - Part b
clc;

% finding best taw using auto correlation of signal
max_corr = 0;
min_taw = 300;
correlations = [];
for taw = 300:700
    shifted_signal = [zeros(8,taw),X(:,1:L-taw)];
    correlation = corrcoef(X,shifted_signal);
    correlation = correlation(1,2);
    if(correlation > max_corr)
        max_corr = correlation;
        min_taw = taw;
    end
end

% GEVD

V = GEVD_periodic(X,min_taw);
s_1_GEVD = V(:,1).'*X;
y_new = [s_1_GEVD;zeros(num_channels-1,L)];
x_1_hat_GEVD = inv(V.')*y_new;
disp('GEVD finished');

%DSS

% whitening
covariance_matrix = cov(X.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
z = D * X;
% algorithm
w_p_0 = randn(1,num_channels)';
w_p_new = DSS_periodic(w_p_0,z,min_taw);
w_p_previous = w_p_new;
flag=0;
counter = 1;
disp(['DSS step :',num2str(counter)]);
counter = counter+1;
while flag == 0
    w_p_new = DSS_periodic(w_p_previous,z,min_taw);
    if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
        flag=1;
    end
    clc;
    disp('GEVD finished');
    disp(['DSS step :',num2str(counter)]);
    w_p_previous = w_p_new;
    counter = counter + 1;
end
w_p_final = w_p_new;

s_1_DSS = w_p_final.'*z;
x_1_hat_DSS = inv(D)*w_p_final*s_1_DSS;

t = 1/fs: 1/fs:L/fs;
figure;
subplot(2,1,1);
plot(t,s_1_GEVD);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_1}(t),GEVD$','Fontsize',14,'Interpreter','latex');

subplot(2,1,2);
plot(t,s_1_DSS);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_1}(t),DSS$','Fontsize',14,'Interpreter','latex');

plotEEG(x_1_hat_GEVD);
title('$\hat{x_1}(t),GEVD$','Fontsize',14,'Interpreter','latex');

plotEEG(x_1_hat_DSS);
title('$\hat{x_1}(t),DSS$','Fontsize',14,'Interpreter','latex');

disp (['RRMSE of X1 estimation using GEVD= ',num2str(RRMSE(x_1_hat_GEVD,X1))]);
disp (['RRMSE of X1 Estimation using DSS= ',num2str(RRMSE(x_1_hat_DSS,X1))]);

%% Q1 - Part c
clc;
%GEVD
V = GEVD_non_stationary(X,T1);

s_2_GEVD = V(:,1).'*X;
y_new = [s_2_GEVD;zeros(7,L)];
x_2_hat_GEVD = inv(V.')*y_new;
disp('GEVD finished');
%DSS
% whitening
covariance_matrix = cov(X.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
z = D * X;
% algorithm
w_p_0 = randn(1,num_channels)';
w_p_new = DSS_non_stationary(w_p_0,z,T1);
w_p_previous = w_p_new;
flag=0;
counter = 1;
clc;
disp('GEVD finished');
disp(['DSS step :',num2str(counter)]);
counter = counter + 1;
while flag == 0
    w_p_new = DSS_non_stationary(w_p_previous,z,T1);
    if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
        flag=1;
    end
    w_p_previous = w_p_new;
    clc;
    disp('GEVD finished');
    disp(['DSS step :',num2str(counter)]);
    counter = counter + 1;
end
w_p_final = w_p_new;

s_2_DSS = w_p_final.'*z;
x_2_hat_DSS = inv(D)*w_p_final*s_2_DSS;

t = 1/fs: 1/fs:L/fs;
figure;
subplot(3,1,1);
plot(t,s_2_GEVD);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_2}(t),GEVD$','Fontsize',14,'Interpreter','latex');

subplot(3,1,2);
plot(t,s_2_DSS);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_2}(t),DSS$','Fontsize',14,'Interpreter','latex');

subplot(3,1,3);
stem(t,T1);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$T_1$','Fontsize',14,'Interpreter','latex');

plotEEG(x_2_hat_GEVD);
title('$\hat{x_2}(t),GEVD$','Fontsize',14,'Interpreter','latex');

plotEEG(x_2_hat_DSS);
title('$\hat{x_2}(t),DSS$','Fontsize',14,'Interpreter','latex');

disp (['RRMSE of X2 estimation using GEVD= ',num2str(RRMSE(x_2_hat_GEVD,X2))]);
disp (['RRMSE of X2 Estimation using DSS= ',num2str(RRMSE(x_2_hat_DSS,X2))]);


%% Q1 - Part d
clc;
%GEVD
T = T2;
for i=1:10
    V = GEVD_non_stationary(X,T);
    
    s_2_GEVD = V(:,1).'*X;
    y_new = [s_2_GEVD;zeros(7,L)];
    x_2_hat_GEVD = inv(V.')*y_new;
    T = find_spikes(x_2_hat_GEVD,3) | T2;
    disp(['iteration number: ',num2str(i)]);
end
disp('GEVD finished');

%DSS
% whitening
covariance_matrix = cov(X.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
z = D * X;

% algorithm
w_p_0 = randn(1,num_channels)';
w_p_new = DSS_non_stationary(w_p_0,z,T);
w_p_previous = w_p_new;
flag=0;
counter = 1;
clc;
disp('GEVD finished');
disp(['DSS step :',num2str(counter)]);
counter = counter + 1;
while flag == 0
    w_p_new = DSS_non_stationary(w_p_previous,z,T);
    if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
        flag=1;
    end
    w_p_previous = w_p_new;
    clc;
    disp('GEVD finished');
    disp(['DSS step :',num2str(counter)]);
    counter = counter + 1;
end
w_p_final = w_p_new;

s_2_DSS = w_p_final.'*z;
x_2_hat_DSS = inv(D)*w_p_final*s_2_DSS;

t = 1/fs: 1/fs:L/fs;
figure; 
subplot(3,1,1);
stem(t,T2);
xlabel('Time(s)','Interpreter','latex');
title('$T_2(t)$','Fontsize',14,'Interpreter','latex');

subplot(3,1,2);
stem(t,T);
xlabel('Time(s)','Interpreter','latex');
title('new estimated T(t)','Fontsize',14,'Interpreter','latex');

subplot(3,1,3);
stem(t,T1);
xlabel('Time(s)','Interpreter','latex');
title('$T_1(t)$','Fontsize',14,'Interpreter','latex');

figure;
subplot(3,1,1);
plot(t,s_2_GEVD);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_2}(t),GEVD$','Fontsize',14,'Interpreter','latex');

subplot(3,1,2);
plot(t,s_2_DSS);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_2}(t),DSS$','Fontsize',14,'Interpreter','latex');

subplot(3,1,3);
stem(t,T);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$T(t)$','Fontsize',14,'Interpreter','latex');

plotEEG(x_2_hat_GEVD);
title('$\hat{x_2}(t),GEVD$','Fontsize',14,'Interpreter','latex');

plotEEG(x_2_hat_DSS);
title('$\hat{x_2}(t),DSS$','Fontsize',14,'Interpreter','latex');

disp (['RRMSE of X2 estimation using GEVD= ',num2str(RRMSE(x_2_hat_GEVD,X2))]);
disp (['RRMSE of X2 Estimation using DSS= ',num2str(RRMSE(x_2_hat_DSS,X2))]);

%% Q1 - Part E
clc; 

% GEVD
desired_band =[10 15];
V = GEVD_spectral(X,desired_band,fs);
s_3_GEVD = V(:,1).'*X;
y_new = [s_3_GEVD;zeros(num_channels-1,L)];
x_3_hat_GEVD = inv(V.')*y_new;
disp('GEVD finished');
% DSS
% whitening
covariance_matrix = cov(X.');
[U,Lambda] = eig(covariance_matrix);
D = diag(diag(Lambda).^(-1/2)) * U.';
z = D * X;
% algorithm
w_p_0 = randn(1,num_channels)';
w_p_new = DSS_spectral(w_p_0,z,desired_band,fs);
w_p_previous = w_p_new;
flag=0;
counter = 1;
clc;
disp('GEVD finished');
disp(['DSS step :',num2str(counter)]);
counter = counter + 1;
while flag == 0
    w_p_new = DSS_spectral(w_p_previous,z,desired_band,fs);
    if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
        flag=1;
    end
    w_p_previous = w_p_new;
    clc;
    disp('GEVD finished');
    disp(['DSS step :',num2str(counter)]);
    counter = counter + 1;
end
w_p_final = w_p_new;

s_3_DSS = w_p_final.'*z;
x_3_hat_DSS = inv(D)*w_p_final*s_3_DSS;

t = 1/fs: 1/fs:L/fs;
figure;
subplot(2,1,1);
plot(t,s_3_GEVD);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_3}(t),GEVD$','Fontsize',14,'Interpreter','latex');

subplot(2,1,2);
plot(t,s_3_DSS);
ylim([-10 10]);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$\hat{s_3}(t),DSS$','Fontsize',14,'Interpreter','latex');

plotEEG(x_3_hat_GEVD);
title('$\hat{x_3}(t),GEVD$','Fontsize',14,'Interpreter','latex');

plotEEG(x_3_hat_DSS);
title('$\hat{x_3}(t),DSS$','Fontsize',14,'Interpreter','latex');

disp (['RRMSE of X3 estimation using GEVD= ',num2str(RRMSE(x_3_hat_GEVD,X3))]);
disp (['RRMSE of X3 Estimation using DSS= ',num2str(RRMSE(x_3_hat_DSS,X3))]);

%% Q2
clc; clear; close all;
load('Ex2.mat')
fs = 250;

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

%% Q2 - Part a
clc;
T1 = find_spikes(X_org,3);
L = length(X_org);

t = 1/fs:1/fs:L/fs;
plotEEG(X_org);
title('$x_{org}(t)$','Fontsize',14,'Interpreter','latex');

figure;
stem(t,T1);
xlabel('Time(s)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex')
title('$T_1$','Fontsize',15,'Interpreter','latex');
%% Q2 - Part b,c,d
clc;
signals = {X_SNR_10_N2,X_SNR_20_N2,X_SNR_10_N4,X_SNR_20_N4};
titles = {'SNR = -10dB , Noise 2$','SNR = -20dB , Noise 2$','SNR = -10dB , Noise 4$','SNR = -20dB , Noise 4$'};
error = cell(2,4);
L = length(signals{1});
num_sources = 2;
num_channels = size(signals{1},1);
for j=1:4
    clc;
    disp(['signal number ',num2str(j)]);
    % GEVD
    V = GEVD_non_stationary(signals{j},T1);
    S = V'*signals{j};
    y_new = [S(1:num_sources,:);zeros(num_channels - num_sources,L)];
    X_den_GEVD = inv(V')*y_new ;
    disp(['signal number ',num2str(j)]);
    disp('GEVD finished');
    %DSS
    % whitening
    covariance_matrix = cov(signals{j}.');
    [U,Lambda] = eig(covariance_matrix);
    D = diag(diag(Lambda).^(-1/2)) * U.';
    z = D * signals{j};
    S= zeros(num_sources,L);
    W = zeros(num_channels,num_sources);
    for k=1:num_sources
    
        w_p_0 = randn(1,num_channels)';
        w_p_new = DSS_non_stationary(w_p_0,z,T1);
        w_p_previous = w_p_new;
        flag=0;
        counter = 1;
        clc;
        disp(['signal number ',num2str(j)]);
        disp('GEVD finished');
        disp(['Extracting Source number ',num2str(k),' using DSS']);
        disp(['DSS step :',num2str(counter)]);
        while flag == 0
            w_p_new = DSS_non_stationary(w_p_previous,z,T1);
            if norm(w_p_new - w_p_previous)<0.01 || counter>=10 % convergence or step limit
                flag=1;
            end
            w_p_previous = w_p_new;
            clc;
            disp(['signal number ',num2str(j)]);
            disp('GEVD finished');
            disp(['Extracting Source number ',num2str(k),' using DSS']);
            disp(['DSS step :',num2str(counter)]);
            counter = counter + 1;
            
        end
        w_p_final = w_p_new;
        
        S(k,:) = w_p_final.'*z;
        W(:,k) = w_p_final;
    end
    X_den_DSS = inv(D)*W*S;
    
    t = 1/fs: 1/fs:L/fs;
    figure;
    subplot(4,2,1);
    plot(t,signals{j}(13,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title(['$X_{noisy,13}(t)',titles{j}],'Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,2);
    plot(t,signals{j}(24,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title(['$X_{noisy,24}(t)',titles{j}],'Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,3);
    plot(t,X_den_GEVD(13,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{denoised,13}(t) , GEVD$','Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,4);
    plot(t,X_den_GEVD(24,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{denoised,24}(t),GEVD$','Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,5);
    plot(t,X_den_DSS(13,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{denoised,13}(t) ,DSS$','Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,6);
    plot(t,X_den_DSS(24,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{denoised,24}(t),DSS$','Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,7);
    plot(t,X_org(13,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{orginal,13}(t)$','Fontsize',15,'Interpreter','latex');
    
    subplot(4,2,8);
    plot(t,X_org(24,:));
    xlim([0 L/fs]);
    xlabel('Time(s)','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex')
    title('$X_{orginal,24}(t)$','Fontsize',15,'Interpreter','latex');

    error{1,j}=RRMSE(X_den_GEVD,X_org);
    error{2,j}=RRMSE(X_den_DSS,X_org);

end
clc;
disp (['RRMSE for (Noise 2, SNR = -10db) using GEVD = ',num2str(error{1,1})]);
disp (['RRMSE for (Noise 2, SNR = -10db) using DSS = ',num2str(error{2,1})]);
disp (['RRMSE for (Noise 2, SNR = -20db) using GEVD = ',num2str(error{1,2})]);
disp (['RRMSE for (Noise 2, SNR = -20db) using DSS = ',num2str(error{2,2})]);
disp (['RRMSE for (Noise 4, SNR = -10db) using GEVD = ',num2str(error{1,3})]);
disp (['RRMSE for (Noise 4, SNR = -10db) using DSS = ',num2str(error{2,3})]);
disp (['RRMSE for (Noise 4, SNR = -20db) using GEVD = ',num2str(error{1,4})]);
disp (['RRMSE for (Noise 4, SNR = -20db) using DSS = ',num2str(error{2,4})]);

%% Functions
function plotEEG(X) 
offset = max(abs(X(:))) ;
feq = 250 ;
disp_eeg(X,offset,feq,[]);
end

function result = RRMSE(X_org,X_den)

temp_num = sum(sum((X_org-X_den).^2,2));
temp_den = sum(sum((X_org.^2),2));

result = sqrt(temp_num)/sqrt(temp_den);

end

function V = GEVD_non_stationary(X,T1)
    
    L = length(X);
    tmp = 0;
    for i=1:L
        if(T1(i)==1)
            tmp = tmp + X(:,i)*X(:,i)';
        end
    end
    
    C_x_tilde = tmp/(sum(T1));
    C_x = cov(X.');
    [V,D] = eig(C_x_tilde,C_x);
    [r, ind] = sort(diag(D),'descend');
    D=diag(r);
    V = V(:, ind);

    
end

function V = GEVD_periodic(X,taw)
    C_x = cov(X.');
    L = length(X);
    tmp = 0;
    for i=1:L-taw
        tmp = tmp + X(:,i)*X(:,i+taw).';
    end
    P_x = tmp/(L-taw);
    % P_x = (X_org(:,1:L-taw))*(X_org(:,taw+1:L)')/(L-taw);
    P_x_tilde = (1/2)*(P_x + P_x.');
    [V,D] = eig(P_x_tilde,C_x);
    [r, ind] = sort(diag(D),'descend');
    D=diag(r);
    V = V(:, ind);
end

function V = GEVD_spectral(X,desired_band,fs)
    C_x = cov(X.');
    L = length(X);
    num_channels = size(X,1);
    F_X = zeros(num_channels,L);
    for i=1:8
        F_X(i,:)= fftshift(fft(X(i,:)));
    end
    
    k = [L*desired_band(1)/fs, L*desired_band(2)/fs];
    tmp = 0;
    for i=k(1):k(2)
        tmp = tmp + F_X(:,i+(L/2+1))*F_X(:,i+(L/2+1))';
        tmp = tmp + F_X(:,-i+(L/2+1))*F_X(:,-i+(L/2+1))';
    end
    
    S_x = tmp / (2*length(k(1):k(2)));
    
    [V,D] = eig(S_x,C_x);
    [r, ind] = sort(diag(D),'descend');
    D=diag(r);
    V = V(:, ind); 
end

function w_p_new = DSS_periodic(w_p,z,taw)
    
    L=length(z);
    tmp = 0;
    r_p = w_p'*z;
    for i=1:(L/taw)
        tmp = tmp + r_p(1,1+(i-1)*400:400*i);
    end
    mean_period = tmp / (L/taw);
    r_p_new = zeros(1,L);
    for i=1:(L/taw)
        r_p_new(1,1+(i-1)*400:400*i) = mean_period;
    end
    tmp=0;
    for i=1:L
        tmp = tmp + z*r_p_new';
    end
    w_p_plus = tmp/L;
    w_p_new = w_p_plus / norm(w_p_plus);
    
end

function w_p_new = DSS_non_stationary(w_p,z,T1)
    
    L=length(z);
    r_p = w_p'*z;
    r_p_new = r_p.*T1;
    tmp=0;
    for i=1:L
        tmp = tmp + z*r_p_new';
    end
    w_p_plus = tmp/L;
    w_p_new = w_p_plus / norm(w_p_plus);
    
end

function w_p_new = DSS_spectral(w_p,z,fpass,fs)
    
    L=length(z);
    r_p = w_p'*z;
    r_p_new = bandpass(r_p,fpass,fs);
    tmp=0;
    for i=1:L
        tmp = tmp + z*r_p_new';
    end
    w_p_plus = tmp/L;
    w_p_new = w_p_plus / norm(w_p_plus);
    
end

function T = find_spikes(X,threshold_coef)
    num_channels = size(X,1);
    max_channel = max(abs(X)')';
    normalized_EEG = X./max_channel;
    mean_channels = sum(normalized_EEG)/num_channels;
    T = mean_channels>(max(mean_channels)/threshold_coef);
end

