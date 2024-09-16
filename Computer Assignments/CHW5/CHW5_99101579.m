%% EEG Signal Processing - Radin Khayyam - 99101579
%% CHW5
%% Part a
clear; close all; clc ;

load('ElecPosXYZ') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ; % Radius of diffetent layers
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3]; 
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;

scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Dipoles Locations','FontSize',14,'Interpreter','latex');
%% Part b
clc;
Elec_numbers = 21;
Elec_positions=zeros(3,Elec_numbers);
Elec_names=strings(1,Elec_numbers);

for i=1:Elec_numbers
    Elec_positions(:,i)=ElecPos{1, i}.XYZ;
    Elec_names(1,i)=ElecPos{1, i}.Name;
end

scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Dipoles and Electrodes Locations','FontSize',14,'Interpreter','latex');
%% Part c - shallow dipole
clc; 

idx_slct = 1299;
loc_slct = LocMat(:,idx_slct);
dir_slct = loc_slct/sqrt(sum(loc_slct.^2));

figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g');
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('One dipole vector','FontSize',14,'Interpreter','latex');

%% Part d - shallow dipole
clc;
Interictal=load('Interictal.mat') ;
interictal_signal=Interictal.Interictal;
signal = interictal_signal(1,:);
Q = zeros(3951,length(signal));
Q((idx_slct-1)*3+1:idx_slct*3,:) = [signal*dir_slct(1);signal*dir_slct(2);signal*dir_slct(3)];
M = GainMat*Q;
disp_eeg(M,max(abs(M(:))),256,Elec_names,'EEG signal');

%% Part e - shallow dipole
clc;

peaks = cell(2, 21);

for i=1:21
    [pks,locs] = findpeaks(M(i,:),MinPeakDistance=10,SortStr="descend");
    locs(pks < max(pks/2)) = [];
    pks(pks < max(pks/2)) = [];
    peaks{1,i} = pks;
    peaks{2,i} = locs;
end

figure()
for i=1:21
    subplot( 6 ,4 ,i);
    plot( M(i,:));
    hold on;
    plot(peaks{2,i},peaks{1,i},"o",'LineWidth',1);
    xlabel('smaple','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex');
    title(Elec_names(1,i),'Fontsize',14,'Interpreter','latex')
end

mean_peaks = zeros(21,1);
for i=1:21
    window = (cell2mat(peaks(2,:))-3):(cell2mat(peaks(2,:))+3);
    mean_peaks(i) = mean(abs(M(i,window)));
end



normalized_mean_peaks = mean_peaks / max(mean_peaks);
for i=1:Elec_numbers
    Elec_positions(:,i)=ElecPos{1, i}.XYZ;
    Elec_names(1,i)=ElecPos{1, i}.Name;
end
colors = parula(100);
tmp = sort(normalized_mean_peaks);
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
for i=1:21
    color_idx = ceil(normalized_mean_peaks(i)*100);
    scatter3(ModelParams.R(3)*Elec_positions(1,i),ModelParams.R(3)*Elec_positions(2,i),ModelParams.R(3)*Elec_positions(3,i),'MarkerEdgeColor','k','MarkerFaceColor',colors((color_idx),:));
    hold on
end

text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Dipoles and Electrodes Locations','FontSize',14,'Interpreter','latex');
colorbar;


%% Part f - shallow dipole
clc;
figure();
Display_Potential_3D(ModelParams.R(3),mean_peaks);

%% Part g - shallow dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
M = mean_peaks;
G = GainMat;
alpha=0.1;
Q_MNE=G.'*inv(G*G.'+alpha*eye(21))*M;

% WMNE
Omega = zeros(1317,1317);
for i = 1:1317
    tmp = 0;
    for j = 1:21
        tmp = tmp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    Omega(i,i) = sqrt(tmp);
end
W = kron(Omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; 

% LORETA
A = zeros(1317,1317);
d = 1;
P=1317;
for i = 1:P
    for j = 1:P
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end
A_tmp = A*ones(P);
A0 = inv(A_tmp.*eye(size(A_tmp)))*A;
B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*P,3*P));
W = kron(Omega, eye(3)) * B'*B*kron(Omega, eye(3));
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M;

% sLORETA
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_sLORETA = zeros(3*P, 1);
for i = 1:P
    Q_sLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%% Part h - shallow dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
moment_MNE=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    moment_MNE(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
idx_est_MNE = min(find(moment_MNE==max(moment_MNE)));
loc_est_MNE = LocMat(:,idx_est_MNE);
dir_est_MNE=Q_MNE((idx_est_MNE-1)*3+1:idx_est_MNE*3)./sqrt(sum(Q_MNE((idx_est_MNE-1)*3+1:idx_est_MNE*3).^2,'all'));

% WMNE
moment_WMNE=zeros(1,length(Q_WMNE)/3);
for i = 1:length(Q_WMNE)/3
    moment_WMNE(1,i) = sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
idx_est_WMNE = min(find(moment_WMNE==max(moment_WMNE)));
loc_est_WMNE = LocMat(:,idx_est_WMNE);
dir_est_WMNE=Q_MNE((idx_est_WMNE-1)*3+1:idx_est_WMNE*3)./sqrt(sum(Q_WMNE((idx_est_WMNE-1)*3+1:idx_est_WMNE*3).^2,'all'));

% LORETA
moment_LORETA=zeros(1,length(Q_LORETA)/3);
for i=1:length(Q_LORETA)/3
    moment_LORETA(1,i)=sqrt(sum(Q_LORETA((i-1)*3+1:i*3).^2,'all'));
end
idx_est_LORETA = min(find(moment_LORETA==max(moment_LORETA)));
loc_est_LORETA = LocMat(:,idx_est_LORETA);
dir_est_LORETA=Q_LORETA((idx_est_LORETA-1)*3+1:idx_est_LORETA*3)./sqrt(sum(Q_LORETA((idx_est_LORETA-1)*3+1:idx_est_LORETA*3).^2,'all'));

% sLORETA
moment_sLORETA=zeros(1,length(Q_sLORETA)/3);
for i=1:length(Q_sLORETA)/3
    moment_sLORETA(1,i)=sqrt(sum(Q_sLORETA((i-1)*3+1:i*3).^2,'all'));
end
idx_est_sLORETA = min(find(moment_sLORETA==max(moment_sLORETA)));
loc_est_sLORETA = LocMat(:,idx_est_sLORETA);
dir_est_sLORETA=Q_sLORETA((idx_est_sLORETA-1)*3+1:idx_est_sLORETA*3)./sqrt(sum(Q_sLORETA((idx_est_sLORETA-1)*3+1:idx_est_sLORETA*3).^2,'all'));

% Plot
figure;
subplot(2,2,1);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_est_MNE),LocMat(2,idx_est_MNE),LocMat(3,idx_est_MNE),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_MNE),LocMat(2,idx_est_MNE),LocMat(3,idx_est_MNE),'r','Linewidth',1);
hold on
plot3( [ loc_est_MNE(1,1),loc_est_MNE(1,1)+dir_est_MNE(1,1) ] ,[ loc_est_MNE(2,1),loc_est_MNE(2,1)+dir_est_MNE(2,1)] , [loc_est_MNE(3,1),loc_est_MNE(3,1)+dir_est_MNE(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('MNE method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,2);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_WMNE),LocMat(2,idx_est_WMNE),LocMat(3,idx_est_WMNE),'r','Linewidth',1);
hold on
plot3( [ loc_est_WMNE(1,1),loc_est_WMNE(1,1)+dir_est_WMNE(1,1) ] ,[ loc_est_WMNE(2,1),loc_est_WMNE(2,1)+dir_est_WMNE(2,1)] , [loc_est_WMNE(3,1),loc_est_WMNE(3,1)+dir_est_WMNE(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('WMNE method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,3);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_LORETA),LocMat(2,idx_est_LORETA),LocMat(3,idx_est_LORETA),'r','Linewidth',1);
hold on
plot3( [ loc_est_LORETA(1,1),loc_est_LORETA(1,1)+dir_est_LORETA(1,1) ] ,[ loc_est_LORETA(2,1),loc_est_LORETA(2,1)+dir_est_LORETA(2,1)] , [loc_est_LORETA(3,1),loc_est_LORETA(3,1)+dir_est_LORETA(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('LORETA method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,4);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_sLORETA),LocMat(2,idx_est_sLORETA),LocMat(3,idx_est_sLORETA),'r','Linewidth',1);
hold on
plot3( [ loc_est_sLORETA(1,1),loc_est_sLORETA(1,1)+dir_est_LORETA(1,1) ] ,[ loc_est_sLORETA(2,1),loc_est_sLORETA(2,1)+dir_est_LORETA(2,1)] , [loc_est_sLORETA(3,1),loc_est_sLORETA(3,1)+dir_est_sLORETA(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('sLORETA method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

figure;

subplot(2,2,1);
plot3([0,dir_est_MNE(1)],[0,dir_est_MNE(2)],[0,dir_est_MNE(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,2);
plot3([0,dir_est_WMNE(1)],[0,dir_est_WMNE(2)],[0,dir_est_WMNE(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,3);
plot3([0,dir_est_LORETA(1)],[0,dir_est_LORETA(2)],[0,dir_est_LORETA(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,4);
plot3([0,dir_est_sLORETA(1)],[0,dir_est_sLORETA(2)],[0,dir_est_sLORETA(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

% print results
disp(['Estimated location of MNE method: ', num2str(loc_est_MNE')]);
disp(['Estimated location of WMNE method: ', num2str(loc_est_WMNE')]);
disp(['Estimated location of LORETA method: ', num2str(loc_est_LORETA')]);
disp(['Estimated location of sLORETA method: ', num2str(loc_est_sLORETA')]);
disp(['Actual location: ', num2str([LocMat(1,idx_slct) LocMat(2,idx_slct) LocMat(3,idx_slct)])]);
disp('----------------------------------------------');
disp(['Estimated direction of MNE method: ', num2str(dir_est_MNE')]);
disp(['Estimated direction of WMNE method: ', num2str(dir_est_WMNE')]);
disp(['Estimated direction of LORETA method: ', num2str(dir_est_LORETA')]);
disp(['Estimated direction of sLORETA method: ', num2str(dir_est_sLORETA')]);
disp(['Actual direction: ', num2str(dir_slct')]);
disp('----------------------------------------------');
%% Part i - shallow dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
loc_error_MNE=sqrt(sum((loc_slct-loc_est_MNE).^2,'all'));
dir_error_MNE = sqrt(sum((dir_est_MNE-dir_slct).^2));
% WMNE
loc_error_WMNE=sqrt(sum((loc_slct-loc_est_WMNE).^2,'all'));
dir_error_WMNE = sqrt(sum((dir_est_WMNE-dir_slct).^2));
%LORETA
loc_error_LORETA=sqrt(sum((loc_est_LORETA-loc_slct).^2,'all'));
dir_error_LORETA = sqrt(sum((dir_est_LORETA-dir_slct).^2));
%sLORETA
loc_error_sLORETA=sqrt(sum((loc_est_sLORETA-loc_slct).^2,'all'));
dir_error_sLORETA = sqrt(sum((dir_est_sLORETA-dir_slct).^2));

% print results
disp(['MSE of location in MNE method: = ',num2str(loc_error_MNE)]);
disp(['MSE of location in WMNE method: = ',num2str(loc_error_WMNE)]);
disp(['MSE of location in LORETA method: = ',num2str(loc_error_LORETA)]);
disp(['MSE of location in sLORETA method: = ',num2str(loc_error_sLORETA)]);
disp('----------------------------------------------------')
disp(['MSE of direction in MNE method: = ',num2str(dir_error_MNE)]);
disp(['MSE of direction in WMNE method: = ',num2str(dir_error_WMNE)]);
disp(['MSE of direction in LORETA method: = ',num2str(dir_error_LORETA)]);
disp(['MSE of direction in sLORETA method: = ',num2str(dir_error_sLORETA)]);
%% Part c - deep dipole
clc; 

idx_slct = 290;
loc_slct = LocMat(:,idx_slct);
dir_slct = loc_slct/sqrt(sum(loc_slct.^2));

figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g');
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('One dipole vector','FontSize',14,'Interpreter','latex');

%% Part d - deep dipole
clc;
Interictal=load('Interictal.mat') ;
interictal_signal=Interictal.Interictal;
signal = interictal_signal(1,:);
Q = zeros(3951,length(signal));
Q((idx_slct-1)*3+1:idx_slct*3,:) = [signal*dir_slct(1);signal*dir_slct(2);signal*dir_slct(3)];
M = GainMat*Q;
disp_eeg(M,max(abs(M(:))),256,Elec_names,'EEG signal');

%% Part e - deep dipole
clc;

peaks = cell(2, 21);

for i=1:21
    [pks,locs] = findpeaks(M(i,:),MinPeakDistance=10,SortStr="descend");
    locs(pks < max(pks/2)) = [];
    pks(pks < max(pks/2)) = [];
    peaks{1,i} = pks;
    peaks{2,i} = locs;
end

figure()
for i=1:21
    subplot( 6 ,4 ,i);
    plot( M(i,:));
    hold on;
    plot(peaks{2,i},peaks{1,i},"o",'LineWidth',1);
    xlabel('smaple','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex');
    title(Elec_names(1,i),'Fontsize',14,'Interpreter','latex')
end

mean_peaks = zeros(21,1);
for i=1:21
    window = (cell2mat(peaks(2,:))-3):(cell2mat(peaks(2,:))+3);
    mean_peaks(i) = mean(abs(M(i,window)));
end



normalized_mean_peaks = mean_peaks / max(mean_peaks);
for i=1:Elec_numbers
    Elec_positions(:,i)=ElecPos{1, i}.XYZ;
    Elec_names(1,i)=ElecPos{1, i}.Name;
end
colors = parula(100);
tmp = sort(normalized_mean_peaks);
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
for i=1:21
    color_idx = ceil(normalized_mean_peaks(i)*100);
    scatter3(ModelParams.R(3)*Elec_positions(1,i),ModelParams.R(3)*Elec_positions(2,i),ModelParams.R(3)*Elec_positions(3,i),'MarkerEdgeColor','k','MarkerFaceColor',colors((color_idx),:));
    hold on
end

text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Dipoles and Electrodes Locations','FontSize',14,'Interpreter','latex');
colorbar;


%% Part f - deep dipole
clc;
figure();
Display_Potential_3D(ModelParams.R(3),mean_peaks);

%% Part g - deep dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
M = mean_peaks;
G = GainMat;
alpha=0.1;
Q_MNE=G.'*inv(G*G.'+alpha*eye(21))*M;

% WMNE
Omega = zeros(1317,1317);
for i = 1:1317
    tmp = 0;
    for j = 1:21
        tmp = tmp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    Omega(i,i) = sqrt(tmp);
end
W = kron(Omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; 

% LORETA
d = 1;
P=1317;
A = zeros(P,P);
for i = 1:P
    for j = 1:P
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end
A_tmp = A*ones(P);
A0 = inv(A_tmp.*eye(size(A_tmp)))*A;
B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*P,3*P));
W = kron(Omega, eye(3)) * B'*B*kron(Omega, eye(3));
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M;

% sLORETA
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_sLORETA = zeros(3*P, 1);
for i = 1:P
    Q_sLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%% Part h - deep dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
moment_MNE=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    moment_MNE(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
idx_est_MNE = min(find(moment_MNE==max(moment_MNE)));
loc_est_MNE = LocMat(:,idx_est_MNE);
dir_est_MNE=Q_MNE((idx_est_MNE-1)*3+1:idx_est_MNE*3)./sqrt(sum(Q_MNE((idx_est_MNE-1)*3+1:idx_est_MNE*3).^2,'all'));

% WMNE
moment_WMNE=zeros(1,length(Q_WMNE)/3);
for i = 1:length(Q_WMNE)/3
    moment_WMNE(1,i) = sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
idx_est_WMNE = min(find(moment_WMNE==max(moment_WMNE)));
loc_est_WMNE = LocMat(:,idx_est_WMNE);
dir_est_WMNE=Q_MNE((idx_est_WMNE-1)*3+1:idx_est_WMNE*3)./sqrt(sum(Q_WMNE((idx_est_WMNE-1)*3+1:idx_est_WMNE*3).^2,'all'));

% LORETA
moment_LORETA=zeros(1,length(Q_LORETA)/3);
for i=1:length(Q_LORETA)/3
    moment_LORETA(1,i)=sqrt(sum(Q_LORETA((i-1)*3+1:i*3).^2,'all'));
end
idx_est_LORETA = min(find(moment_LORETA==max(moment_LORETA)));
loc_est_LORETA = LocMat(:,idx_est_LORETA);
dir_est_LORETA=Q_LORETA((idx_est_LORETA-1)*3+1:idx_est_LORETA*3)./sqrt(sum(Q_LORETA((idx_est_LORETA-1)*3+1:idx_est_LORETA*3).^2,'all'));

% sLORETA
moment_sLORETA=zeros(1,length(Q_sLORETA)/3);
for i=1:length(Q_sLORETA)/3
    moment_sLORETA(1,i)=sqrt(sum(Q_sLORETA((i-1)*3+1:i*3).^2,'all'));
end
idx_est_sLORETA = min(find(moment_sLORETA==max(moment_sLORETA)));
loc_est_sLORETA = LocMat(:,idx_est_sLORETA);
dir_est_sLORETA=Q_sLORETA((idx_est_sLORETA-1)*3+1:idx_est_sLORETA*3)./sqrt(sum(Q_sLORETA((idx_est_sLORETA-1)*3+1:idx_est_sLORETA*3).^2,'all'));

% Plot
figure;
subplot(2,2,1);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_MNE),LocMat(2,idx_est_MNE),LocMat(3,idx_est_MNE),'r','Linewidth',1);
hold on
plot3( [ loc_est_MNE(1,1),loc_est_MNE(1,1)+dir_est_MNE(1,1) ] ,[ loc_est_MNE(2,1),loc_est_MNE(2,1)+dir_est_MNE(2,1)] , [loc_est_MNE(3,1),loc_est_MNE(3,1)+dir_est_MNE(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('MNE method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,2);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_WMNE),LocMat(2,idx_est_WMNE),LocMat(3,idx_est_WMNE),'r','Linewidth',1);
hold on
plot3( [ loc_est_WMNE(1,1),loc_est_WMNE(1,1)+dir_est_WMNE(1,1) ] ,[ loc_est_WMNE(2,1),loc_est_WMNE(2,1)+dir_est_WMNE(2,1)] , [loc_est_WMNE(3,1),loc_est_WMNE(3,1)+dir_est_WMNE(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('WMNE method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,3);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_LORETA),LocMat(2,idx_est_LORETA),LocMat(3,idx_est_LORETA),'r','Linewidth',1);
hold on
plot3( [ loc_est_LORETA(1,1),loc_est_LORETA(1,1)+dir_est_LORETA(1,1) ] ,[ loc_est_LORETA(2,1),loc_est_LORETA(2,1)+dir_est_LORETA(2,1)] , [loc_est_LORETA(3,1),loc_est_LORETA(3,1)+dir_est_LORETA(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('LORETA method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

subplot(2,2,4);
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
scatter3(LocMat(1,idx_slct),LocMat(2,idx_slct),LocMat(3,idx_slct),'g','Linewidth',1);
hold on
plot3( [ loc_slct(1,1),loc_slct(1,1)+dir_slct(1,1) ] ,[ loc_slct(2,1),loc_slct(2,1)+dir_slct(2,1)] , [loc_slct(3,1),loc_slct(3,1)+dir_slct(3,1)] ,'Color','black','LineWidth',2)
hold on
scatter3(LocMat(1,idx_est_sLORETA),LocMat(2,idx_est_sLORETA),LocMat(3,idx_est_sLORETA),'r','Linewidth',1);
hold on
plot3( [ loc_est_sLORETA(1,1),loc_est_sLORETA(1,1)+dir_est_LORETA(1,1) ] ,[ loc_est_sLORETA(2,1),loc_est_sLORETA(2,1)+dir_est_LORETA(2,1)] , [loc_est_sLORETA(3,1),loc_est_sLORETA(3,1)+dir_est_sLORETA(3,1)] ,'Color','black','LineWidth',2)
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('sLORETA method','FontSize',14,'Interpreter','latex');
legend('dipoles','electrodes','selected dipole','direction','estimated dipole');

figure;

subplot(2,2,1);
plot3([0,dir_est_MNE(1)],[0,dir_est_MNE(2)],[0,dir_est_MNE(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,2);
plot3([0,dir_est_WMNE(1)],[0,dir_est_WMNE(2)],[0,dir_est_WMNE(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,3);
plot3([0,dir_est_LORETA(1)],[0,dir_est_LORETA(2)],[0,dir_est_LORETA(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

subplot(2,2,4);
plot3([0,dir_est_sLORETA(1)],[0,dir_est_sLORETA(2)],[0,dir_est_sLORETA(3)],'LineWidth',2);
hold on;
plot3([0,dir_slct(1)],[0,dir_slct(2)],[0,dir_slct(3)],'LineWidth',2);
legend('Estimated direction','Original direction');

% print results
disp(['Estimated location of MNE method: ', num2str(loc_est_MNE')]);
disp(['Estimated location of WMNE method: ', num2str(loc_est_WMNE')]);
disp(['Estimated location of LORETA method: ', num2str(loc_est_LORETA')]);
disp(['Estimated location of sLORETA method: ', num2str(loc_est_sLORETA')]);
disp(['Actual location: ', num2str([LocMat(1,idx_slct) LocMat(2,idx_slct) LocMat(3,idx_slct)])]);
disp('----------------------------------------------');
disp(['Estimated direction of MNE method: ', num2str(dir_est_MNE')]);
disp(['Estimated direction of WMNE method: ', num2str(dir_est_WMNE')]);
disp(['Estimated direction of LORETA method: ', num2str(dir_est_LORETA')]);
disp(['Estimated direction of sLORETA method: ', num2str(dir_est_sLORETA')]);
disp(['Actual direction: ', num2str(dir_slct')]);
disp('----------------------------------------------');
%% Part i - deep dipole - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
loc_error_MNE=sqrt(sum((loc_slct-loc_est_MNE).^2,'all'));
dir_error_MNE = sqrt(sum((dir_est_MNE-dir_slct).^2));
% WMNE
loc_error_WMNE=sqrt(sum((loc_slct-loc_est_WMNE).^2,'all'));
dir_error_WMNE = sqrt(sum((dir_est_WMNE-dir_slct).^2));
%LORETA
loc_error_LORETA=sqrt(sum((loc_est_LORETA-loc_slct).^2,'all'));
dir_error_LORETA = sqrt(sum((dir_est_LORETA-dir_slct).^2));
%sLORETA
loc_error_sLORETA=sqrt(sum((loc_est_sLORETA-loc_slct).^2,'all'));
dir_error_sLORETA = sqrt(sum((dir_est_sLORETA-dir_slct).^2));

% print results
disp(['MSE of location in MNE method: = ',num2str(loc_error_MNE)]);
disp(['MSE of location in WMNE method: = ',num2str(loc_error_WMNE)]);
disp(['MSE of location in LORETA method: = ',num2str(loc_error_LORETA)]);
disp(['MSE of location in sLORETA method: = ',num2str(loc_error_sLORETA)]);
disp('----------------------------------------------------')
disp(['MSE of direction in MNE method: = ',num2str(dir_error_MNE)]);
disp(['MSE of direction in WMNE method: = ',num2str(dir_error_WMNE)]);
disp(['MSE of direction in LORETA method: = ',num2str(dir_error_LORETA)]);
disp(['MSE of direction in sLORETA method: = ',num2str(dir_error_sLORETA)]);
%% Part j

%% part k
clc;

center_dipole_idx = 550;
patch_slct =[];

for i = 1:length(LocMat)
    if(norm(LocMat(:,i) - LocMat(:,center_dipole_idx)) <= 1.5)
        patch_slct = [patch_slct i];
    end
end

loc_patch = LocMat(:, patch_slct);
dir_slct = loc_patch./sqrt(sum(loc_patch.^2));

Elec_numbers = 21;
Elec_positions=zeros(3,Elec_numbers);
Elec_names=strings(1,Elec_numbers);

for i=1:Elec_numbers
    Elec_positions(:,i)=ElecPos{1, i}.XYZ;
    Elec_names(1,i)=ElecPos{1, i}.Name;
end

scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
scatter3(ModelParams.R(3)*Elec_positions(1,:),ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),'MarkerEdgeColor','k','MarkerFaceColor','red');
hold on
text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
hold on
for i=1:length(patch_slct)
    scatter3(LocMat(1,patch_slct(i)),LocMat(2,patch_slct(i)),LocMat(3,patch_slct(i)),'g','Linewidth',1);
    hold on
    plot3( [ loc_patch(1,i),loc_patch(1,i)+dir_slct(1,i) ] ,[ loc_patch(2,i),loc_patch(2,i)+dir_slct(2,i)] , [loc_patch(3,i),loc_patch(3,i)+dir_slct(3,i)] ,'Color','black','LineWidth',2)
    hold on
end

xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Selected patch with r = 1.5cm around dipole number 1142','FontSize',14,'Interpreter','latex');
%% Part l
clc;
Interictal=load('Interictal.mat') ;
interictal_signal=Interictal.Interictal;
Q = zeros(3951,size(interictal_signal,2));
for i = 1:length(patch_slct)
    signal = interictal_signal(i,:);
    dipole_idx = patch_slct(i);
    loc_dipole = LocMat(:,dipole_idx);
    dir_dipole = loc_dipole/sqrt(sum(loc_dipole.^2));
    Q((dipole_idx-1)*3+1:dipole_idx*3,:) = [signal*dir_dipole(1);signal*dir_dipole(2);signal*dir_dipole(3)];
end
M = GainMat*Q;
disp_eeg(M,max(abs(M(:))),256,Elec_names,'EEG signal');

%% repeat part e
clc;

peaks = cell(2, 21);

for i=1:21
    [pks,locs] = findpeaks(M(i,:),MinPeakDistance=10,SortStr="descend");
    locs(pks < max(pks/2)) = [];
    pks(pks < max(pks/2)) = [];
    peaks{1,i} = pks;
    peaks{2,i} = locs;
end

figure()
for i=1:21
    subplot( 6 ,4 ,i);
    plot( M(i,:));
    hold on;
    plot(peaks{2,i},peaks{1,i},"o",'LineWidth',1);
    xlabel('smaple','Interpreter','latex');
    ylabel('Amplitude','Interpreter','latex');
    title(Elec_names(1,i),'Fontsize',14,'Interpreter','latex')
end

mean_peaks = zeros(21,1);
for i=1:21
    window = (cell2mat(peaks(2,:))-3):(cell2mat(peaks(2,:))+3);
    mean_peaks(i) = mean(abs(M(i,window)));
end



normalized_mean_peaks = mean_peaks / max(mean_peaks);
for i=1:Elec_numbers
    Elec_positions(:,i)=ElecPos{1, i}.XYZ;
    Elec_names(1,i)=ElecPos{1, i}.Name;
end
colors = parula(100);
tmp = sort(normalized_mean_peaks);
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'.');
hold on
for i=1:21
    color_idx = ceil(normalized_mean_peaks(i)*100);
    scatter3(ModelParams.R(3)*Elec_positions(1,i),ModelParams.R(3)*Elec_positions(2,i),ModelParams.R(3)*Elec_positions(3,i),'MarkerEdgeColor','k','MarkerFaceColor',colors((color_idx),:));
    hold on
end

text(ModelParams.R(3)*Elec_positions(1,:)+0.5,ModelParams.R(3)*Elec_positions(2,:),ModelParams.R(3)*Elec_positions(3,:),Elec_names(1,:),'FontSize',14);
xlabel('x','Fontsize',14,'Interpreter','latex')
ylabel('y','Fontsize',14,'Interpreter','latex')
zlabel('z','Fontsize',14,'Interpreter','latex')
title('Dipoles and Electrodes Locations','FontSize',14,'Interpreter','latex');
colorbar;
%% repeat part f
clc;
figure();
Display_Potential_3D(ModelParams.R(3),mean_peaks);
%% repeat part g - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
M = mean_peaks;
G = GainMat;
alpha=0.1;
Q_MNE=G.'*inv(G*G.'+alpha*eye(21))*M;

% WMNE
Omega = zeros(1317,1317);
for i = 1:1317
    tmp = 0;
    for j = 1:21
        tmp = tmp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    Omega(i,i) = sqrt(tmp);
end
W = kron(Omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; 

% LORETA
A = zeros(1317,1317);
d = 1;
P=1317;
for i = 1:P
    for j = 1:P
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end
A_tmp = A*ones(P);
A0 = inv(A_tmp.*eye(size(A_tmp)))*A;
B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*P,3*P));
W = kron(Omega, eye(3)) * B'*B*kron(Omega, eye(3));
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M;

% sLORETA
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_sLORETA = zeros(3*P, 1);
for i = 1:P
    Q_sLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%% Part m - MNE & WMNE & LORETA & sLORETA
clc;
% MNE
moment_MNE=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    moment_MNE(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
% WMNE
moment_WMNE=zeros(1,length(Q_WMNE)/3);
for i = 1:length(Q_WMNE)/3
    moment_WMNE(1,i) = sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
% LORETA
moment_LORETA=zeros(1,length(Q_LORETA)/3);
for i=1:length(Q_LORETA)/3
    moment_LORETA(1,i)=sqrt(sum(Q_LORETA((i-1)*3+1:i*3).^2,'all'));
end
% sLORETA
moment_sLORETA=zeros(1,length(Q_sLORETA)/3);
for i=1:length(Q_sLORETA)/3
    moment_sLORETA(1,i)=sqrt(sum(Q_sLORETA((i-1)*3+1:i*3).^2,'all'));
end
%% Part n - MNE & WMNE & LORETA & sLORETA
clc;
ground_truth = zeros(length(LocMat),1);
ground_truth(patch_slct) = 1;

figure();
plotroc(ground_truth',moment_MNE);
[tpr_MNE,fpr_MNE,thresholds_MNE] = roc(ground_truth',moment_MNE);
title('ROC of MNE method','interpreter','latex','Fontsize',14);

figure;
plotroc(ground_truth',moment_WMNE);
[tpr_WMNE,fpr_WMNE,thresholds_WMNE] = roc(ground_truth',moment_WMNE);
title('ROC of WMNE method','interpreter','latex','Fontsize',14);

figure;
plotroc(ground_truth',moment_LORETA);
[tpr_LORETA,fpr_LORETA,thresholds_LORETA] = roc(ground_truth',moment_LORETA);
title('ROC of LORETA method','interpreter','latex','Fontsize',14);

figure;
plotroc(ground_truth',moment_sLORETA);
[tpr_sLORETA,fpr_sLORETA,thresholds_sLORETA] = roc(ground_truth',moment_sLORETA);
title('ROC of sLORETA method','interpreter','latex','Fontsize',14);

figure;
subplot(4,2,1);
plot(thresholds_MNE,tpr_MNE,'Linewidth',1);
title('TPR of MNE method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('TPR','interpreter','latex');
grid on;
subplot(4,2,3);
plot(thresholds_MNE,fpr_MNE,'Linewidth',1);
title('FPR of MNE method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('FPR','interpreter','latex');
grid on;
subplot(4,2,2);
plot(thresholds_WMNE,tpr_WMNE,'Linewidth',1);
title('TPR of WMNE method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('TPR','interpreter','latex');
grid on;
subplot(4,2,4);
plot(thresholds_WMNE,fpr_WMNE,'Linewidth',1);
title('FPR of WMNE method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('FPR','interpreter','latex');
grid on;
subplot(4,2,5);
plot(thresholds_LORETA,tpr_LORETA,'Linewidth',1);
title('TPR of LORETA method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('TPR','interpreter','latex');
grid on;
subplot(4,2,7);
plot(thresholds_LORETA,fpr_LORETA,'Linewidth',1);
title('FPR of LORETA method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('FPR','interpreter','latex');
grid on;
subplot(4,2,6);
plot(thresholds_sLORETA,tpr_sLORETA,'Linewidth',1);
title('TPR of sLORETA method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('TPR','interpreter','latex');
grid on;
subplot(4,2,8);
plot(thresholds_sLORETA,fpr_sLORETA,'Linewidth',1);
title('FPR of sLORETA method','interpreter','latex','Fontsize',14);
xlabel('Thresholds','interpreter','latex');
ylabel('FPR','interpreter','latex');
grid on;







