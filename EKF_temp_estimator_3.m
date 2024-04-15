close all;
clc;
clear all;

% Data will be re-tuned for this dataset, different file will be used for
% different data file
load("good data knock on wood\ONOFF_UNKNOWN_PROBABLY_0_25HZ\interpolated_data.mat")

T504(1) = 27.4;

resistance = voltage./current;
resit_prev = 7;
for i = 2:length(current)
    if(current(i) > 2)
        resistance(i) = 11.84/current(i);
        %if(resistance(i) > 10)
         %   resistance(i) = resit_prev;
        %end
        resit_prev = resistance(i);
    else
        resistance(i) = resit_prev;
    end
end

figure(2)
subplot(3,2,1)
ax = plot(t, T501);
title('T501 Interpolated')

subplot(3,2,2)
ax = plot(t, T502);
title('T502 Interpolated')

subplot(3,2,3)
ax = plot(t, T504);
title('T504 Interpolated')

subplot(3,2,4)
ax = plot(t, voltage);
title("Voltage interpolated")

subplot(3,2,5)
ax = plot(t, current);
title("Current interpolated")

subplot(3,2,6)
ax = semilogy(t, resistance);
title("Resistance interpolated")


L1 = 5*10^-3;
L2 = 60*10^-3;
L3 = 5*10^-3;
L4 = 50*10^-3;
L5 = 5*10^-3;
L6 = 500*10^-3;

%actual energy
%Q = voltage*current;
%conduction heat transfer
L = L1 +L2  + L3 + L4 + L5 + L6; %0.3302; %legnth in m
%A = 1.3*10^(-7); %https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.bartleby.com/questions-and-answers/a-nichrome-wire-has-a-cross-sectional-area-1.3-x-10-7m2.-resistivity-of-nichrome-1.1-x-10-8wmand-res/e3d695f0-cbbd-43cb-adea-8670792888cb&ved=2ahUKEwiV5P75yrWFAxXMMDQIHVc3AjYQFnoECBcQAw&usg=AOvVaw1lSSXpHmuGMqKW2MptLixl
%Q_cond = k*A*t/L;
%convection heat transfer
%Q_conv = h*

T_initial = 27;
Cp = 0.480*40; % joules/kgdegc
Ar = 53.4751*10^-6; %m^2
m = 8.20*1000*Ar*L; %kg/m^3
K = 13.0; %w/m-K 
temp_est_conv = t; %% just forcing it to be the right length
temp_est_conv(1) = T501(1); %starting point

temp_est_basic = t; %% just forcing it to be the right length
temp_est_basic(1) = T501(1); %starting point

temp_est_1 = t; %% just forcing it to be the right length
temp_est_1(1) = T_initial; %starting point

temp_est_2 = t; %% just forcing it to be the right length
temp_est_2(1) = T502(1); %starting point

temp_est_3 = t; %% just forcing it to be the right length
temp_est_3(1) = T502(1); %starting point

temp_est_4 = t; %% just forcing it to be the right length
temp_est_4(1) = T501(1); %starting point


temp_est_5 = t; %% just forcing it to be the right length
temp_est_5(1) = T504(1); %starting point

temp_est_6 = t; %% just forcing it to be the right length
temp_est_6(1) = T504(1); %starting point

tdiv = 0.1;
hc = 50; % W/m^2k convective heat transferfor air

Ac = 2*8.2515*pi*10^-6;
As = 2*8.2515*pi*10^-6;

%sectionT504 = 
K1 = 60;
e = 5.67*10^-8;
emissivity  = 0.5;


for i = 2:length(t)
    pwr_i = current(i)*voltage(i);
    %section 1 power terminal
    %L = 5*10*-3; 
    temp_est_1(i) = (pwr_i - (hc*As +K1*Ac/L1)*(temp_est_1(i-1) - T_initial) - ...
        K*Ac*(temp_est_1(i-1) - temp_est_2(i-1))/L1)*tdiv/(Cp*m*500) + temp_est_1(i-1);

    %section 2 T502
    %L = 20*10*-3; 
    temp_est_2(i) = (pwr_i - (hc*As)*(temp_est_2(i-1) - T_initial) - ...
        K*Ac*((temp_est_2(i-1) - temp_est_1(i-1))/L1 + (temp_est_2(i-1) - temp_est_3(i-1))/L2))*tdiv/(Cp*m*20) + temp_est_2(i-1);

    %section 3 mounting clip
    %L = 20*10*-3; 
    temp_est_3(i) = (pwr_i - (hc*As +K1*Ac/L3)*(temp_est_3(i-1) - T_initial) ...
        - K*Ac*((temp_est_3(i-1) - temp_est_2(i-1))/L2 + (temp_est_3(i-1) - temp_est_4(i-1))/L3))*tdiv/(m*Cp*20) + temp_est_3(i-1);

    %section 4 T501
    %L = 20*10*-3; 
    temp_est_4(i) = (pwr_i - (hc*As)*(temp_est_4(i-1) - T_initial) ...
        - K*Ac*((temp_est_4(i-1) - temp_est_3(i-1))/L3 + (temp_est_4(i-1) - temp_est_5(i-1))/L4))*tdiv/(m*Cp) + temp_est_4(i-1);


    %section 5 mounting clip
    %L = 20*10*-3; 
    temp_est_5(i) = (pwr_i - (hc*As +K1*Ac/L5)*(temp_est_5(i-1) - T_initial) - ...
        K*Ac*((temp_est_5(i-1) - temp_est_4(i-1))/L4+(temp_est_5(i-1) - temp_est_6(i-1))/L5))*tdiv/(m*Cp*200) + temp_est_5(i-1);

    %power 2 terminal and T504
    %L = 20*10*-3; 
    temp_est_6(i) = (pwr_i - (hc*As +K1*Ac/L6)*(temp_est_6(i-1) - T_initial) - ...
        K*Ac*(temp_est_6(i-1) - temp_est_5(i-1))/L5)*tdiv/(m*Cp*20) + temp_est_6(i-1);

    %no losses
    %temp_est(i) = current(i)*voltage(i)*tdiv/(m*Cp) + temp_est(i-1);
    temp_est_conv(i) = (current(i)*voltage(i) - ...
        hc*Ac*(temp_est_conv(i-1) - T_initial))*tdiv/(m*Cp) + temp_est_conv(i-1);
    temp_est_basic(i) = (current(i)*voltage(i) - ...
        hc*Ac*(temp_est_conv(i-1) - T_initial) - Ac*e*emissivity*((temp_est_conv(i-1)+273) - T_initial))*tdiv/(m*Cp) + temp_est_conv(i-1);
end


figure(3)

subplot(3,2,1)
plot(t, temp_est_1)
title('section 1 estimated, basic')

subplot(3,2,2)
plot(t, temp_est_2)
hold on;
plot(t, T502)
title('section 2 estimated + T502 overlaid')

subplot(3,2,3);
plot(t, temp_est_3);
title('section 3 estimated, basic')

subplot(3,2,4)
plot(t, temp_est_4)
hold on
plot(t, T501);
title('section 4 estimated + T501 overlaid')

subplot(3,2,5)
plot(t, temp_est_5);
title('section 5 estimated, basic')

subplot(3,2,6)
plot(t, temp_est_6);
hold on;
plot(t, T504)
title('section 6 estimated + T504 overlaid')


dc = 0.6;
Vo = 11.84;

sigmat_meas = (15*10^-3)^2;
sigmat_adc = (15*10^-3)^2;

sig_voltage = (15*10^-3);

Cp = Cp*1000/40;
Ro = 4;
alpha = 0.00017; %temperature coefficeint of resistance nichrome wire
sigmav = ((sigmat_meas^(1/2))*Ro*alpha)^2;

A = [1 - tdiv*hc*As      , 0;
     -tdiv*alpha*Ro*hc*As, 1];

% B = [(Vo^2)*Ro/(m*Cp);  -(hc*As*ti*tdiv)/(m*Cp); 0];
B = [((Vo^2)/(Ro*m*Cp))*tdiv; 
     ((Vo^2)/(Ro*m*Cp))*tdiv*Ro*alpha;];


H = [0, Vo^(-2)];
% X = [T, R]
X_predicted_array = [[0; 0;]];
X_corrected_array = [[0; Ro]];

% Pk|k-1 will not be stored explicitly
P_kk_array = cell([1 length(t)]);
Kk = cell([1 length(t)]);

% R and Q matrices

% 2 state variables, Q is 2x2
sig_T = tdiv*(sig_voltage^2 / Ro)/(m*Cp);
sig_R = tdiv*alpha*Ro*(sig_voltage^2 / Ro)/(m*Cp);

Q = [sig_T^2    , sig_T*sig_R;
     sig_T*sig_R, sig_R^2    ;];

% one measurement, R is 1x1
R = [((sigmat_adc^2)/Ro)^2]*0.4


P_kk_array{1} = Q;
y_hat_k = cell([1 length(t)]);

BU = [];

for i = 2:(length(t))

    % classification of whether relay is on or not
    power_on = 1;
    if(current(i-1) < 2)
        power_on = 0;
    end

    %prediction
    X_k_given_km1 = A*X_corrected_array(:,i-1) + B*[power_on];
    BU = [BU, B*[power_on]];
    P_k_given_km1 = A*P_kk_array{i - 1}*A' + Q;

    X_predicted_array = [X_predicted_array, X_k_given_km1];
    
    
    %correction
    zk = (current(i)*Vo)^(-1);
    y_hat_k{i} = zk - H*X_k_given_km1;

    S_k = H*P_k_given_km1*H' + R;
    Kk{i} = (P_k_given_km1*H')/(S_k);
    
    X_k_given_k = X_k_given_km1 + Kk{i}*y_hat_k{i};
    P_kk_array{i} = (eye(2) - Kk{i}*H)*P_k_given_km1;

    
    X_corrected_array = [X_corrected_array, X_k_given_k];
    
    %T(i) = T(i) + B;
end

yhat_mat = cell2mat(y_hat_k);
Kk_mat = cell2mat(Kk);

figure(4)

hold on
plot(t, current*Vo)
% plot(t, X_corrected_array(1, :));
% plot(t, X_corrected_array(1, :));
% plot(t, X_predicted_array(1, :));
plot(t, [0 yhat_mat]);
plot(t, [0 Kk_mat(1, :)]);
% plot(t, X_corrected_array(2, :)*Vo^2);
plot(t, [0 BU(1, :)]);
legend("zk", "T corrected", "T predicted", "yhat", "Kk", "HXk", "BU")
% legend("zk", "Kk")

figure(5)
hold on
plot(t, temp_est_conv)
plot(t, X_corrected_array(1, :) + T_initial);
plot(t, X_predicted_array(1, :) + T_initial);
legend("temp_est_conv", "T corrected", "T predicted")

% legend("temp_est_conv", "T corrected")


autoArrangeFigures(2,2,3)

% disp('')


%state 
% X =  T1, T2, T3, T4, T5, T6, R


K = 13.0; %w/m-K 
hc = 50;
Cp = Cp*10;

% Diagonal elements
A11 = 1 - tdiv*hc*As*(L1/L);
A22 = 1 - tdiv*hc*As*(L2/L);
A33 = 1 - tdiv*hc*As*(L3/L);
A44 = 1 - tdiv*hc*As*(L4/L);
A55 = 1 - tdiv*hc*As*(L5/L);
A66 = 1 - tdiv*hc*As*(L6/L);
A77 = 1;

% Resistance calculation elements are based on an averga of all
% temperatures

A71 = -tdiv*hc*As*(L1/L)*alpha*Ro;
A72 = -tdiv*hc*As*(L2/L)*alpha*Ro;
A73 = -tdiv*hc*As*(L3/L)*alpha*Ro;
A74 = -tdiv*hc*As*(L4/L)*alpha*Ro;
A75 = -tdiv*hc*As*(L5/L)*alpha*Ro;
A76 = -tdiv*hc*As*(L6/L)*alpha*Ro;

A = ...
[
A11,0  ,0  ,0  ,0  ,0  ,0  ;
0  ,A22,0  ,0  ,0  ,0  ,0  ;
0  ,0  ,A33,0  ,0  ,0  ,0  ;
0  ,0  ,0  ,A44,0  ,0  ,0  ;
0  ,0  ,0  ,0  ,A55,0  ,0  ;
0  ,0  ,0  ,0  ,0  ,A66,0  ;
A71,A72,A73,A74,A75,A76,A77;
]

B_base = ((Vo^2)*Ro/(m*Cp))*tdiv
B1 = B_base*(L1/L);
B2 = B_base*(L2/L);
B3 = B_base*(L3/L);
B4 = B_base*(L4/L);
B5 = B_base*(L5/L);
B6 = B_base*(L6/L);
B7 = ((Vo^2)*Ro/(m*Cp))*tdiv*Ro*alpha;
B = [B1 B2 B3 B4 B5 B6 B7]';

H = ...
[0 0 0 0 0 0 Vo^(-2)]

sig_T_base = tdiv*(sig_voltage^2 * Ro)/(m*Cp);

Q11 = sig_T_base*(L1/L);
Q22 = sig_T_base*(L2/L);
Q33 = sig_T_base*(L3/L);
Q44 = sig_T_base*(L4/L);
Q55 = sig_T_base*(L5/L);
Q66 = sig_T_base*(L6/L);
Q77 = tdiv*alpha*Ro*(sig_voltage^2 * Ro)/(m*Cp);

% Noise will be worse but a proper matrix would be too difficult to
% compute
Q = diag([Q11 Q22 Q33 Q44 Q55 Q66 Q77]);

R = [(Ro*sigmat_adc^2)^2]; % one measruement TF 1x1 array

% Process Arrays

X_kk = cell([1 length(t)]);

X_kkm1 = cell([1 length(t)]);

P_kkm1 = cell([1 length(t)]);
P_kk = cell([1 length(t)]);

S_k = cell([1 length(t)]);
K_k = cell([1 length(t)]);
yhat_k = cell([1 length(t)]);

% Initial conditions
X_kk{1} = [0; 0; 0; 0; 0; 0; Ro];
P_kk{1} = Q;


for i = 2:(length(t))

    % classification of whether relay is on or not
    power_on = 1;
    if(current(i-1) < 2)
        power_on = 0;
    end

    %prediction
    X_k_given_km1 = A*X_kk{i - 1} + B.*power_on;
    P_k_given_km1 = A*P_kk{i - 1}*A' + Q;

    X_kkm1{i} = X_k_given_km1;
    P_kkm1{i} = P_k_given_km1;
    
    %correction
    zk = (current(i)*Vo)^(-1);
    yhat_k{i} = zk - H*X_k_given_km1;

    S_k = H*P_k_given_km1*H' + R;
    K_k{i} = (P_k_given_km1*H')/(S_k);
    
    X_k_given_k = X_k_given_km1 + K_k{i}*y_hat_k{i};
    X_kk{i} = X_k_given_k;
    P_kk{i} = (eye(7) - K_k{i}*H)*P_k_given_km1;
end

X_mat = cell2mat(X_kk);

figure(6);
subplot(3,2,1);
title('Temperature estimate at Section 1');
plot(t, X_mat(2,:));


subplot(3, 2, 2); hold on;
plot(t, X_mat(3,:))
plot(t, T502)
title('T502 at Section 2, with filtered estimate overlaid');
legend('Filtered Estimate', "True Temperature at Node");

e_node_502 = abs(X_mat(3,:) - T502);
RMS_e_node_502 = sqrt((1/(length(t))*sum(e_node_502)))

subplot(3,2,3); hold on;

plot(t, X_mat(4,:))

subplot(3,2,4); hold on;
plot(t, X_mat(5,:))
plot(t, T501);
title('T501 at section 4, with filtered estimate overlaid');
legend('Filtered Estimate', "True Temperature at Node");

e_node_501 = abs(X_mat(5,:) - T501);
RMS_e_node_501 = sqrt((1/(length(t))*sum(e_node_501)))


subplot(3,2,5);
plot(t, X_mat(6,:))

subplot(3,2,6); hold on; 
% title("temp est")

plot(t, X_mat(7,:))
plot(t, T504)

e_node_504 = abs(X_mat(7,:) - T504);
RMS_e_node_504 = sqrt((1/(length(t))*sum(e_node_504)))

title('T504 at section 6 with filtered estimate overlaid');
legend('Filtered Estimate', "True Temperature at Node");



