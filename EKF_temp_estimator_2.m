

%import(temp1)

%cross multiply power to figure out how much power we got
%Power = V*i;

%load("good data knock on wood\ONOFF_TEST2HZ\ONOFF_TEST_0_25HZ_HUMAN_INTERFERENCE_LATER\all_messages_2024_04_04-07_10_28_PM_DATA_SAVE")

load("good data knock on wood\ONOFF_TEST0_5HZ\interpolated_data.mat")
T504(1) = 27.4;

resistance = voltage./current;
resit_prev = 7;
for i = 2:length(current)
    if(current(i) > 0.003/(50.0*0.005))
        resistance(i) = voltage(i)/current(i);
        if(resistance(i) > 10)
            resistance(i) = resit_prev;
        end
        resit_prev = resistance(i);
    else
        resistance(i) = resit_prev;
    end
end

figure(2)
subplot(3,2,1)
title('T501 Interpolated')
plot(t, T501);
hold on;
subplot(3,2,2)
title('T502 Interpolated')
plot(t, T502);
subplot(3,2,3)
title('T504 Interpolated')
plot(t, T504);
subplot(3,2,4)
title("Voltage interpolated")
plot(t, voltage)
subplot(3,2,5)
title("Current interpolated")
plot(t, current)
subplot(3,2,6)
title("Resistance interpolated")
plot(t, resistance)


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

ti = 27;
Cp = 0.480*40; % joules/kgdegc
Ar = 53.4751*10^-6; %m^2
m = 8.20*1000*Ar*L; %kg/m^3
K = 13.0; %w/m-K 
temp_est_conv = t; %% just forcing it to be the right length
temp_est_conv(1) = T501(1); %starting point

temp_est_basic = t; %% just forcing it to be the right length
temp_est_basic(1) = T501(1); %starting point

temp_est_1 = t; %% just forcing it to be the right length
temp_est_1(1) = ti; %starting point

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
%sectionT504 = 
K1 = 60;
e = 5.67*10^-8;
emissivity  = 0.5;

% temp_array = [t1; t2; t3; t4; t5; t6];
% coefft1t1 = -(hc*As +K*Ac/L1)*tdiv/(m*Cp);
% coefft1t2 = K1*Ac*tdiv/(L1*m*Cp);
% coefft1t0 = ti*K*Ac*tdiv/(L1*m*cp);
% 
% 
% coefft2t1 = K1*Ac*tdiv/(L1*m*Cp);
% coefft2t2 = -(hc*As +2*K1*Ac/L2)*tdiv/(m*Cp);
% coefft2t3 = K1*Ac*tdiv/(L2*m*Cp);
% coefft2t0 = ti*K*Ac*tdiv/(L2*m*cp);
% 
% coefft3t2 = K1*Ac*tdiv/(L2*m*Cp);
% coefft3t3 = -(hc*As +2*K1*Ac/L2)*tdiv/(m*Cp);
% coefft3t4 = K1*Ac*tdiv/(L3*m*Cp);
% coefft3t0 = ti*K*Ac*tdiv/(L3*m*cp);
% 
% coefft4t3 = K1*Ac*tdiv/(L3*m*Cp);
% coefft4t4 = -(hc*As +2*K1*Ac/L4)*tdiv/(m*Cp);
% coefft4t5 = K1*Ac*tdiv/(L4*m*Cp);
% coefft4t0 = ti*K*Ac*tdiv/(L4*m*cp);
% 
% coefft5t4 = K1*Ac*tdiv/(L5*m*Cp);
% coefft5t5 = -(hc*As +2*K1*Ac/L5)*tdiv/(m*Cp);
% coefft5t6 = K1*Ac*tdiv/(L5*m*Cp);
% coefft5t0 = t*K*Ac*tdiv/(L5*m*cp);
% 
% coefft6t6 = -(hc*As +K*Ac/L6)*tdiv/(m*Cp);
% coefft6t5 = K1*Ac*tdiv/(L6*m*Cp);
% coefft6t0 = ti*K*Ac*tdiv/(L6*m*cp);
% 
% A = [1, 0, 0, 0, 0, 0, 0;
%     coefft1t0, coefft1t1,  coefft1t2, 0, 0, 0 0;
%     coefft1t0, coefft2t1, coefft2t2, coefft2t3, 0, 0, 0;
%     coefft3t0, 0, coefft3t2, coefft3t3, coefft3t4, 0, 0;
%     coefft4t0, 0, 0, coefft4t3, coefft4t4, coefft4t5, 0;
%     coefft5t0, 0, 0, 0, coefft5t4, coefft5t5, coefft6t5;
%     coefft6t0, 0, 0, 0, 0, coefft6t5, coefft6t6;];





for i = 2:length(t)
    %section 1 power terminal
    %L = 5*10*-3; 
    temp_est_1(i) = (current(i)*voltage(i) - (hc*As +K1*Ac/L1)*(temp_est_1(i-1) - ti) - K*Ac*(temp_est_1(i-1) - temp_est_2(i-1))/L1)*tdiv/(Cp*m*500) + temp_est_1(i-1);

    %section 2 T502
    %L = 20*10*-3; 
    temp_est_2(i) = (current(i)*voltage(i) - (hc*As)*(temp_est_2(i-1) - ti) - K*Ac*((temp_est_2(i-1) - temp_est_1(i-1))/L1 + (temp_est_2(i-1) - temp_est_3(i-1))/L2))*tdiv/(Cp*m*20) + temp_est_2(i-1);

    %section 3 mounting clip
    %L = 20*10*-3; 
    temp_est_3(i) = (current(i)*voltage(i) - (hc*As +K1*Ac/L3)*(temp_est_3(i-1) - ti) - K*Ac*((temp_est_3(i-1) - temp_est_2(i-1))/L2 + (temp_est_3(i-1) - temp_est_4(i-1))/L3))*tdiv/(m*Cp*20) + temp_est_3(i-1);


    %section 4 T501
    %L = 20*10*-3; 
    temp_est_4(i) = (current(i)*voltage(i) - (hc*As)*(temp_est_4(i-1) - ti) - K*Ac*((temp_est_4(i-1) - temp_est_3(i-1))/L3 + (temp_est_4(i-1) - temp_est_5(i-1))/L4))*tdiv/(m*Cp) + temp_est_4(i-1);


    %section 5 mounting clip
    %L = 20*10*-3; 
    temp_est_5(i) = (current(i)*voltage(i) - (hc*As +K1*Ac/L5)*(temp_est_5(i-1) - ti) - K*Ac*((temp_est_5(i-1) - temp_est_4(i-1))/L4+(temp_est_5(i-1) - temp_est_6(i-1))/L5))*tdiv/(m*Cp*200) + temp_est_5(i-1);

    %power 2 terminal and T504
    %L = 20*10*-3; 
    temp_est_6(i) = (current(i)*voltage(i) - (hc*As +K1*Ac/L6)*(temp_est_6(i-1) - ti) - K*Ac*(temp_est_6(i-1) - temp_est_5(i-1))/L5)*tdiv/(m*Cp*20) + temp_est_6(i-1);

    %no losses
    %temp_est(i) = current(i)*voltage(i)*tdiv/(m*Cp) + temp_est(i-1);
    temp_est_conv(i) = (current(i)*voltage(i) - hc*Ac*(temp_est_conv(i-1) - ti))*tdiv/(m*Cp) + temp_est_conv(i-1);
    temp_est_basic(i) = (current(i)*voltage(i) - hc*Ac*(temp_est_conv(i-1) - ti) - Ac*e*emissivity*((temp_est_conv(i-1)+273) - ti))*tdiv/(m*Cp) + temp_est_conv(i-1);
end


figure(3)
subplot(3,2,1)
title('section 1 estimated, basic')
plot(t, temp_est_1)
hold on;
subplot(3, 2, 2)
plot(t, temp_est_2)
hold on;
plot(t, T502)
subplot(3,2,3);
plot(t, temp_est_3);
subplot(3,2,4)
plot(t, temp_est_4)
hold on
plot(t, T501);
subplot(3,2,5)
plot(t, temp_est_5);
hold on 
subplot(3,2,6)
title("temp est")
plot(t, temp_est_6);
hold on 
plot(t, T504)


hold on;
%qc + qr + mCp*dTc/dt = qs + (I^2)*R
A = [1, -hc*Ac*tdiv/(m*Cp)+1];
B = ti*hc*Ac*tdiv/(m*Cp);
H = [1, 0];
T = [ti];
Tk = [ti];
Pk = [0];
Kk = [0];
Q = [0];
Kk = [0];
%H = [1];
R = [0];
for i = 2:(length(t))

    Tnew = A*[current(i)*voltage(i); T(i-1)] + B;
    T = [T,Tnew];
    %prediction
    Tk = [Tk,Tnew];
    Pk = [Pk, A*Pk(i-1)*A' + Q];
    Zk = [current(i)*voltage(i); Tk(i)];
    %correction
    Kk = [Kk, (Pk(i)*H')*(H*Pk(i)*H' + R)'];
    Tk(i) = Tk(i) +Kk*(Zk - H*Tk(i))
    %T(i) = T(i) + B;
end

figure(4)
plot(t, temp_est_conv)
hold on
plot(t, temp_est_basic)
plot(t, T)
