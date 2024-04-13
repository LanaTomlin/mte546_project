

%import(temp1)

%cross multiply power to figure out how much power we got
%Power = V*i;

load("good data knock on wood\ONOFF_TEST2HZ\all_messages_2024_04_04-06_33_12_PM_DATA_SAVE.mat")

% t501_time = [0, time_T501{:,1}, time_T501{:,2}];
% t502_time = [0, time_T502{:,1}, time_T502{:,2}];
% t504_time = [0, time_T504{:,1}, time_T504{:,2}];
% 
% t501_temp = [96.25, temp_T501{:,1}, temp_T501{:,2}];
% t502_temp = [48, temp_T502{:,1}, temp_T502{:,2}];
% t504_temp = [29.75, temp_T504{:,1}, temp_T504{:,2}];
% 
% voltage_time = [time_voltage{:,1}, time_voltage{:,2}];
% voltage = [voltage{:,1}, voltage{:,2}];
% 
% current_time = [time_current{:,1}, time_current{:,2}];
% current = [current{:,1}, current{:,2}];

% for i = 2:185835
%     if(voltage_time(i-1) >= time_voltage(i));
%         time_voltage(i) = [];
%         voltage(i) = [];
%     end
% end

[time_voltage, ia, ic] = unique(time_voltage);
voltage = voltage(ia);



[time_current, ia, ic] = unique(time_current);
current = current(ia);

figure(1)
plot(time_current, current)

figure(2)
plot(time_voltage, voltage)

curveFitter

t = 0:0.1:600;
T501 = [T501_4000fittedmodel(t(1:4000))', T501_650fittedmodel(t(4001:6001))'];
T502 = [T502_450fittedmodel(t(1:4500))', T502_650fittedmodel(t(4501:6001))'];
T504 = T504_fittedmodel(t)';
voltage = voltage_fittedmodel(t)';
current = current_fittedmodel(t)';

figure(1)
subplot(3,2,1)
title('T501 Interpolated')
plot(time_T501, temp_T501)
hold on;
plot(t, T501);
subplot(3,2,2)
title('T502 Interpolated')
plot(time_T502, temp_T502)
hold on;
plot(t, T502);
subplot(3,2,3)
title('T504 Interpolated')
plot(time_T504, temp_T504)
hold on;
plot(t, T504);
subplot(3,2,4)
title("Voltage interpolated")
plot(t, voltage)
subplot(3,2,5)
title("Current interpolated")
plot(t, current)

for i = 1:length(voltage)
    if(voltage(i) > 1.24)
    voltage(i) = voltage(i) - 1.24;
    end
end


save("good data knock on wood\ONOFF_TEST2HZ\interpolated_data", "T501_4000fittedmodel", "T501_650fittedmodel", "T502_450fittedmodel", "T502_650fittedmodel", "T504_fittedmodel", "current_fittedmodel", "voltage_fittedmodel", 'T501', "T502", 'T504', "current", "voltage", 't')