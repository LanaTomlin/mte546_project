

%import(temp1)

%cross multiply power to figure out how much power we got
%Power = V*i;

load("good data knock on wood\ONOFF_TEST_0_25HZ_HUMAN_INTERFERENCE_LATER\all_messages_2024_04_04-07_10_28_PM_DATA_SAVE.mat")

t501_time = [0, time_T501{:,1}, time_T501{:,2}, time_T501{:,3}, time_T501{:,4}, time_T501{:,5}];
t502_time = [0, time_T502{:,1}, time_T502{:,2}, time_T502{:,3}, time_T502{:,4}, time_T502{:,5}];
t504_time = [0, time_T504{:,1}, time_T504{:,2}, time_T504{:,3}, time_T504{:,4}, time_T504{:,5}];

t501_temp = [96.25, temp_T501{:,1}, temp_T501{:,2}, temp_T501{:,3}, temp_T501{:,4}, temp_T501{:,5}];
t502_temp = [48, temp_T502{:,1}, temp_T502{:,2}, temp_T502{:,3}, temp_T502{:,4}, temp_T502{:,5}];
t504_temp = [29.75, temp_T504{:,1}, temp_T504{:,2}, temp_T504{:,3}, temp_T504{:,4}, temp_T504{:,5}];

voltage_time = [time_voltage{:,1}, time_voltage{:,2}, time_voltage{:,3}, time_voltage{:,4}, time_voltage{:,5}];
voltage = [voltage{:,1}, voltage{:,2}, voltage{:,3}, voltage{:,4}, voltage{:,5}];

current_time = [time_current{:,1}, time_current{:,2}, time_current{:,3}, time_current{:,4}, time_current{:,5}];
current = [current{:,1}, current{:,2}, current{:,3}, current{:,4}, current{:,5}];

% for i = 2:185835
%     if(voltage_time(i-1) >= time_voltage(i));
%         time_voltage(i) = [];
%         voltage(i) = [];
%     end
% end

[voltage_time, ia, ic] = unique(voltage_time);
voltage = voltage(ia);



[current_time, ia, ic] = unique(current_time);
current = current(ia);

figure(1)
plot(current_time, current)

figure(2)
plot(voltage_time, voltage)

curveFitter

t = 0:0.1:600;
T501 = [T501_350_Lfittedmodel(t(1:3500))' T501_350_Hfittedmodel(t(3501:6001))'];
T502 = T502_fittedmodel(t)';
T504 = T504_fittedmodel(t)';
voltage = voltage_fittedmodel(t)';
current = current_fittedmodel(t)';

figure(1)
subplot(3,2,1)
title('T501 Interpolated')
plot(t501_time, t501_temp)
hold on;
plot(t, T501);
subplot(3,2,2)
title('T502 Interpolated')
plot(t502_time, t502_temp)
hold on;
plot(t, T502);
subplot(3,2,3)
title('T504 Interpolated')
plot(t504_time, t504_temp)
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

curveFitter

save("good data knock on wood\ONOFF_TEST_0_25HZ_HUMAN_INTERFERENCE_LATER\interpolated_data", "T501_350_Hfittedmodel", "T501_350_Lfittedmodel", "T502_fittedmodel", "T504_fittedmodel", "current_fittedmodel", "voltage_fittedmodel", 'T501', "T502", 'T504', "current", "voltage", 't')