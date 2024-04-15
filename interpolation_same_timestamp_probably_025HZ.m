

%import(temp1)

%cross multiply power to figure out how much power we got
%Power = V*i;

load("good data knock on wood\ONOFF_UNKNOWN_PROBABLY_0_25HZ\all_messages_2024_04_04-07_04_52_PM_DATA_SAVE.mat")


t501_time = [0, time_T501{:,1}, time_T501{:,2}];
t502_time = [0, time_T502{:,1}, time_T502{:,2}];
t504_time = [0, time_T504{:,1}, time_T504{:,2}];

t501_temp = [96.25, temp_T501{:,1}, temp_T501{:,2}];
t502_temp = [48, temp_T502{:,1}, temp_T502{:,2}];
t504_temp = [29.75, temp_T504{:,1}, temp_T504{:,2}];

voltage_time = [time_voltage{:,1}, time_voltage{:,2}];
voltage = [voltage{:,1}, voltage{:,2}];

current_time = [time_current{:,1}, time_current{:,2}];
current = [current{:,1}, current{:,2}];

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

% curveFitter
load("good data knock on wood\ONOFF_UNKNOWN_PROBABLY_0_25HZ\interpolated_data.mat")

t = 0:0.1:200;
T501 = [T501_fittedmodel(t(1:600))', T501_60Hfittedmodel(t(601:2001))'];
T502 = T502_fittedmodel(t)';
T504 = T504_fittedmodel(t)';
voltage = voltage_fittedmodel(t)';
current = current_fittedmodel(t)';

figure(1)
subplot(3,2,1)
plot(t501_time, t501_temp); hold on;
plot(t, T501);
title('T501 Interpolated')

subplot(3,2,2)
plot(t502_time, t502_temp); hold on;
plot(t, T502);
title('T502 Interpolated')

subplot(3,2,3)
plot(t504_time, t504_temp); hold on;
plot(t, T504);
title('T504 Interpolated')

subplot(3,2,4)
plot(t, voltage)
title("Voltage interpolated")

subplot(3,2,5)
plot(t, current)
title("Current interpolated")


for i = 1:length(voltage)
    if(voltage(i) > 1.24)
    voltage(i) = voltage(i) - 1.24;
    end
end

% curveFitter

% save("good data knock on wood\ONOFF_UNKNOWN_PROBABLY_0_25HZ\interpolated_data", "T501_fittedmodel", "T501_60Hfittedmodel", "T502_fittedmodel", "T504_fittedmodel", "current_fittedmodel", "voltage_fittedmodel", "T501_goodness", 'T501', "T502", 'T504', "current", "voltage", 't')