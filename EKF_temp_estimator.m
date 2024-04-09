

%import(temp1)

%cross multiply power to figure out how much power we got
%Power = V*i;

%load("good data knock on wood\ONOFF_TEST2HZ\ONOFF_TEST_0_25HZ_HUMAN_INTERFERENCE_LATER\all_messages_2024_04_04-07_10_28_PM_DATA_SAVE")

load("good data knock on wood\ONOFF_TEST0_5HZ\interpolated_data.mat")


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

%actual energy
Q = voltage*current;
%conduction heat transfer
L = 0.3302;
A = 1.3*10^(-7); %https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.bartleby.com/questions-and-answers/a-nichrome-wire-has-a-cross-sectional-area-1.3-x-10-7m2.-resistivity-of-nichrome-1.1-x-10-8wmand-res/e3d695f0-cbbd-43cb-adea-8670792888cb&ved=2ahUKEwiV5P75yrWFAxXMMDQIHVc3AjYQFnoECBcQAw&usg=AOvVaw1lSSXpHmuGMqKW2MptLixl
Q_cond = k*A*t/L
%convection heat transfer
%Q_conv = h*

%qc + qr + mCp*dTc/dt = qs + (I^2)*R