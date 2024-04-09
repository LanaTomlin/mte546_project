import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import savemat
import math

if __name__ == '__main__':

    parser = argparse.ArgumentParser()    
    parser.add_argument('filepath')
    arg = parser.parse_args()

    temp_T501 = [[]];
    temp_T502 = [[]];
    temp_T504 = [[]];
    voltage = [[]];
    current = [[]];

    time_T501 = [[]];
    time_T502 = [[]];
    time_T504 = [[]];
    time_voltage = [[]];
    time_current = [[]];

    status = ['ACTUATOR_OFF'];
    status_num = 0;
    BREAK_OUTPUT = ['ACTUATOR_ON', 'ACTUATOR_OFF']
    sensor_ids = {}
    ##sorting data into various piles
    with open(arg.filepath, 'r', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=':', quotechar='|') 

        for row in spamreader:
            #print(row)
            if(row[0].find('Health') != -1):
                continue
            if(row[2].find('ACTUATOR_STATUS') != -1):
                new_status = row[6].strip().split(',')[0];
                if(new_status != status[-1]):
                    status_num = status_num + 1;
                    status.append(new_status);
                    temp_T501.append([])
                    temp_T502.append([])
                    temp_T504.append([])
                    voltage.append([])
                    current.append([])
                    time_T501.append([])
                    time_T502.append([])
                    time_T504.append([])
                    time_voltage.append([])
                    time_current.append([])
                #print(len(temp_T501))
                #print(len(temp_T501[0]))
            if(row[5].find('T501') != -1):
                #print(row)
                #print(len(temp_T501))
                #print(len(temp_T501[0]))
                temp_T501[status_num].append(float(row[6].split('}')[0].strip()));
                time_T501[status_num].append(float(row[4].split(',')[0].strip()));
            if(row[5].find('T502') != -1):
                temp_T502[status_num].append(float(row[6].split('}')[0].strip()));
                time_T502[status_num].append(float(row[4].split(',')[0].strip()));
            if(row[5].find('T504') != -1):
                temp_T504[status_num].append(float(row[6].split('}')[0].strip()));
                time_T504[status_num].append(float(row[4].split(',')[0].strip()));
            if(row[5].find('GOG_ISENSE') != -1):
                current[status_num].append(float(row[6].split('}')[0])/(50.0*0.005));
                time_current[status_num].append(float(row[4].split(',')[0].strip()));
            if(row[5].find('ISENSE_24V') != -1):
                voltage[status_num].append(float(row[6].split('}')[0])*(64.9 + 1)/(8.2));
                time_voltage[status_num].append(float(row[4].split(',')[0].strip()));


    def print_time_val(time, val):
        f.write(str(time))
        f.write(',')
        f.write(str(val))
        f.write(',')
    status_num = 0;
    def small_len(status_num):
        return min(len(time_current[status_num]), len(time_T501[status_num]), len(time_T502[status_num]), len(time_T504[status_num]), len(time_voltage[status_num]))
    #putting piled data into labelled csv files seperated by status
    for state_found in status:
        if(small_len(status_num) != 0):
            with open(f'{arg.filepath.split(".txt")[0]}_{status[status_num]}_{status_num}.csv', 'w') as f:
                #print(small_len(status_num))    
                for i in range(small_len(status_num)):
                    print_time_val(time_voltage[status_num][i], voltage[status_num][i])
                    print_time_val(time_current[status_num][i], current[status_num][i])
                    print_time_val(time_T501[status_num][i], temp_T501[status_num][i])
                    print_time_val(time_T502[status_num][i], temp_T502[status_num][i])
                    print_time_val(time_T504[status_num][i], temp_T504[status_num][i])
        status_num = status_num + 1

    print(status)
    #find all indexes where actuator is on, make a plot for each of these
    i = []
    for ii in range(len(status)):
        if(status[ii] == "'ACTUATOR_ON'" or status[ii] == "'ACTUATOR_OFF'"):
            i.append(ii)
    fig, ax = plt.subplots(2, 3);
    for ii in i:
        #print(ii)
        #fig, ax = plt.subplots(2, 3);
        fig.suptitle(f'Plot for ACTUATOR_ON {ii}')
        axs = ax[0,0];
        axs.plot((time_T501[ii]), (temp_T501[ii]))
        axs.set_title('Temp T501 (degC) vs Time')
        axs = ax[0,1];
        axs.plot((time_T502[ii]), (temp_T502[ii]))
        axs.set_title('Temp T502 (degC)vs Time')
        axs = ax[0,2];
        axs.plot((time_T504[ii]), (temp_T504[ii]))
        axs.set_title('Temp T504 (degC) vs Time')
        axs = ax[1,0];
        axs.plot((time_voltage[ii]), (voltage[ii]))
        axs.set_title('Voltage Nicrome vs Time')
        axs = ax[1,1];
        axs.plot((time_current[ii]), (current[ii]))
        axs.set_title('Nicrome current (A_) vs Time')
    plt.savefig(f'{arg.filepath.split(".txt")[0]}_ACTUATOR_ALL.png')

    ##doing some shit to get rid of fucked data
    ##remove the values i know to be caused by tc shorting, not by actual measurements
    #and make a list of the old times and temps, so we can use those to plot later
    index_deleted = [[] for _ in range(len(status))]
    def delete_shorted(val, time):
        #these are values caused by tc shorting, gonna wanna delete these
        array_num_to_delete = [63.75, 255.75, 127.75]
        time_T = []
        temp_T = []
        #make big 1D array of all the shit before we delete it
        for status_num in range(len(status)):
            for i in range(len(val[status_num])):
                time_T.append(time[status_num][i])
                temp_T.append(val[status_num][i])
        #delete the tc shorted vals and times
        for status_num in range(len(status)):
            for num_delete in array_num_to_delete:
                while(num_delete in val[status_num]):
                    time[status_num].pop(val[status_num].index(num_delete))
                    #index_deleted[status_num].append(array_new_val[status_num].index(63.75))
                    val[status_num].remove(num_delete)
        return val, time, time_T, temp_T
    temp_T501_new, time_T501_new, time_T501_S, temp_T501_S= delete_shorted(temp_T501, time_T501)
    temp_T502_new, time_T502_new, time_T502_S, temp_T502_S= delete_shorted(temp_T502, time_T502)
    temp_T504_new, time_T504_new, time_T504_S, temp_T504_S= delete_shorted(temp_T504, time_T504)

        #correct the voltage data so we read about 0 when current is low (to indicate that we are disconnected)
    for status_num in range(len(status)):
        voltage_h = 13.139817073170734
        for i in range(len(time_current[status_num])):
            if(current[status_num][i]<= 0.003/(50.0*0.005) and i < len(voltage[status_num])):
                voltage[status_num][i] = abs(voltage[status_num][i])
                voltage[status_num][i] = abs(voltage[status_num][i] - voltage_h)

    def straighten_data(time, val):
        #we put lms coefficients in this array
        opt = [[] for _ in range(len(status))];
        #1D array for all data with bad vals removed
        time_T_n = [];
        temp_T_n = [];
        for status_num in range(len(status)):
            #val[status_num] = np.reshape(val[status_num], (-1)) # = (val[status_num][0])
            #time[status_num] = np.reshape(val[status_num], (-1)) # = (val[status_num][0])
            for i in range(len(val[status_num])):
                time_T_n.append((time[status_num][i]))
                temp_T_n.append(val[status_num][i])
        return time_T_n, temp_T_n
    
    T501_time, T501_temp = straighten_data(time_T501, temp_T501);
    T502_time, T502_temp = straighten_data(time_T502, temp_T502);
    T504_time, T504_temp = straighten_data(time_T504, temp_T504);

    voltage_time, voltage_s = straighten_data(time_voltage, voltage);
    current_time, current_s = straighten_data(time_current, current);
    
    mdic = {'temp_T501': T501_temp, 'time_T501': T501_time, 'time_T501_S': time_T501_S, 'temp_T501_S': temp_T501_S, 'temp_T502': T502_temp, 'time_T502': T502_time, 'time_T502_S': time_T502_S, 'temp_T502_S': temp_T502_S, 'temp_T504': T504_temp, 'time_T504': T504_time, 'time_T504_S': time_T504_S, 'temp_T504_S': temp_T504_S, 'status': status, 'time_current': current_time, 'current': current_s, 'time_voltage': voltage_time, 'voltage': voltage_s}
    savemat(f'{arg.filepath.split(".txt")[0]}_DATA_SAVE.mat', mdic)
    
 