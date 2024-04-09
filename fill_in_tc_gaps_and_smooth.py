
import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import loadmat
import math

if __name__ == '__main__':

    parser = argparse.ArgumentParser()    
    parser.add_argument('filepath')
    arg = parser.parse_args()

    t_div = 0.1;

    override_6 = []
    override_1 = []
    override_4 = []
    override_4_state = []
    time_stamp_good = []
    if(arg.filepath.find('ONOFF_TEST0_5HZ') != -1):
        override_6 = ['T501']
        override_6_state = [0, 1, 2]
        override_1 = []
        time_stamp_good = [0, 180]
    if(arg.filepath.find('HUMAN_INTERFERENCE') != -1):
        override_6 = []
        override_1_state = []
        override_6_state = []
        override_4 = ['T501', 'T502']
        override_4_state = [2, 3, 4]
        override_1 = []
        time_stamp_good = [0, 500]
    if(arg.filepath.find('ONOFF_TEST2HZ') != -1):
        override_6 = ['T501']
        override_6_state = [3, 4]
        override_1 = []
        override_4 = ['T502']
        override_4_state = [0, 1]
        time_stamp_good = [0, 600]

    #mdic = {'temp_T501': mat['temp_T501'], 'time_T501': mat['time_T501'], 'mat['time_T501_S']': mat['time_T501_S'], 'mat['temp_T501_S']': mat['temp_T501_S'], 'temp_T502': mat['temp_T502'], 'time_T502': mat['temp_T502'], 'mat['time_T502_S']': mat['time_T502_S'], 'mat['temp_T502_S']': mat['temp_T502_S'], 'mat['temp_T504']': mat['mat['temp_T504']'], 'mat['temp_T504']': mat['mat['temp_T504']'], 'mat['mat['temp_T504']_S']': mat['mat['temp_T504']_S'], 'mat['mat['temp_T504']_S']': mat['mat['temp_T504']_S'], 'mat['status']': mat['status'], 'time_mat['current']': time_mat['current'], 'mat['current']': mat['current'], 'mat['time_mat['voltage']']': mat['time_mat['voltage']'], 'mat['voltage']': mat['voltage']}
    mat = loadmat(f'{arg.filepath.split(".txt")[0]}_DATA_SAVE.mat')
    
    mat['time_T501'] = np.reshape(mat['time_T501'], (13, -1))
    mat['temp_T501'] = np.reshape(mat['temp_T501'], (13, -1))
    print(np.asarray(list(mat['time_T501'])).shape)
    print(mat['status'])
    print(mat['time_T501'])


    #function for least mean square fitting
    def func6(x, a, b, c, d, e, f, g):
        #return a * np.power(x, 4) + b*np.power(x, 3) + c*np.power(x, 2) + d*np.power(x, 1) + e
        #return c*np.power(x, 2) + d*np.power(x, 1) + e
        return a * np.power(x, 6) + b*np.power(x, 5) + c*np.power(x, 4) + d*np.power(x, 3) + e*np.power(x, 2) + f*np.power(x, 1) + g

        #function for least mean square fitting
    def func4(x, a, b, c, d, e):
        return a * np.power(x, 4) + b*np.power(x, 3) + c*np.power(x, 2) + d*np.power(x, 1) + e
        #return c*np.power(x, 2) + d*np.power(x, 1) + e
        #return a * np.power(x, 6) + b*np.power(x, 5) + c*np.power(x, 4) + d*np.power(x, 3) + e*np.power(x, 2) + f*np.power(x, 1) + g

        #function for least mean square fitting
    def func2(x, a, b, c):
        #return a * np.power(x, 4) + b*np.power(x, 3) + c*np.power(x, 2) + d*np.power(x, 1) + e
        return a*np.power(x, 2) + b*np.power(x, 1) + c
        #return a * np.power(x, 6) + b*np.power(x, 5) + c*np.power(x, 4) + d*np.power(x, 3) + e*np.power(x, 2) + f*np.power(x, 1) + g

    def func1(x, a, b, c):
        #return a * np.power(x, 4) + b*np.power(x, 3) + c*np.power(x, 2) + d*np.power(x, 1) + e
        return a*np.exp(-b*x) + c
        #return a * np.power(x, 6) + b*np.power(x, 5) + c*np.power(x, 4) + d*np.power(x, 3) + e*np.power(x, 2) + f*np.power(x, 1) + g

    #generate least means squares to estimate the data
    #takes the data with bad vals removed, and makes least means squares estimate
    #from mat['current'] mat['status'] to next non-zero mat['status'] (ie, we didnt remove all the fucking data)
    def calc_lmsq(time, val, id_override):
        #we put lms coefficients in this array
        opt = [[] for _ in range(len(mat['status']))];
        func_order = [func6 for _ in range(len(mat['status']))];
        #1D array for all data with bad vals removed
        time_T_n = [];
        temp_T_n = [];
        for status_num in range(len(mat['status'])):
            #val[status_num] = np.reshape(val[status_num], (-1)) # = (val[status_num][0])
            #time[status_num] = np.reshape(val[status_num], (-1)) # = (val[status_num][0])
            for i in range(len(val[status_num][0])):
                time_T_n.append((time[status_num][0][i]))
                temp_T_n.append(val[status_num][0][i])
                print(time_T_n)
                print(np.asarray(time_T_n).shape)
        #np.reshape(time[status_num][0], (1, -1)) # = (val[status_num][0])
            #time[status_num][0] = np.transpose(time[status_num][0])
        #print(time_T_n)
        #print(val[status_num])
        #print(np.asarray(val[status_num][0]).shape)
        #print(range(len(val[status_num][0])))
            #print(range(len(val[status_num])))
            #print(time_T_n)

        #lms estimate for whole dataset(generallly bad estimate)
        opt_full, pcov = curve_fit(func6, time_T_n, temp_T_n)     
        time_temp = [0]
        temp_temp = [0]
        for status_num in range(len(mat['status'])):

            if(len(time[status_num])>0 and status_num < len(mat['status'])-1):
                iii = 1;
                while(len(time[status_num+iii])<2 and status_num+iii < len(mat['status'])-1):
                    iii=1+iii
                #if we reached the end and its zero, go back till its non-zero
                while(len(time[status_num+iii])<1):
                    iii=iii-1
                #print(f'{status_num}: len -1 = {len(time[status_num])} len +iii = {len(time[status_num+iii])}')
                time_temp = time_T_n[time_T_n.index(time[status_num][0]):time_T_n.index(time[status_num+iii][-1])]
                #print(len(time_temp));
                temp_temp = temp_T_n[time_T_n.index(time[status_num][0]):time_T_n.index(time[status_num+iii][-1])]
            #if we only have one status_num of non-zero length
            elif(len(time[status_num])>0 ):
                time_temp = time[status_num]
                #print(len(time_temp));
                temp_temp = val[status_num]
                #print(len(temp_temp));
            func_order[status_num] = func6;
            if(len(time_temp)>1):
                    #override_6 = ['T501']
                    #override_1 = ['T502']
                    #override_6_l = [-1.00210016*pow(1, -11)]
                #print(f"what the fuck len time {len(time_temp)}, time {len(temp_temp)}")
                opt[status_num], pcov = curve_fit(func6, time_temp, temp_temp)
                new_func = 0
                if(abs(opt[status_num][0]) < pow(10,-9)):
                    if not ((id_override) in override_6 and (status_num in override_6_state)):
                        #print(f'hello 6 was too high {opt[status_num]}')
                        opt[status_num], pcov = curve_fit(func4, time_temp, temp_temp)
                        #print(f'hello new 4 {opt[status_num]}')
                        func_order[status_num] = func4;
                        new_func = 1
                override_4_found = ((id_override in override_4) and (status_num in override_4_state))

                if((new_func and abs(opt[status_num][0]) < pow(10,-5)) or override_4_found):
                    opt[status_num], pcov = curve_fit(func2, time_temp, temp_temp)
                    func_order[status_num] = func2;
                if(new_func and not override_4_found and abs(opt[status_num][0]) < pow(10,-5) or id_override in override_1):
                    try:
                        opt[status_num], pcov = curve_fit(func1, time_temp, temp_temp)
                        func_order[status_num] = func1;
                    except: 
                        print("failed, exponenential, keeping quadratic")
                        opt[status_num], pcov = curve_fit(func2, time_temp, temp_temp)
                        func_order[status_num] = func2;
                print(f'status_num {status_num} opt: {opt[status_num]}, zero: {abs(opt[status_num][0])} {id_override} {(id_override in override_4)}  {(status_num in override_4_state)}')
        return opt_full, opt, func_order

    t501_opt_full, t501_opt, func_fit_T501 = calc_lmsq(mat['time_T501'], mat['temp_T501'], 'T501')
    t502_opt_full, t502_opt, func_fit_T502 = calc_lmsq(mat['time_T502'], mat['temp_T502'], 'T502')
    t504_opt_full, t504_opt, func_fit_T504 = calc_lmsq(mat['time_T504'], mat['temp_T504'], 'T504')

    #take our lms function and use it to generate a line for all the original timestamps,
    #filling in our missing data, and smoothing the stuff we have
    def plot_holey_data_and_estimate(time_hole, val_hole, time_whole, status_num, opt, func_order):
        axs.plot((time_hole[status_num]), (val_hole[status_num]))
        if(len(time_hole[status_num])>10 and status_num < len(mat['status'])-1):
            iii = 1;
            time_start_index = time_whole.index(time_hole[status_num][0]);
            while(len(time_hole[status_num+iii])<2 and iii + status_num < len(mat['status'])-1):
                iii=1+iii
            ##this'll happen if the last mat['status'] num is empty, so we gotta force it to 
            ##plot to the end of the time
            if(len(time_hole[status_num+iii])<1):
                time_stop_index = time_whole.index(time_whole[-1])
            else:
                time_stop_index = time_whole.index(time_hole[status_num+iii][-1])
                print(f"end of time {time_stop_index} mat['status'] num: {status_num+iii}")
            time_temp = time_whole[time_start_index:time_stop_index]
            #print(f'ext num: {ii} len {len(time_temp)}')
            axs.plot(time_temp, func_order[status_num](time_temp, *opt[ii]), linestyle='dashed')

    #really probably dont need this anymore but WHATEVER
    i = []
    for status_num in range(len(mat['status'])):
        if(mat['status'][status_num] == "'ACTUATOR_ON'" or mat['status'][status_num] == "'ACTUATOR_OFF'"):
            i.append(status_num)
    fig, ax = plt.subplots(2, 3);

    for ii in i:
        #print(ii)
        #fig, ax = plt.subplots(2, 3);
        fig.suptitle(f'Plot for ACTUATOR_ON {ii}')
        axs = ax[0,0];
        plot_holey_data_and_estimate(mat['time_T501'], mat['temp_T501'], mat['time_T501_S'], ii, t501_opt, func_fit_T501)
        axs.set_title('Temp T501 (degC)vs Time')
        axs = ax[0,1];
        plot_holey_data_and_estimate(mat['temp_T502'], mat['temp_T502'], mat['time_T502_S'], ii, t502_opt, func_fit_T502)
        axs.set_title('Temp T502 (degC)vs Time')
        axs = ax[0,2];
        axs.plot((mat['temp_T504'][ii]), (mat['temp_T504'][ii]))
        plot_holey_data_and_estimate(mat['time_T504'], mat['temp_T504'], mat['time_T504_S'], ii, t504_opt, func_fit_T504)
        axs.set_title('Temp T504 (degC) vs Time')
        axs = ax[1,0];
        axs.plot((mat['time_voltage'][ii]), (mat['voltage'][ii]))
        axs.set_title('Nicrome Voltage vs Time')
        axs = ax[1,1];
        axs.plot((mat['time_current'][ii]), (mat['current'][ii]))
        axs.set_title('Nicrome current (A_) vs Time')
    #axs = ax[0,0];
    #axs.plot(time_T, func(time_T, *t501_opt_full))

    plt.savefig(f'{arg.filepath.split(".txt")[0]}_ACTUATOR_ALL_FILTERED.png')
    t_est = np.linspace(time_stamp_good[0], time_stamp_good[1], int(abs((time_stamp_good[1]-time_stamp_good[0])/t_div)))

    def find_nearest(array,value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return idx-1
        else:
            return idx
        
    def val_estimate(time, val, opt, func_order, time_hole):
        time_stop_index = 0;
        time_temp = []
        for status_num in range(len(mat['status'])):
            iii = 0
            #make sure we predict for the whole range
            if(status_num == 0):
                time_start_index = find_nearest(time, time_stamp_good[0]);
            else:
                time_start_index = time_stop_index
            if(len(time_hole[status_num])>0 and status_num < len(mat['status'])-1):
                iii = 1;
                while(len(time_hole[status_num+iii])<0 and iii + status_num < len(mat['status'])-1):
                    iii=1+iii
                
                if(len(time_hole[status_num+iii])<1):
                    time_stop_index = len(time)-1
                else:
                    time_stop_index = find_nearest(time,time_hole[status_num+iii][0])
            
            ##this'll happen if the last mat['status'] num is empty, so we gotta force it to 
            ##plot to the end of the time
            elif(len(time_hole[status_num])>0 and status_num >= len(mat['status'])-1):
                time_stop_index = find_nearest(time,time_hole[status_num][-1])

            if(time_stop_index > time_start_index):
                time_temp = time[time_start_index:time_stop_index]

                #print(f'ext num: {ii} len {len(time_temp)}')
                #axs.plot(time_temp, func_order[status_num](time_temp, *opt[ii]), linestyle='dashed')
                #for x in func_order[status_num](time_temp, *opt[status_num]):
                val[time_start_index:time_stop_index] = func_order[status_num](time_temp, *opt[status_num])
                print(val[0])
                print(f"start of time {time_start_index} end of time {time_stop_index} len time {len(time_temp)} len vals {len(val)} mat['status'] num: {status_num+iii}")


    val_t501_est = [0 for _ in t_est]
    val_t502_est = [0 for _ in t_est]
    val_t504_est = [0 for _ in t_est]

    val_estimate(t_est, val_t501_est, t501_opt, func_fit_T501, mat['time_T501'])
    val_estimate(t_est, val_t502_est, t502_opt, func_fit_T502, mat['time_T502'])
    val_estimate(t_est, val_t504_est, t504_opt, func_fit_T504, mat['time_T504'])
    print(f'final len 501 {len(t_est)} {len(val_t501_est)}')
    print(f'final len 502 {len(t_est)} {len(val_t502_est)}')
    print(f'final len 504 {len(t_est)} {len(val_t504_est)}')

    def print_time_val(time, val):
        f.write(str(time))
        f.write(',')
        f.write(str(val))
        f.write(',')
    with open(f'{arg.filepath.split(".txt")[0]}_CLEAN_EST.csv', 'w') as f:
        #print(small_len(status_num))    
        for i in range(len(t_est)):
            print_time_val(t_est[i], val_t501_est[i])
            print_time_val(t_est[i], val_t502_est[i])
            print_time_val(t_est[i], val_t504_est[i])
            f.write('\n')

    for status_num in range(len(mat['status'])):
        #print(ii)
        #fig, ax = plt.subplots(2, 3);
        fig.suptitle(f'Plot for ACTUATOR_ON ACTUATOR_ALL_ESTIMATED')
        axs = ax[0,0];
        axs.plot((mat['time_T501'][status_num]), (mat['temp_T501'][status_num]))
        axs.set_title('Temp T501 (degC)vs Time')
        axs = ax[0,1];
        axs.plot((mat['temp_T502'][status_num]), (mat['temp_T502'][status_num]))
        axs.set_title('Temp T502 (degC)vs Time')
        axs = ax[0,2];
        axs.plot((mat['temp_T504'][status_num]), (mat['temp_T504'][status_num]))
        axs.set_title('Temp T504 (degC) vs Time')
        axs = ax[1,0];
        axs.plot((mat['time_voltage'][status_num]), (mat['voltage'][status_num]))
        axs.set_title('Voltage Nicrome vs Time')
        axs = ax[1,1];
        axs.plot((mat['time_current'][status_num]), (mat['current'][status_num]))
        axs.set_title('Nicrome current (A_) vs Time')
    axs = ax[0,0];
    axs.plot(t_est[0:len(t_est)], val_t501_est, linestyle='dashed')
    axs = ax[0,1];
    axs.plot(t_est[0:len(t_est)], val_t502_est, linestyle='dashed')
    axs = ax[0,2];
    axs.plot(t_est[0:len(t_est)], val_t504_est, linestyle='dashed')


    plt.savefig(f'{arg.filepath.split(".txt")[0]}_ACTUATOR_ALL_ESTIMATED.png')