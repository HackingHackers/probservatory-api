import numpy as np
import serial
import struct
import time
import subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime
from scipy.optimize import fsolve
import virgo

def get_serial():
    '''This function get serial numbers of three connection ports, 
    output str elements.'''
    addresses = subprocess.run("ls /dev/tty.usb*", shell=True, capture_output=True, text=True).stdout.splitlines() #to get return and split the output to lines
    if len(addresses)<=3:
        return addresses
    else:
        return 'incorrect amount of arduino connected'

# serial_nums = get_serial()

def start_serial(a, b, c):
    ser1 = serial.Serial(a, 9600)  
    ser2 = serial.Serial(b, 9600)  
    ser3 = serial.Serial(c, 9600) 
    return ser1, ser2, ser3

# ser1, ser2, ser3 = start_serial(serial_nums[0], serial_nums[1], serial_nums[2])
ser1, ser2, ser3 = 1,2,3

def send_position(positions = [], serials = [ser1, ser2, ser3], time = 5):
    '''This function send serial information to the arduinos
    it should be only used once and do not change the port number, then, the arduino won't be messed up
    position: list with 3 numbers defining the length of linear actuators
    original: if or not the numbers are initial points
    timer: bool type to determine whether the serial represent the time arduino take to reach a position
    serials: a list of started serials'''
    # config the ports
    data = [0,0,0]
    for i in positions:
        ser1.write(str(i[0]).encode(encoding='utf8'))
        ser2.write(str(i[1]).encode(encoding='utf8'))
        ser3.write(str(i[2]).encode(encoding='utf8'))
        time.sleep(time+1)
        output1 = float(ser1.readline().decode('utf-8').strip())  # Decode from bytes to string
        output2 = float(ser2.readline().decode('utf-8').strip()) # Decode from bytes to string
        output3 = float(ser3.readline().decode('utf-8').strip())  # Decode from bytes to string
        data = [output1, output2, output3]
        print(data)
    return data # return the current real position of the linear actuators'

def send_init_position(positions = [], serials = [ser1, ser2, ser3]):
    '''send position of linear actuator'''
    # config the ports
    data = [0,0,0]
    ser1.write('000'+str(positions[0]).encode(encoding='utf8'))
    ser2.write('000'+str(positions[1]).encode(encoding='utf8'))
    ser3.write('000'+str(positions[2]).encode(encoding='utf8'))
    output1 = float(ser1.readline().decode('utf-8').strip())  # Decode from bytes to string
    output2 = float(ser2.readline().decode('utf-8').strip()) # Decode from bytes to string
    output3 = float(ser3.readline().decode('utf-8').strip())  # Decode from bytes to string
    data = [output1, output2, output3]
    return data # return the current real position of the linear actuators'

def send_timer(positions = [], serials = [ser1, ser2, ser3]):
    ser1.write('00'+str(positions[0]).encode(encoding='utf8'))
    ser2.write('00'+str(positions[1]).encode(encoding='utf8'))
    ser3.write('00'+str(positions[2]).encode(encoding='utf8'))
    return 'successfully sent'


def save_length(len1, len2, len3, file_name = 'leg_lengths'):
    data = [len1, len2, len3]
    np.savetxt(file_name, data, delimiter=',') 
    return 'success'

def read_length(filename = 'leg_lengths'):
    data = np.genfromtxt(filename, delimiter=',', dtype=float, encoding='utf-8')
    return data

def read_position():
    '''this function returns the lengths of three linear actuators'''
    ser1, ser2, ser3 = start_serial()
    output1 = float(ser1.readline().decode('utf-8').strip())  # Decode from bytes to string
    output2 = float(ser2.readline().decode('utf-8').strip()) # Decode from bytes to string
    output3 = float(ser3.readline().decode('utf-8').strip())  # Decode from bytes to string
    return output1, output2, output3

# moving_positions = np.radians(np.genfromtxt('psrt_obs.csv', delimiter=','))

# leg_lengs = []
# for i in moving_positions:
#     leg_lengs.append(lc.calculate_leg_lengths(i[0], i[1], 0))

# leg_lengs = np.array(leg_lengs)
# # data_result = []
# # for i in leg_lengs:
# #     result = send_serial(positions=i, original=False, serials=moving_positions)
# #     data_result.append(result)

# # print(data_result)
# num = ['0006',90,80,70,60,50]
# print(leg_lengs)


# # ser1 = serial.Serial(serial_nums[0], 9600)  

# # time.sleep(2) # delay to sutup the machine, need 2 seconds to setup
# # for i in num:
# #     ser1.write(str(i).encode(encoding='utf8'))
# #     time.sleep(3)

# # print('finished')