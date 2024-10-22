import linear_actuator_coordinate as lc
# import arduino_control as ac
import numpy as np
import serial
import struct
import time
import subprocess
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime
from scipy.optimize import fsolve
import virgo

#-------------------------------------------------------------------------------------------------
# test the whole code structure here: should generate series of (n,3) arrays
#-------------------------------------------------------------------------------------------------
star_name = 'Regulus'
position = lc.celestial_obj_position(star_name) # find the place of the celestial body in sky
goal_orientation = lc.find_p_r(position) # find the orientation the telescope need to be
orientation_path = lc.move_to_point(goal_position=np.degrees(goal_orientation), current_angle=[0,0,0]) # input for this function is in degree
print(orientation_path)
print('Benjamin is the smartest')
print(goal_orientation)
leg_lengths = lc.regional_scan_points([goal_orientation[0],goal_orientation[1]], obs_name='orien_data/'+star_name)

#-------------------------------------------------------------------------------------------------
# Test code for interacting with arduinos
#-------------------------------------------------------------------------------------------------
# serials = ac.start_serial()
# for i in leg_lengths:
#     result = ac.send_position(leg_lengths)
