import numpy as np
import serial
import subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime


# Constants
BASE_RADIUS = 91.61  # cm
TOP_RADIUS = 47  # cm
INITIAL_HEIGHT = 166.5  # cm

def rotation_matrix(pitch, roll, yaw):
    """Create 3D rotation matrix from pitch, roll, and yaw angles (in radians)."""
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(pitch), -np.sin(pitch)],
                   [0, np.sin(pitch), np.cos(pitch)]])
    
    Ry = np.array([[np.cos(roll), 0, np.sin(roll)],
                   [0, 1, 0],
                   [-np.sin(roll), 0, np.cos(roll)]])
    
    Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                   [np.sin(yaw), np.cos(yaw), 0],
                   [0, 0, 1]])
    
    return Rz @ Ry @ Rx

def calculate_leg_lengths(pitch, roll, yaw):
    """Calculate leg lengths for given orientation angles."""
    
    # Base attachment points (120 degrees apart)
    base_points = np.array([
        [BASE_RADIUS * np.cos(0), BASE_RADIUS * np.sin(0), 0],
        [BASE_RADIUS * np.cos(2*np.pi/3), BASE_RADIUS * np.sin(2*np.pi/3), 0],
        [BASE_RADIUS * np.cos(4*np.pi/3), BASE_RADIUS * np.sin(4*np.pi/3), 0]
    ])
    
    # Top attachment points (120 degrees apart)
    top_points = np.array([
        [TOP_RADIUS * np.cos(0), TOP_RADIUS * np.sin(0), INITIAL_HEIGHT],
        [TOP_RADIUS * np.cos(2*np.pi/3), TOP_RADIUS * np.sin(2*np.pi/3), INITIAL_HEIGHT],
        [TOP_RADIUS * np.cos(4*np.pi/3), TOP_RADIUS * np.sin(4*np.pi/3), INITIAL_HEIGHT]
    ])
    
    # Calculate rotation matrix
    R = rotation_matrix(pitch, roll, yaw)
    
    # Transform top attachment points
    transformed_top_points = (R @ (top_points - [0, 0, INITIAL_HEIGHT]).T).T + [0, 0, INITIAL_HEIGHT]
    
    # Calculate leg vectors
    leg_vectors = transformed_top_points - base_points
    
    # Calculate leg lengths
    leg_lengths = np.linalg.norm(leg_vectors, axis=1)
    
    if np.all(leg_lengths<223) and np.all(leg_lengths> 122):
        return leg_lengths
    else:
        return ValueError
    
def celestial_obj_position(object_name):
    '''alculate the azimuth and altitude of the celestial object accoridng to 
    the object name. 
    object_name: str
    alt, az: float64'''
    # Define the location (lat, lon, and height)
    location = EarthLocation(lat=40.36 * u.deg, lon=-74.67 * u.deg, height=61.87 * u.m)
    
    # Get the current time in UTC (Astropy Time object)
    time = Time(datetime.now())
    
    # Create an AltAz frame with the current time and location
    altaz = AltAz(obstime=time, location=location)
    
    # Get the object's position based on name and transform it to AltAz coordinates
    try:
        obj = get_body(object_name, time, location)  # Use for solar system objects (e.g., 'mars')
    except:
        obj = SkyCoord.from_name(object_name)  # Use for other objects (e.g., 'Andromeda Galaxy')
    
    # Transform to AltAz frame
    obj_altaz = obj.transform_to(altaz)
    
    # Extract the altitude and azimuth in degrees
    alt = obj_altaz.alt.degree
    az = obj_altaz.az.degree
    
    # Output or use the altitude and azimuth values
    print(f"{object_name.title()}'s altitude: {alt:.2f} degrees")
    print(f"{object_name.title()}'s azimuth: {az:.2f} degrees")

    return [az,alt]

def turn_to_star(position):
    '''the function turn the radio telescope to point at certain position through
    calculating the length of each leg with respect to the direction of the celestial
    body.
    arg: 
    position: list with 2 elements
    return list with 3 elements'''
    az, alt = position
    x, y = np.arctan(np.sin(az)/np.sin(alt)), np.arctan(np.sin(az)/np.cos(alt))
    return x, y, 0
    
def get_serial():
    '''This function is aimed at getting the serial number of three connection ports, 
    therfore, the output will be a list of str elements.'''
    addresses = subprocess.run("ls /dev/tty*", shell=True, capture_output=True, text=True).stdout.splitlines() #to get return and split the output to lines
    arduinos = []
    for i in addresses:
        if i[0:15] == "/dev/tty.usbmod":
            arduinos.append(i)
    if len(arduinos) == 3:
        return arduinos
    else:
        return "Incorrect Arduino number"
    
def send_serial(arduinos, final_position):
    '''This function send serial information to the arduinos'''
    
    
def move_by_step(original_position, goal_position, step_time = 5):
    '''This function moves a the orientation of the radio telescope from pointing at
    the last place to the next place
    args: 
    original_position: the legs' lengths at the start, list
    goal_position: the legs' lengths at the final place, list
    step_time: time taken by the linear actuator to move forward in each step
    return:
    result: finished or not; boolean
    final_position: the position calculated by hall effect sensor, list'''
    #Calculate the maximum angle move in time range by implementing the speed limit 2.16inches/s
    distance_max = np.max(goal_position-original_position)
    total_time = distance_max/2.16
    while (distance_max>100):
        
        send_serial(arduinos=get_serial, )
    return 'rishi is so handsome'


def reach_aiming_length():
    '''this function return the aiming length for linear actuators
    to focus on the correct celestial object in the Universe'''

    return

def regional_sky_scan(frequency, center_azalt, averaged_time):
    '''make a regional sky scan that later return an tuple of data, 
    connecting the orientation with the intensity in certain frequency'''
    return
# Example usage
pitch = np.radians(55)  # 20 degrees pitch
roll = np.radians(10)   # 10 degrees roll
yaw = np.radians(0)     # 5 degrees yaw

leg_lengths = calculate_leg_lengths(pitch, roll, yaw)[0]
print(leg_lengths)
celestial_obj_position('Cygnus A')
position1 = celestial_obj_position('sun')
a, b, c = turn_to_star(position1)
print(calculate_leg_lengths(a, b, c))