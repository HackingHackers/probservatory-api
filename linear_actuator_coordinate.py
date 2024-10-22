import numpy as np
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime
from scipy.optimize import fsolve
import arduino_control as ac

# Constants
BASE_RADIUS = 96.61  # cm; back up, 96.61cm
TOP_RADIUS = 47  # cm
INITIAL_HEIGHT = 170.7  # cm


'''calculation of leg lengths according to celestial objects' location'''
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

def find_p_r(position):
    '''the function return the pitch, roll and etc.
    arg: 
    position: list with 2 elements
    return list with 3 elements'''
    az, alt = position
    x, y = np.arctan(np.sin(az)/np.sin(alt)), np.arctan(np.sin(az)/np.cos(alt))
    return x, y, 0

def intepret_position(location1, location2, location3):
    '''this one returns the angles of the location in radian'''
    return fsolve(calculate_leg_lengths(), [location1, location2, location3])


'''Movement functions turn objects to places'''
def move_to_point(goal_position = [], current_angle = [], step_time = 5, step_scale = 1):
    '''This function moves a the orientation of the radio telescope from pointing at
    the last place to the next place -> this is a maximum speed moving function
    args: 
    original_position: the legs' lengths at the start, list

    goal_position: angle of the final position --- in DEGREE!

    step_time: time taken by the linear actuator to move forward in each step
    return:
    result: finished or not; boolean
    final_position: the position calculated by hall effect sensor, list'''
    #Calculate the maximum angle move in time range by implementing the speed limit 2.16inches/s
    goal_position, current_angle = np.array(np.radians(goal_position)), np.array(np.radians(current_angle))
    step_num = np.radians(step_scale*step_time)
    # to a universal minimum, 11 degree of total change maximum it can achieve. Just divide the goal to number of steps in each smaller than 11 degrees
    num_turns = int(np.max(np.abs(goal_position-current_angle))/(step_num)+2)
    pitch1 = np.linspace(current_angle[0], goal_position[0], num_turns)[1::]
    roll1 = np.linspace(current_angle[1], goal_position[1], num_turns)[1::]
    yaw1 = np.zeros(num_turns-1)
    print(pitch1, roll1)
    # pitch1, roll1, yaw1 = np.radians(pitch1), np.radians(roll1), np.radians(yaw1)
    leg_lss = []
    if len(pitch1)>1:
        for i in range(len(pitch1)):
            leg_lss.append(calculate_leg_lengths(pitch1[i], roll1[i], yaw1[i]))
    if leg_lss[:][:]!=ValueError:
        return leg_lss
    return ValueError



def regional_scan_points(center_point = [0,0,0], az_range = 15, alt_range = 15 , resolution = [15,15], obs_name = 'psrt_obs'):
    '''calculate the route of a sky scan one by one

    Changing to input in degree
    change input parameters to the az and alt of the star

    input:
    center_point: a celestial object axis position
    figure_range: four values defining the degree of the figure
    resolution: list, 2 elements specifying horizontal and vertical resolution'''
    
    # first, check if the orientation of telescope is to the object; orient to it if not
    try:
        start_lenths = ac.read_position() # need to think of another way to insert position here
    except:
        start_lenths= [0,0,0]
    #test the angular flexibility of the telescope; make it able to turn to the point
    center_point = find_p_r(np.radians(np.array(center_point)))
    leg_lses = calculate_leg_lengths(center_point[0], center_point[1], 0)
    angle_diff = max(abs(leg_lses-start_lenths))
    if(angle_diff>0.1):
        move_to_point(center_point, start_lenths)
    
    # start to plan out the matrix of movement goals to skan in sky; I
    resolution = np.int16(np.array(resolution)/2)*2+1
    az_step = (az_range)/resolution[0]
    alt_step = (alt_range)/resolution[1]
    az_nsteps = np.linspace(0, resolution[0], resolution[0]+np.int16(az_step))
    alt_nsteps = np.linspace(0, resolution[1], resolution[1]+np.int16(alt_step))

    point_bottomleft = [center_point[0]+resolution[0]/2*az_step, center_point[1]+resolution[1]/2*alt_step]
    # print('buttomleft:',point_bottomleft)
    test_point = np.radians(point_bottomleft)

    if calculate_leg_lengths(test_point[0], test_point[1], 0).any() == ValueError:
        return ValueError + 'out of scannable region'
    
    # it should make the process easier if the scanning range is following a turning path (up,down,up while moving to right) 
    az_lins = point_bottomleft[0]-az_step*az_nsteps
    alt_lins = point_bottomleft[1]-alt_nsteps*alt_step
    obs_points = np.ones((len(az_nsteps), len(alt_nsteps), 2)) #make a n*m*2 array to store all the locations on the path
    obs_points[:,:,0] = obs_points[:,:,0]*az_lins+np.degrees(center_point[0])
    obs_points = np.transpose(obs_points, axes=[1,0,2])
    obs_points[:,:, 1] = obs_points[:,:,1]*alt_lins+np.degrees(center_point[1])
    obs_points[1::2] = obs_points[1::2, ::-1, :] #change the order of the points to make them zigzagging across rows
    obs_points = np.float16(np.reshape(obs_points, (len(az_nsteps)*len(alt_nsteps), 2)))
    # convert the az/alt into leg lengths
    obs_points = np.radians(obs_points)
    leg_lss = np.zeros([len(obs_points), 3])
    for i in range(len(obs_points)):
        leg_lss[i,:] = calculate_leg_lengths(obs_points[i,0],obs_points[i,1],0)
    print(np.shape(leg_lss))
    np.savetxt(obs_name+'.csv', leg_lss, delimiter=",") 
    return leg_lss

def star_oriented(obj_name):
    '''make a regional sky scan that later return an tuple of data, 
    connecting the orientation with the intensity in certain frequency'''
    obj_pos = celestial_obj_position(obj_name)
    turn_angles = find_p_r(obj_pos)
    final_length = move_to_point(turn_angles, fsolve(calculate_leg_lengths, ac.read_position))
    print("oriented to star")
    return "oriented to star"


#-------------------------------------------------------------------------------------------------
# Test code 
#-------------------------------------------------------------------------------------------------

# pitch = np.radians(55)  # 20 degrees pitch
# roll = np.radians(10)   # 10 degrees roll
# yaw = np.radians(0)     # 5 degrees yaw

# leg_lengths = calculate_leg_lengths(pitch, roll, yaw)[0]
# print(leg_lengths)
# celestial_obj_position('Cygnus A')
# position1 = celestial_obj_position('sun')
# a, b, c = turn_to_star(position1)
# print(calculate_leg_lengths(a, b, c))
# print(move_to_point([pitch, roll, yaw],[0,0,0]))


