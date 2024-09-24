import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from astropy.time import Time
from astropy.coordinates import get_sun, AltAz, EarthLocation
import astropy.units as u

# Constants
BASE_RADIUS = 91.61 / 100  # meters
TOP_RADIUS = 47 / 100  # meters
INITIAL_HEIGHT = 166.5 / 100  # meters
SUN_DISTANCE = 3  # Scaled distance of the Sun

# Define the positions of the base joints
base_joints = np.array([
    [BASE_RADIUS * np.cos(np.radians(angle)), BASE_RADIUS * np.sin(np.radians(angle)), 0]
    for angle in np.arange(0, 360, 120)
])

# Get Sun's positional data
location = EarthLocation(lat=40.36 * u.deg, lon=-74.67 * u.deg, height=61.87 * u.m)
times = Time('2023-07-19 10:00:00') + np.arange(0, 10, 0.05) * u.hour  # from 10 AM to 8 PM
altaz = AltAz(obstime=times, location=location)
sun = get_sun(times).transform_to(altaz)
altitudes = sun.alt.degree
azimuths = sun.az.degree

# Initialize the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([0, 3.5])

# Create the sun and manipulator
sun_plot, = ax.plot([], [], [], 'yo', markersize=10)
manipulator_base, = ax.plot(base_joints[:, 0], base_joints[:, 1], base_joints[:, 2], 'bo', markersize=8)
manipulator_platform, = ax.plot([], [], [], 'ro', markersize=6)
legs = [ax.plot([], [], [], 'g-')[0] for _ in range(3)]
top_platform_lines = [ax.plot([], [], [], 'b-')[0] for _ in range(3)]  # Changed color to blue

# Initialize text annotations for leg lengths
leg_length_texts = [ax.text2D(0.05, 0.95 - 0.05 * i, '', transform=ax.transAxes) for i in range(3)]

# Initialize the objects
def init():
    sun_plot.set_data([], [])
    sun_plot.set_3d_properties([])
    manipulator_platform.set_data([], [])
    manipulator_platform.set_3d_properties([])
    for leg in legs:
        leg.set_data([], [])
        leg.set_3d_properties([])
    for line in top_platform_lines:
        line.set_data([], [])
        line.set_3d_properties([])
    for text in leg_length_texts:
        text.set_text('')
    return sun_plot, manipulator_base, manipulator_platform, *legs, *top_platform_lines, *leg_length_texts

# Calculate rotation matrix from Sun direction
def rotation_matrix_from_sun(az, alt):
    """Create a rotation matrix to orient the top platform towards the Sun."""
    sun_direction = np.array([
        np.cos(np.radians(alt)) * np.sin(np.radians(az)),
        np.cos(np.radians(alt)) * np.cos(np.radians(az)),
        np.sin(np.radians(alt))
    ])
    sun_direction /= np.linalg.norm(sun_direction)
    
    z_axis = np.array([0, 0, 1])
    rotation_axis = np.cross(z_axis, sun_direction)
    rotation_angle = np.arccos(np.dot(z_axis, sun_direction))
    
    K = np.array([
        [0, -rotation_axis[2], rotation_axis[1]],
        [rotation_axis[2], 0, -rotation_axis[0]],
        [-rotation_axis[1], rotation_axis[0], 0]
    ])
    
    R = np.eye(3) + np.sin(rotation_angle) * K + (1 - np.cos(rotation_angle)) * np.dot(K, K)
    
    return R

def calculate_top_joints(az, alt):
    R = rotation_matrix_from_sun(az, alt)
    platform_center = np.array([0, 0, INITIAL_HEIGHT])
    top_joints = np.array([
        platform_center + R @ [TOP_RADIUS * np.cos(np.radians(angle)), TOP_RADIUS * np.sin(np.radians(angle)), 0]
        for angle in np.arange(0, 360, 120)
    ])
    return top_joints

def calculate_leg_lengths(top_joints):
    leg_vectors = top_joints - base_joints
    leg_lengths = np.linalg.norm(leg_vectors, axis=1)
    return leg_lengths

# Animation function
def animate(i):
    az = azimuths[i]
    alt = altitudes[i]
    
    # Update sun position
    sun_x = SUN_DISTANCE * np.cos(np.radians(az)) * np.cos(np.radians(alt))
    sun_y = SUN_DISTANCE * np.sin(np.radians(az)) * np.cos(np.radians(alt))
    sun_z = SUN_DISTANCE * np.sin(np.radians(alt)) + 1.0  # Increase z-coordinate for visibility
    sun_plot.set_data(sun_x, sun_y)
    sun_plot.set_3d_properties(sun_z)
    
    # Update manipulator platform
    top_joints = calculate_top_joints(az, alt)
    manipulator_platform.set_data(top_joints[:, 0], top_joints[:, 1])
    manipulator_platform.set_3d_properties(top_joints[:, 2])
    
    # Update legs
    for j, leg in enumerate(legs):
        leg.set_data([base_joints[j, 0], top_joints[j, 0]], [base_joints[j, 1], top_joints[j, 1]])
        leg.set_3d_properties([base_joints[j, 2], top_joints[j, 2]])
    
    # Link the top joints (forming a triangle)
    top_platform_lines[0].set_data([top_joints[0, 0], top_joints[1, 0]], [top_joints[0, 1], top_joints[1, 1]])
    top_platform_lines[0].set_3d_properties([top_joints[0, 2], top_joints[1, 2]])

    top_platform_lines[1].set_data([top_joints[1, 0], top_joints[2, 0]], [top_joints[1, 1], top_joints[2, 1]])
    top_platform_lines[1].set_3d_properties([top_joints[1, 2], top_joints[2, 2]])

    top_platform_lines[2].set_data([top_joints[2, 0], top_joints[0, 0]], [top_joints[2, 1], top_joints[0, 1]])
    top_platform_lines[2].set_3d_properties([top_joints[2, 2], top_joints[0, 2]])
    
    # Update leg lengths
    leg_lengths = calculate_leg_lengths(top_joints)
    for j, text in enumerate(leg_length_texts):
        text.set_text(f'Leg {j+1} Length: {leg_lengths[j]:.2f} m')
    
    return sun_plot, manipulator_base, manipulator_platform, *legs, *top_platform_lines, *leg_length_texts

# Create the animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(times), interval=50, blit=True)

plt.show()
