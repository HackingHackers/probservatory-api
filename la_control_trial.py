import serial
import time
import pyautogui
import PySimpleGUI as sg

# Initialize the Arduino serial connection (use the correct port and baudrate)
arduino = serial.Serial(port='/dev/cu.usbmodem101', baudrate=9600, timeout=.1)  # Adjust port to your actual port

# Function to send data to the Arduino via serial
def write_read(x):
    arduino.write(bytes(x, 'utf-8'))  # Write the data to the Arduino
    time.sleep(0.05)  # Small delay to ensure data is sent
    data = arduino.readline()  # Read any response from Arduino (optional)
    return data

# Function to simulate key presses
def press_key(key, press_time):
    pyautogui.keyDown(key)  # Simulate key press down
    time.sleep(press_time)  # Hold the key for a certain duration
    pyautogui.keyUp(key)  # Release the key

# Define the layout for the GUI
layout = [[sg.Text("Control the Actuator")], [sg.Button("RISE")], [sg.Button("RETRACT")]]

# Create the window
window = sg.Window("Actuator Control", layout)

# Main event loop
while True:
    event, values = window.read()

    # Break the loop if user closes the window
    if event == sg.WIN_CLOSED:
        break

    # Handle the "RISE" button press
    elif event == "RISE":
        press_key('w', 0.1)  # Simulate 'w' key press for 0.1 seconds (you can adjust the time)
        write_read('w')  # Send 'w' to Arduino over serial

    # Handle the "RETRACT" button press
    elif event == "RETRACT":
        press_key('s', 0.1)  # Simulate 's' key press for 0.1 seconds
        write_read('s')  # Send 's' to Arduino over serial

# Close the GUI window
window.close()
