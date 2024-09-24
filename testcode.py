import serial
import subprocess

addresses = subprocess.run("ls /dev/tty*", shell=True, capture_output=True, text=True).stdout.splitlines()

arduinos = []
for i in addresses:
    if i[0:15] == "/dev/tty.usbmod":
        arduinos.append(i)


