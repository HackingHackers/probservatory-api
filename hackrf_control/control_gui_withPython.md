This is a mark down file notes how to control GUI or GNUradio and pass parameteres through python
## Notice, the Virgo(radio telescope analysis coding system); helped the system use command line to pass argument
## This method of passing through command line can be done through Argparse() function in Python

Penn radio telescope:
1. How data is stored as files?
The data is saved in npz file type. 
# example of command
def cmdScan():
    ser.write(SCAN)
    deg = input("Enter number of degrees to turn: ")
    print("Sending " + deg)
    ser.write(str.encode(deg))
    print("Reading data")
    ndata = readStream(ser)
    current_state = readState(ser)
    PrintState()
    # Convert
    # ndata = numpyState(ndata)
    # Save
    np.savez(file=time.ctime().replace(' ', '_') + '.npz',
             ndata=ndata)
    # Plot
    PlotData(ndata)

RTL_SDR can be used to scan through some easy functions: 
# sdr.center_freq = freq
Upon is the line that make the obervation

2. How to analysis data and create 2D file?
2D file is created through analyzing the signal intensity relative to the position;
therefore, storing the az/alt and average signal intensity can be helpful on this survey

MIT Radio Telescope:
1. What graphs are created by the telescope?
a 21 cm hydrogen survey of the sky at galaxy corss section plain
2. How are the graphs created?
The graph is created through value correction from C++ and later being made through matlab
3. How can we create similar graph?
We need to figure out how to make hackrf give up precise value on this range, get the positional value from arduino feedback
and later put the peak of data and positional value into matlab and make such a graph.

Current Questions:
1. Control the Hackrf with python or GNU radio
2. Python code to control the three linear actuators
s