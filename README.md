# eject3d
3d version of eject with a wrapper that generates a large number of block ejections in a Monte Carlo style

This python code calculates the 3d trajectories of multiple blocks being ejected from a volcanic crater.  The user can specify a range of initial velocities, block sizes, and trajectory angles in the horizontal and vertical direction.  The model writes out the location of the landing point (x,y,z), block travel time, and, optionally, for each ejection, a table of the flight path along with atmospheric properties along that path.

Run the model as follows:
1)  In eject3dfun.py, modify parameters in Block 1.  These include the thermal properties in the atmosphere, wind, block density, and whether or not to write out separate trajectory files for each block.  
2)  Modify Block 1 of the file call_eject3dfun.py, to specify the number of block to simulate, and the range of sizes, velocities, and trajectory angles simulated.  Block 1 of call_eject3dfun.py is given below.
3)  Run the model by typing "python call_eject3dfun.py" at the command line.
4)  Examine the output files.  The file call_eject3dfun.py should write out a table showing the x,y,z landing location of each block.  Each row in the table is a different simulation.
     
The model uses a solution method and block shapes that are described in Mastin (2001).

Reference:
Mastin, L.G., 2001, A simple calculator of ballistic trajectories for blocks ejected during volcanic eruptions, U.S. Geological Survey Open-File Report 01-45, U.S. Geological Survey, p. 26. (https://vhub.org/resources/455)
