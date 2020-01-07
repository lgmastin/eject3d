# eject3d
3d version of eject with a wrapper that generates a large number of block ejections in a Monte Carlo style

This python code calculates the 3d trajectories of multiple blocks being ejected from a volcanic crater.  The user can specify a range of initial velocities, block sizes, and trajectory angles in the horizontal and vertical direction.  The model writes out the location of the landing point (x,y,z), block travel time, and, optionally, for each ejection, a table of the flight path along with atmospheric properties along that path.

Run the model as follows:
1)  In eject3dfun.py, modify parameters in Block 1.  These include the thermal properties in the atmosphere, wind, block density, and whether or not to write out separate trajectory files for each block.  Block 1 is shown below.

  #############################################################################################
      # BLOCK 1.  VOLCANO-SPECIFIC PARAMETERS.
      # These parameters should be reviewed for each series of runs

      dragred = 0                     #zone of reduced drag, meters from vent
      elev = 1000                     #elevation of ejection point
      lapse = 6.5                     #thermal lapse rate, deg. K/km
      rhor = 2500                     #density of ballistic, kg/m3
      wind = 0.                       #wind velocity, m/s
      winddir_deg = 0.                #direction toward which wind is blowing, degrees CW from north
      xi = 0                          #vertical distance of landing point above ejection point, meters
                                      ## Note: Only negative values (<0) are allowed for xi.
      TzeroC = 25                     #Temperature at sealevel, Celsius
      write_outfile = 'yes'           #='yes' if an output file is to be written for each run:
                                      #if write_outfile='yes', the output files have the format:
                                      # eject_out_YY-MM-DD_HH_mm_run###.txt, where YY, MM, DD, HH and mm
                                      # are the year, month, day, hour, and minute of the run series
      filename_stem = 'output/Halemaumau'    #If write_outfile='yes', this gives the stem of filename to be 
                                    #written out.  The full filename will be:
                                    # "filename_stem_run###.txt", where ### is the run number (inow)
  #############################################################################################
                                    
2)  Modify Block 1 of the file call_eject3dfun.py, to specify the number of block to simulate, and the range of sizes, velocities, and trajectory angles simulated.  Block 1 of call_eject3dfun.py is given below.

  #########################################################################################
  #BLOCK 1.  PARAMETERS THAT DEFINE THE RANGE AND DISTRIBUTION OF INITIAL BLOCK SIZES, SHAPES,
  #EJECTION VELOCITIES, AND TRAJECTORY ANGLES

  #The values below define the mean and standard
  #deviation of a Gaussian distribution of values.

  #NUMBER OF BLOCKS TO EJECT
  n_ejects  = 10

  #STEM OF OUTPUT FILENAME
  #The full filename will be filename_stem_YY-MM-DD_HHmm.txt, where
  #YY,MM,DD,HH, and mm are the year, month, day, hour and minute of the run.
  filename_stem = 'output/call_eject3d_out_'

  #INITIAL ANGLE FROM HORIZONTAL (Gaussian distribution)
  #The resulting values are filtered to ensure that 0 < theta < 90 degrees
  thetadeg_mean = 90.      #Mean ejection angle (measured vertically from horizontal, in degrees)
  thetadeg_std  = 45.      #standard deviation of ejection angle, degrees

  #INITIAL VELOCITY (Gaussian distribution)
  vi_mean = 200.                #mean initial velocity (m/s) for vertical blocks
  vi_std  = 100.                #standard deviation of velocity (m/s)

  #BLOCK DIAMETER (LOGNORMAL)
  logdiam_mean = -0.3               #log of mean diameter, meters
  logdiam_std  =  0.5                  #log of std. deviation in diameter

  #BLOCK DIRECTION (degrees east of north; uniform distribution)
  phimin_degrees = 0.
  phimax_degrees = 360.

  #BLOCK TYPE
  #The model randomly chooses between the four types listed in the 'dragtype' variable.
  dragtypes = ['hicube', 'locube', 'sphere', 'shell']
  #To change the drag types considered, modify the elements of dragtypes.
  #For example, to consider only "shell" types, uncomment the line below, and comment out
  #the line above.
  #dragtypes = ['shell', 'shell', 'shell', 'shell']

  ###################################################################################################
  
3)  Run the model by typing "python call_eject3dfun.py" at the command line.
4)  Examine the output files.  The file call_eject3dfun.py should write out a file that looks like the following.  Each row in the table is a different simulation.

  ###################################################################################################
  Output from call_eject3d run on 2020-01-07-1130

     i      vi    diam   theta     phi dragtype    xfinal    yfinal    zfinal  distance    tfinal
           m/s       m     deg     deg               m          m         m        m         s
     0   281.4    1.47    89.3   153.7  locube      42.5     -86.0       0.0      96.0      44.3
     1   352.1    0.14    80.9    77.5   shell    1179.8     261.4       0.0    1208.4      45.5
     2   451.0    0.37    79.4    86.6  sphere    1788.5     106.5       0.0    1791.7      48.1
     3   124.0    0.66    52.1   100.3  locube     988.4    -179.0       0.0    1004.5      17.6
     4   195.2    1.83    59.0   331.4   shell   -1541.0    2828.1       0.0    3220.7      33.5
     5    68.8    0.24    21.8   296.9  sphere    -287.3     145.5       0.0     322.0       5.2
     6   150.2    1.12    61.9   103.9  hicube    1358.2    -335.7       0.0    1399.1      24.5
     7     5.5    1.87     5.0     9.6  sphere       0.1       0.5       0.0       0.5       0.1
     8   185.8    0.98    75.1   119.8  locube     871.0    -499.6       0.0    1004.1      30.4
     9   129.3    1.40    44.5   320.2   shell   -1052.9    1263.2       0.0    1644.4      18.3
  ###################################################################################################
     
The model uses a solution method and block shapes that are described in Mastin (2001).

Reference:
Mastin, L.G., 2001, A simple calculator of ballistic trajectories for blocks ejected during volcanic eruptions, U.S. Geological Survey Open-File Report 01-45, U.S. Geological Survey, p. 26. (https://vhub.org/resources/455)
