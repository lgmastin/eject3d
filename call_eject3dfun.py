#  This script, written by Larry G. Mastin (USGS, lgmastin@usgs.gov) generates a large number 
#  of inputs for ballistic block ejection and then calls the python function "eject3dfunc", which
#  calculates block trajectories and returns the final time of flight and landing location of each block.

#  The inputs  that define the range and distribution of block size, initial velocity, and
#  trajectory angles are specified in Block 1 of this script

#  For a specific volcano, some other parameters may have to be set in the module eject3dfun: for example,
    # elev         = elevation of vent, in meters above sea level
    # lapse        = thermal lapse rate in the atmosphere, in degrees C per kilometer
    # wind         = the wind speed in meters per second
    # winddir_deg  = the direction toward which the wind is blowing, in degrees east of north
    # TzeroC       = the temperature at sea level, in Celsius
    # rhor         = block density, kg/m3
    # dragred      = extent of zone of reduced drag, in meters around the vent.
    # write_outfile = 'yes' if a trajectory output file is to be written for every block ejection.
    # filename_stem = stem of file names to be written out

###################################################################################################
#Python modules required for this script
import numpy as np
import datetime
import scipy
from scipy import random
import sys
import math
import eject3dfun
from eject3dfun import eject3dfunc
from math import sqrt, sin, cos

###################################################################################################
#  BLOCK 1.  PARAMETERS THAT DEFINE THE RANGE AND DISTRIBUTION OF INITIAL BLOCK SIZES, SHAPES,
#            EJECTION VELOCITIES, AND TRAJECTORY ANGLES

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
# the line above.
#dragtypes = ['shell', 'shell', 'shell', 'shell']

###################################################################################################
#BLOCK 2.  SET PHYSICAL CONSTANTS
pi         =    3.14159

###################################################################################################
#BLOCK 3.  CALCULATE RANDOM DISTRIBUTIONS

#Vertical angle
thetadeg  = random.normal(loc = thetadeg_mean, scale=thetadeg_std, size=n_ejects)
for i in range(0,n_ejects):
    thetadeg[i] = max(5.,thetadeg[i])          #set minimum angle to 5 degrees
    if thetadeg[i] > 90.:
        thetadeg[i] = thetadeg[i] - 90.        #Make sure everything is between zero and 90
theta = pi*thetadeg/180.                       #convert to radians

#Initial velocity
vi =  random.normal(loc=vi_mean, scale=vi_std, size=n_ejects)    #Generate initial velocities
for i in range(0,n_ejects):
    vi[i] = max(2.,vi[i])               #make sure nothing has zero or negative ejection velocity
    vi[i] =  vi[i] * sin(theta[i])**(3/2)   #account for higher velocities at steeper angles

#Block diameter (set a lognormal distribution)
logdiam   = random.normal(loc=logdiam_mean, scale=logdiam_std, size=n_ejects)
diam      = 10**logdiam

#Initial horizontal direction
phideg  = phimin_degrees + (phimax_degrees-phimin_degrees)*random.random(n_ejects)

#Block type
typenow   = 4.*random.random(n_ejects)                 #generate random numbers between 0 and 4
typenow   = typenow.astype(int)                        #convert them to integers

###################################################################################################
#BLOCK 4.  DECLARE OUTPUT VARIABLES AND CALL EJECT3DFUNC

#Declare output variables
xfinal = np.zeros((n_ejects,1))
yfinal = np.zeros((n_ejects,1))
zfinal = np.zeros((n_ejects,1))
tfinal = np.zeros((n_ejects,1))

#Open output file, Write table header
now = datetime.datetime.now()
datestring = now.strftime('%Y-%m-%d-%H%M')
filename ="{0:s}{1:s}.txt".format(filename_stem,datestring)
outfile = open(filename,'w')                #open output file
print('Output from call_eject3d run on %s\n' % (datestring), file=outfile)
print('   i      vi    diam   theta     phi dragtype    xfinal    yfinal    zfinal  distance    tfinal')
print('         m/s       m     deg     deg               m          m         m        m         s')
print('   i      vi    diam   theta     phi dragtype    xfinal    yfinal    zfinal  distance    tfinal', file=outfile)
print('         m/s       m     deg     deg               m          m         m        m         s', file=outfile)

#Call eject3dfunc
for i in range(0,n_ejects):
    dragtypenow = dragtypes[typenow[i]]                                   #set dragtypenow
    [xfinal[i], yfinal[i], zfinal[i], tfinal[i]] = \
            eject3dfunc(i,diam[i],vi[i],thetadeg[i],phideg[i],dragtypenow)
    distance = sqrt(xfinal[i]**2 + yfinal[i]**2)
    print('%4d%8.1f%8.2f%8.1f%8.1f%8s%10.1f%10.1f%10.1f%10.1f%10.1f' \
            % (i,vi[i],diam[i],thetadeg[i],phideg[i],dragtypenow, \
            xfinal[i],yfinal[i],zfinal[i],distance,tfinal[i]))
    print('%4d%8.1f%8.2f%8.1f%8.1f%8s%10.1f%10.1f%10.1f%10.1f%10.1f' \
            % (i,vi[i],diam[i],thetadeg[i],phideg[i],dragtypenow, \
            xfinal[i],yfinal[i],zfinal[i],distance,tfinal[i]), file=outfile)

outfile.close()
print('All done')
