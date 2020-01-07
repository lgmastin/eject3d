# This module, called eject3dfun, contains one function, eject3dfunc, that runs ballistic trajectory calculation

def eject3dfunc(inow,diam,vi,thetadeg,phideg,dragtype):

    import sys
    import scipy
    import numpy as np
    import drag
    import rk3d
    import datetime
    from numpy import zeros, nan
    from scipy import sin, cos, sqrt
    from drag import drag_hicube, drag_locube, drag_sphere, drag_shell
    from rk3d import rk3d

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
    # BLOCK 2.  PHYSICAL CONSTANTS 

    pi   = 3.14159                 #pi
    grav = 9.80665                 #gravitational constant, m/s2
    R_air =286.98                  #Specific gas constant for air (J/kg K)

    #############################################################################################
    ## BLOCK 3.  PERFORM INITIAL CALCULATIONS AND ERROR CHECKS

    #dragtype = 'shell'             #can choose 'hicube', 'locube', 'sphere', 'shell', or 'const'
    if dragtype=='const':
        cd_const   = 1              #if dragtype='const', specify value of Cd

    #check whether xi is negative
    if xi > 0:
       xi = input('Please enter a negative value for xi.')

    rad = diam / 2                              #convert block diameter to radius, m
    theta=pi*thetadeg/180                       #convert theta to radians
    phi2 =pi*(0.5-phideg/180)                   #convert phi to radians, degrees CCW from east
    xwind = wind*cos(pi*(0.5-winddir_deg/180))  #x (east) component of wind
    ywind = wind*sin(pi*(0.5-winddir_deg/180))  #y (north) component of wind

    TzeroK = TzeroC + 273.15                    #convert temperature to Kelvin
    ttotcd0 = (-vi * sin(theta) - \
       sqrt(vi ** 2 * (sin(theta) ** 2) - \
       2 * grav * xi))/ (-grav)                  #estimate time of travel, assuming cd=0
    dt=ttotcd0/500                               #Give initial time step as a small fraction of total time

    #############################################################################################
    ## BLOCK 4.  DECLARE OUTPUT ARRAYS

    time=zeros((5000,1)); time[:]=nan 
    x=zeros((5000,1)); x[:]=nan                  #vector for horizontal position from vent (m)
    y=zeros((5000,1)); y[:]=nan                  #vector for horizontal position from vent (m)
    z=zeros((5000,1)); z[:]=nan                  #vector for vertical position above vent (m)
    vx=zeros((5000,1)); vx[:]=nan                #vector for x-(east) components of velocity (m/s)
    vy=zeros((5000,1)); vy[:]=nan                #vector for y-(north) components of velocity (m/s)
    vz=zeros((5000,1)); vz[:]=nan                #vector for z-components of velocity (m/s)
    mach=zeros((5000,1)); mach[:]=nan            #mach number now
    temp=zeros((5000,1)); temp[:]=nan            #air temperature now
    rhoa=zeros((5000,1)); rhoa[:]=nan            #air density now
    reynolds = zeros((5000,1)); reynolds[:]=nan  #Reynolds number now
    cdnow=zeros((5000,1)); cdnow[:]=nan          #cd now

    #Provide values to first elements in arrays
    x[0] = 0.                                    #initialize x
    y[0] = 0.                                    #initialize y
    z[0] = 0.                                    #initialize z
    vx[0]=vi*cos(theta)*cos(phi2)                #initial vx
    vy[0]=vi*cos(theta)*sin(phi2)                #initial vy
    vz[0]=vi*sin(theta)                          #initial vz
    #v=vi                                         #initial velocity, m/s
    time[0]=0                                    #initial time (seconds)
    zmax=0                                       #maximum point in block's trajectory
    po=101300 * ((TzeroK - elev * lapse / 1000)/ TzeroK) ** \
        (-grav / (R_air ** 2 * lapse / 1000))  #air pressure at vent

    #############################################################################################
    ## BLOCK 5.  CALCULATE BLOCK MASS AND FRONTAL AREA DEPENDING ON DRAG TYPES

    if dragtype=='hicube':
        mass = 8 * rad ** 3 * rhor  #mass of clast calculated based on the option "high cube"
        area = 4 * rad * rad       #area of clast calculated based on the option "high cube"
    elif dragtype=='locube':
        mass = 8 * rad **3 * rhor   #mass of a cube
        area = 4 * rad **2 * sqrt(3)   #frontal area of a cube pointing forward
    elif dragtype=='sphere':
        mass = (4/3)*pi*rad**3*rhor   #mass of a sphere
        area = pi*rad**2              #frontal area of a sphere
    elif dragtype=='shell':
        mass = pi * 2.4107 * rad**3 * rhor   #mass of a G1 shaped ogive cylinder
        area = pi * rad**2                   #frontal area of an ogive cylinder
    elif dragtype=='const':
        mass = (4/3)*pi*rad**3*rhor   #assume the mass of a sphere
        area = pi*rad**2              #and the frontal area of a sphere

    #############################################################################################
    ## BLOCK 6.  START INTEGRATING THROUGH TIME

    xnow=x[0]
    ynow=y[0]
    znow=z[0]
    vxnow=vx[0]
    vynow=vy[0]
    vznow=vz[0]
    i=0

    while znow >= xi: #run while block elevation is above the ground elevation (xi) 
    
        #Specify current properties
        x[i]=xnow 
        y[i]=ynow 
        z[i]=znow
        vx[i]=vxnow
        vy[i]=vynow
        vz[i]=vznow
        temp[i] = TzeroK - (znow + elev) * lapse / 1000              #temperature at znow
        visc = (0.000172 * (390 / (temp[i] + 117)) * \
            (temp[i] / 273) ** 1.5) / 10                         #air viscosity at znow
        pressure = po *((TzeroK - (znow + elev) * \
            lapse / 1000) / TzeroK) ** \
            (grav / (R_air * lapse / 1000))                      #air pressure at znow
        rhoa[i] = pressure / (R_air * temp[i])                       #air density at znow
        c_sound    = 20.116 * sqrt(temp[i])                          #sound speed at znow
        vwind   = sqrt((vxnow-xwind)**2 + (vynow-ywind)**2 \
                + vznow**2)                                      #velocity minus wind component
        reynolds[i] = rhoa[i] * vwind * 2 * rad / visc               #Reynolds number
        mach[i] = vwind / c_sound                                    #Mach number

        #Calculate drag coefficient
        if dragtype=='hicube':
            cd = drag_hicube(mach[i])
        elif dragtype=='locube':
            cd = drag_locube(mach[i])
        elif dragtype=='sphere':
            cd = drag_sphere(reynolds[i],mach[i])
        elif dragtype=='shell':
             cd = drag_shell(mach[i])
        elif dragtype=='const':
            cd=cd_const

        #calculate drag reduction in reduced-drag zone
        if dragred>0 and sqrt(xnow**2 + znow**2) < dragred:
            cdnow[i] = cd * (sqrt(xnow**2 + znow**2)/dragred)**2
        else:
            cdnow[i]=cd

        #############################################################################################
        ## Runge Kutta to determine new position and velocities
           
        RK = zeros((4,1))

        RK=rk3d(xnow,ynow,znow,vxnow,vynow,vznow,xwind,ywind,dt,rhoa[i],cdnow[i],area,mass,grav,rhor)

        #update position and velocity
        xnow  = RK[0]
        ynow  = RK[1]
        znow  = RK[2]
        vxnow = RK[3]
        vynow = RK[4]
        vznow = RK[5]
                                       
        #update i so that next value is stored in the next position of storage vector
        i=i+1

        #update time
        time[i] = time[i-1] + dt

        #print('%4d    %4.2f    %6.1f    %5.1f    %6.2f    %5.2f' % (i,time[i],xnow,znow,vxnow,vznow))
                                                                   

        #END OF ITERATIONS

    ##################################################################################################
    #BLOCK 7:  INTERPOLATE TO FIND THE FINAL VALUES

    xfinal = x[i-1] + (xnow - x[i-1]) * (xi - z[i-1]) / (znow - z[i-1])           #final x
    yfinal = y[i-1] + (ynow - y[i-1]) * (xi - z[i-1]) / (znow - z[i-1])           #final y
    zfinal = xi                                                                   #final z
    tfinal = time[i-1] + dt * (xi - z[i-1]) / (znow - z[i-1])                          #final time (s)
    vxfinal = vx[i-1] + (vxnow - vx[i-1]) * (xi - z[i-1]) / (znow - z[i-1])       #final vx
    vyfinal = vy[i-1] + (vynow - vy[i-1]) * (xi - z[i-1]) / (znow - z[i-1])       #final vy
    vzfinal = vz[i-1] + (vznow - vz[i-1]) * (xi - z[i-1]) / (znow - z[i-1])       #final vz
    vfinal = sqrt(vxfinal ** 2 + vzfinal ** 2)                                    #final velocity

    #Make these the final values of the arrays
    x[i]    = xfinal 
    y[i]    = yfinal 
    z[i]    = zfinal 
    vx[i]   = vxfinal 
    vy[i]   = vyfinal 
    vz[i]   = vzfinal 
    time[i] = tfinal

    #convert to list values, not arrays
    xfinal  = xfinal.item()
    yfinal  = yfinal.item()
    vxfinal = vxfinal.item()
    vyfinal = vyfinal.item()
    vzfinal = vzfinal.item()
    tfinal  = tfinal.item()

    #Calculate other properties at the final elevation
    distance = sqrt(xfinal**2 + yfinal**2)
    temp[i] = TzeroK - (zfinal + elev) * lapse / 1000              #temperature at zfinal
    visc = (0.000172 * (390 / (temp[i] + 117)) * \
            (temp[i] / 273) ** 1.5) / 10                     #air viscosity at zfinal
    pressure = po *((TzeroK - (zfinal + elev) * \
         lapse / 1000) / TzeroK) ** \
         (grav / (R_air * lapse / 1000))                      #air pressure at zfinal
    rhoa[i] = pressure / (R_air * temp[i])                        #air density at zfinal
    c_sound    = 20.116 * sqrt(temp[i])                          #sound speed at zfinal
    velocity = sqrt(vxnow*vxnow + vynow*vynow + vznow*vznow)
    vwind   = sqrt((vxnow-xwind)**2 + (vynow-ywind)**2 \
                + vznow**2)                                  #velocity minus wind component
    reynolds[i] = rhoa[i] * vwind * 2 * rad / visc               #Reynolds number
    mach[i] = vwind / c_sound                                    #Mach number
    cdnow[i]=cdnow[i-1]                                          #Cd (using the last value is good enough)

    #Trim arrays
    ifinal = i;
    x=x[0:ifinal+1]; 
    y=y[0:ifinal+1]; 
    z=z[0:ifinal+1]; 
    vx=vx[0:ifinal+1];
    vy=vy[0:ifinal+1];
    vz=vz[0:ifinal+1];
    time=time[0:ifinal+1]
    mach=mach[0:ifinal+1]
    cdnow=cdnow[0:ifinal+1]
    temp=temp[0:ifinal+1]
    rhoa=rhoa[0:ifinal+1]
    reynolds=reynolds[0:ifinal+1]

    ##################################################################################################
    #WRITE OUT RESULT

    if write_outfile == 'yes':
        #Get current date and time for filename
        now = datetime.datetime.now()
        datestring = now.strftime('%Y-%m-%d-%H%M')
        filename = "{0:s}_run{1:03d}.txt".format(filename_stem,inow)

        outfile = open(filename, 'w')

        print('Output %s from Eject, Python version, run on %s\n\n' % (filename,datestring), file=outfile)
        print('Block parameters:', file=outfile)
        print('                       Block diameter, m:  %5.2f' % (diam), file=outfile)
        print('                    Block density, kg/m3:   %4.0f' % (rhor), file=outfile)
        print('                             Block shape:  %s' % (dragtype), file=outfile)
        print('                   Initial velocity, m/s: %6.1f' % (vi), file=outfile)
        print(' Ejection angle from horizontal, degrees:   %4.1f' % (thetadeg), file=outfile)
        print('      Ejection direction, degrees E of N:  %5.1f' % (phideg), file=outfile)
        print('                         Wind speed, m/s:   %4.1f' % (wind), file=outfile)
        print('             Wind direction, deg. E of N:   %4.1f' % (winddir_deg), file=outfile)
        print('landing point meters below takeoff point:   %4.0f' % (xi), file=outfile)
        print('     extent of reduced drag zone, meters:   %4.1f\n' % (dragred), file=outfile)

        print('Meteorologic conditions:', file=outfile)
        print('            Temperature at sea level, Celsius:  %4.1f' % (TzeroC), file=outfile)
        print('            Thermal lapse rate, deg. C per km: %5.2f' % (lapse), file=outfile)
        print('Elevation of takeoff point above sea level, m:  %4.0f\n' % (elev), file=outfile)

        print('TRAJECTORY CALCULATIONS', file=outfile)
        print('###################################################################################################################################', file=outfile)
        print('    i   t (s)     x (m)     y (m)     z (m)   vx(m/s)  vy (m/s)  vz (m/s)   Mach #     Reynolds #     Cd   Temp (K) rho_air (kg/m3)', file=outfile)

        for i in range(0,ifinal+1):
            print('%5d  %6.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %8.3f %14.3f %6.3f %10.2f %10.3f' % \
               (i,time[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],mach[i],reynolds[i],cdnow[i],temp[i],rhoa[i]), \
               file=outfile)

        print('###################################################################################################################################', file=outfile)
        print('Range=%7.1f meters, flight time=%7.1f seconds, final location (x,y,z)=(%7.1f,%7.1f,%7.1f)' % (distance,tfinal,xfinal,yfinal,zfinal), file=outfile)

    return [xfinal, yfinal, zfinal, tfinal]

    #outfile.close()

