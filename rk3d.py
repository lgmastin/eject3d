def rk3d(x,y,z,vx,vy,vz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor):

    import math
    import numpy as np
    from math import sqrt
 
    vv = sqrt(vx*vx + vy*vy + vz*vz) #?? what is this for (unused in script)??
    vwind  = sqrt((vx - xwind)**2 + (vy-ywind)**2 + vz*vz)
    vvxwd  = vx - xwind
    vxwind = vx - xwind
    vvywd  = vy - ywind
    vywind = vy - ywind
    vvx = vx #?? what is this for (unused in script)??
    vvz = vz

    derivs1 = np.zeros((6))
    derivs2 = np.zeros((6))
    derivs3 = np.zeros((6))
    derivs4 = np.zeros((6))

    derivs1=derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor)       

    vwind = sqrt((vx+0.5*derivs1[3])**2 + (vy+0.5*derivs1[4])**2 + (vz+0.5*derivs1[5])**2)
    vvxwd = vx + 0.5 * derivs1[3]
    vvywd = vy + 0.5 * derivs1[4]
    vvz   = vz + 0.5 * derivs1[5]

    derivs2=derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor)       

    vwind = sqrt((vx+0.5*derivs2[3])**2 + (vy+0.5*derivs2[4])**2 + (vz+0.5*derivs2[5])**2)
    vvxwd = vx + 0.5 * derivs2[3]
    vvywd = vy + 0.5 * derivs2[4]
    vvz   = vz + 0.5 * derivs2[5]

    derivs3=derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor)       

    vwind = sqrt((vx+0.5*derivs3[3])**2 + (vy+0.5*derivs3[4])**2 + (vz+0.5*derivs3[5])**2)
    vvxwd = vx + 0.5 * derivs3[3]
    vvywd = vy + 0.5 * derivs3[4]
    vvz   = vz + 0.5 * derivs3[5]

    derivs4=derivs(vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor)       

    r = np.zeros((6,4))
    for i in range(0,6):
        r[i,0] = derivs1[i]
        r[i,1] = derivs2[i]
        r[i,2] = derivs3[i]
        r[i,3] = derivs4[i]

    d=np.zeros((6))
    rknow = np.zeros((6))

    for i in range(0,6):
        d[i] = r[i, 0] / 6 + r[i, 1] / 3 + r[i, 2] / 3 + r[i, 3] / 6
        rknow[0] = x + d[0]                #new x
        rknow[1] = y + d[1]                #new x
        rknow[2] = z + d[2]                #new z
        rknow[3] = vxwind + d[3] + xwind   #new vx
        rknow[4] = vywind + d[4] + ywind   #new vy
        rknow[5] = vz + d[5]               #new vz

    return rknow

def derivs( vwind,vvxwd,vvywd,vvz,xwind,ywind,dt,rhoa,cdnow,area,mass,grav,rhor ):
 
    import math
    import numpy as np
    derivs = np.zeros((6))

    derivs[0] = (vvxwd + xwind) * dt                                                 #distance travelled in x
    derivs[1] = (vvywd + ywind) * dt                                                 #distance travelled in y
    derivs[2] = vvz * dt                                                             #distance travelled in z
    derivs[3] = -rhoa * cdnow * area * vwind ** 2 * vvxwd / (2 * mass * vwind) * dt  #change in vx
    derivs[4] = -rhoa * cdnow * area * vwind ** 2 * vvywd / (2 * mass * vwind) * dt  #change in vy
    derivs[5] = (-grav * (rhor - rhoa) / rhor - rhoa * cdnow * area * vwind ** 2 * vvz /(2 * mass * vwind)) * dt  #in vz

    return derivs
