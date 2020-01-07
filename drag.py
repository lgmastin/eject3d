#Drag functions

def drag_hicube(mach):

    #  Drag of a hi cube
    import sys
    import scipy
    from scipy import ones
    from scipy import interpolate
    import numpy as np

    hicube_mach = np.array([-0.9624766, -0.8649155, -0.6435272, -0.5159475, \
                            -0.3921201, -0.2645403, -0.1932458, -0.1181989, \
                            -0.07317073,-0.03564728, 0.001876173,0.0619137, \
                            0.1332083,   0.2345216,  0.3245779,  0.4333959, \
                            0.4934334,   0.5234522,  0.9362102])
    hicube_cd  = np.array([1.064151,    1.056604,   1.064151,   1.086792,  \
                           1.116981,    1.184906,   1.237736,   1.320755,  \
                           1.396226,    1.486792,   1.607547,   1.728302,  \
                           1.788679,    1.811321,   1.811321,   1.773585,  \
                           1.743396,    1.728302,   1.720755])
    hicube_mach = 10**hicube_mach
    #Create interpolation function
    if mach < hicube_mach[0]:
        drag_hicube = hicube_cd[0]
    elif mach > hicube_mach[18]:
        drag_hicube = hicube_cd[18]
    else:
        #Create interpolation function
        f1 = interpolate.interp1d(hicube_mach,hicube_cd, kind = 'cubic')
        #interpolate using the given mach number
        drag_hicube=f1(mach)

    return drag_hicube
    
def drag_locube(mach):

    #  Drag of a lo cube
    import sys
    import scipy
    from scipy import ones
    from scipy import interpolate
    import math
    import numpy as np

    #Digitized curve of drag coefficient versus Mach number
    locube_mach= np.array([-0.9623352, -0.7702448, -0.6007533, -0.4387947, \
                  -0.3182674, -0.2165725, -0.1412429, -0.08097929, \
                  -0.03954802, 0.01318267, 0.06214689, 0.126177, \
                   0.1939736,   0.2730697,  0.3559322, 0.4274953, \
                   0.5028248, 0.5969868, 0.9397364]);
    locube_cd = np.array([0.7698113, 0.7698113, 0.7773585, 0.7849057, \
                  0.8150944,   0.8528302,  0.8981132,   0.9358491,  \
                  0.9962264,   1.109434,   1.184906,    1.222641,   \
                  1.237736,    1.230189,   1.207547,    1.177359,   \
                  1.14717,     1.139623,   1.139623]);
    locube_mach = 10**locube_mach
    if mach < locube_mach[0]:
        drag_locube = locube_cd[0]
    elif mach > locube_mach[18]:
        drag_locube = locube_cd[18]
    else:
        #Create interpolation function
        f1 = interpolate.interp1d(locube_mach,locube_cd, kind = 'cubic')
        #interpolate using the given mach number
        drag_locube=f1(mach)

    return drag_locube
    
def drag_shell(mach):

    #  Drag of a shell
    import sys
    import scipy
    from scipy import ones
    from scipy import interpolate

    #Digitized curve of drag coefficient versus Mach number
    shell_mach= ([0.024530622, 0.163537451, 0.359782374, \
                  0.556027386, 0.686857275, 0.760449177, \
                  0.866748552, 0.940340275, 0.989401662, \
                  1.071170343, 1.13658493,  1.242884305, \
                  1.373714372, 1.51272104,  1.594490614, \
                  1.741673882, 2.011510618, 2.363116035, \
                  2.706544851, 2.99273568,  3.205334428]);
    shell_cd =  ([0.2307692,   0.2261539,   0.2261539, 0.2261539, \
                  0.2307692,   0.2446154,   0.2861539, 0.3507693, \
                  0.4015385,   0.5123077,   0.5953847, 0.6369231, \
                  0.6461539,   0.6461539,   0.6461539, 0.6276923, \
                  0.6000000,   0.5630770,   0.5307692, 0.5076923, \
                  0.4892308]);
    if mach < shell_mach[0]:
        drag_shell = shell_cd[0]
    elif mach > shell_mach[20]:
        drag_shell = shell_cd[20]
    else:
        #Create interpolation function
        f1 = interpolate.interp1d(shell_mach,shell_cd, kind = 'linear')
        #interpolate using the given mach number
        drag_shell=f1(mach)

    return drag_shell
    
def drag_sphere(reynolds,mach):

    #  Drag of a sphere
    import sys
    import scipy
    from scipy import ones
    from scipy import interpolate
    from math import log10

    if mach<0.5:
    #Calculate Cd as a function of Reynolds number, if M<0.5
      if reynolds<20:
          cd = (24 / reynolds) * (1 + 0.1315 * reynolds ** (0.82 - 0.05 * log10(reynolds)))
      elif reynolds<260:
          cd = (24 / reynolds) * (1 + 0.1935 * reynolds ** 0.6305)
      elif reynolds < 1500:
          cd = 10**(1.6435 - 1.1242 * log10(reynolds) + 0.1558 * log10(reynolds)**2)
      elif reynolds < 12000:
          cd = 10**(-2.4571 + 2.5558 * log10(reynolds) - 0.9295 * log10(reynolds)**2 + 0.1049 * log10(reynolds)**3)
      elif reynolds < 44000:
          cd = 10**(-1.9181 + 0.637 * log10(reynolds) - 0.0636 * log10(reynolds)**2)
      elif reynolds < 338000:
          cd = 10**(-4.339 + 1.5809 * log10(reynolds) - 0.1546 * log10(reynolds)**2)
      elif reynolds < 400000:
          cd = 29.78 - 5.3 * log10(reynolds)
      elif reynolds < 1000000:
          cd = 0.1 * log10(reynolds) - 0.49
      else:
          cd = 0.19 - 80000/reynolds
    else:
      if reynolds < 300000:
          #Calculate Cd as a function of Mach number if Re<3.0e+05
          #using data are digitized from Hoerner, 1965
          losphere_Mach=([0.1000,   0.1278,   0.1714,   0.2382,   0.3151, \
                          0.3728,   0.4140,   0.4897,   0.5364,   0.5834, \
                          0.6257,   0.6710,   0.7298,   0.7665,   0.7938, \
                          0.8163,   0.8336,   0.8513,   0.8755,   0.9455, \
                          0.9792,   1.0141,   1.0876,   1.1996,   1.3606, \
                          1.5112,   1.7505,   2.0135,   2.2679,   2.6454, \
                          3.1293,   3.5246,   4.0542,   4.7957,   5.6333, \
                          9.8610]);
          losphere_cd =([0.4953,   0.4953,   0.4953,   0.4828,   0.4828,   \
                         0.4953,   0.4953,   0.5016,   0.5078,   0.5204,   \
                         0.5266,   0.5455,   0.5956,   0.6395,   0.6771,   \
                         0.6771,   0.6646,   0.6207,   0.6207,   0.6708,   \
                         0.7398,   0.8213,   0.9028,   0.9530,   0.9906,   \
                         1.0094,   1.0219,   1.0282,   1.0157,   0.9906,   \
                         0.9655,   0.9467,   0.9342,   0.9342,   0.9342,   \
                         0.9216]);
          #use cubic spline interpolation; set cd=0.9216 for all M>9.861
          f1 = interpolate.interp1d(losphere_Mach,losphere_cd, kind = 'cubic')
          cd=f1(mach);
      else:
         #Calculate Cd for M>0.5 when Re>3.0e+05
         #using data are digitized from Hoerner, 1965
         hisphere_Mach=([0.0650,   0.1240,   0.1752,   0.2402,   0.3051, \
                         0.3642,   0.4311,   0.5020,   0.5610,   0.6122, \
                         0.6555,   0.6713,   0.6909,   0.7008,   0.7106, \
                         0.7146,   0.7205,   0.7323,   0.7421,   0.7539, \
                         0.7736,   0.7933,   0.8228,   0.8445,   0.8740, \
                         0.8898,   0.9154,   0.9331,   0.9587,   0.9803, \
                         0.9980,   1.0141,   1.0876,   1.1996,   1.3606, \
                         1.5112,   1.7505,   2.0135,   2.2679,   2.6454, \
                         3.1293,   3.5246,   4.0542,   4.7957,   5.6333, \
                         9.8610]);
         hisphere_cd=([ 0.1076,   0.1076,   0.1076,   0.1076,   0.1076, \
                         0.1116,   0.1195,   0.1315,   0.1474,   0.1673, \
                         0.1952,   0.2191,   0.2709,   0.3307,   0.4263, \
                         0.4861,   0.5697,   0.6375,   0.6693,   0.6853, \
                         0.6853,   0.6693,   0.6255,   0.6175,   0.6215, \
                         0.6454,   0.6932,   0.7251,   0.7689,   0.8088, \
                         0.8367,   0.8213,   0.9028,   0.9530,   0.9906, \
                         1.0094,   1.0219,   1.0282,   1.0157,   0.9906, \
                         0.9655,   0.9467,   0.9342,   0.9342,   0.9342, \
                         0.9216]);           
         #interpolate to get cd.  If mach>9.8,cd=0.92
         f2 = interpolate.interp1d(hisphere_Mach,hisphere_cd, kind = 'cubic')
         cd=f2(mach);

    #interpolate using the given mach number
    drag_sphere=cd

    return drag_sphere
    
