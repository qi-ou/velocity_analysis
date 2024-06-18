#!/usr/bin/env python3
'''
## lab2_lib.py

Library of python functions to be used with lab 2.

'''

# packages
import numpy as np
import scipy.special as sc
import scipy.stats as st
import matplotlib.pyplot as plt


def psvelo2rot(in_file, sig_max, std_scale):
    '''
    Estimates plate angular rotation and associated statistics
    from psvelo input file.

     Usage:
       OM, COM, Eu, EL, chi2, chi2r, dof, rmsNEH, R = psvelo2rot(in_file, sig_max, std_scale)

    Input:
        in_file = input velocity file in psvelo format (mm/yr)
        sig_max = minimum velocity uncertainty (mm/yr)
        std_scale = std. dev. scaling factor
    Output:
        OM = Wx  Wy   Wz     (deg/My)
        COM = [XX  XY  XZ
              XY  YY  YZ
              XZ  YZ  ZZ]   (deg/My)**2
        Eu = Lat (deg)  Lon (deg)  Ang (deg/My) - Euler pole
        EL = Semi-max  Semi-min   Azim (deg)  Ang (deg/My)
        chi2, chi2, dof
        rmsNEH  rms N E H  (mm/yr)
                wrms N E H  (mm/yr)
        R = residuals (lon, lat, ve, vn, sve, svn, cne)
    '''
    
    # Read input data file
    lon, lat, ve, vn, sve, svn, cne = np.loadtxt(in_file, dtype={'names': ('long', 'lat', 'N', 'E', 'N_std', 'E_std', 'correlation'), 'formats': (np.float, np.float, np.float, np.float, np.float, np.float, np.float)},unpack=True)
    
    # Verticals
    h = np.zeros_like(lon)
    vu = np.zeros_like(lat)
    svu = np.zeros_like(lat)
    cnu = np.zeros_like(cne)
    ceu = np.zeros_like(cne)
    nsites = lat.size;
    
    print(nsites,'sites provided')
    
    # Select sites
    I = np.where(np.sqrt(sve**2 + svn**2) < sig_max)[0]
    lon = lon[I]; lat = lat[I]; h = h[I]
    ve = ve[I]; vn = vn[I]; vu = vu[I]
    sve = sve[I]; svn = svn[I]; svu = svu[I]
    cne = cne[I]; cnu = cnu[I]; ceu = ceu[I]
    #site = site[I]
    nsites = lat.size
    
    print('Selected ' + str(nsites) + ' sites with uncertainty less than ' + str(sig_max) + ' mm/yr');
    
    # Rescale variances
    sve = sve * std_scale
    svn = svn * std_scale
    svu = svu * std_scale
    
    # Convert site positions and velocities to cartesian and m/yr
    print('Converting positions and velocities to cartesian')
    x, y, z = wgs2xyz(lon, lat, h)
    Pxyz = np.vstack((x, y, z)).T
    O = np.vstack((lat, lon, h)).T
        
    # Make sure velocities are converted to m/yr
    Vneu = np.vstack((vn, ve, vu)).T / 1000.0
    SVneu = np.vstack((svn, sve, svu)).T / 1000.0
    CORneu = np.stack([cne, cnu, ceu]).T # rjw / (1000.0^2)
    Vxyz, CVxyz = neu2xyz(O, Vneu, SVneu, CORneu)
    
    # Compute angular rotation and statistics
    print('Computing angular velocity and statistics')
    OM, COM, Eu, EL, STAT, RES = vel2rot(Pxyz, Vxyz, CVxyz)
    COM = np.array([COM[0,0], COM[0,1], COM[0,2], COM[1,1], COM[1,2], COM[2,2]]) # rjw .*(1000.0^2);
    chi2 = STAT[0]
    chi2r = STAT[1]
    dof = STAT[2]

    # Extract residuals and covariance matrix of residuals
    Vres = RES[:,0:3]
    CVres = np.vstack((RES[:,3], RES[:,4], RES[:,5], RES[:,6], RES[:,7], RES[:,8])).T
    
    # Sum variance of observations + variance of residuals
    Ctot = CVres # rjw +CVxyz already been summed in vel2rot, no need to do again
    SV = np.sqrt(np.vstack((Ctot[:,0], Ctot[:,3], Ctot[:,5])).T)
    COR = np.vstack((np.array([Ctot[:,1] / (Ctot[:,0]**0.5 * Ctot[:,3]**0.5)]), 
                     np.array([Ctot[:,2] / (Ctot[:,0]**0.5 * Ctot[:,5]**0.5)]), 
                     np.array([Ctot[:,4] / (Ctot[:,3]**0.5 * Ctot[:,5]**0.5)]))).T # rjw
    
    # Convert residuals to NEU
    VRneu, CRneu, LLH = xyz2neu(Pxyz, Vres, SV, COR)
    VRneu = VRneu * 1000
    CRneu = CRneu * (1000**2)
    SVRneu = np.sqrt(np.vstack((CRneu[:,0], CRneu[:,3], CRneu[:,5])).T)
    R = np.vstack((LLH[:,1], LLH[:,2], VRneu[:,1], VRneu[:,0], SVRneu[:,1], SVRneu[:,0], CRneu[:,1]/(SVRneu[:,1]*SVRneu[:,0]))).T
    
    # Compute weighted mean of residuals
    wmean_N, wrms_N = get_wrms(VRneu[:,1], SVRneu[:,1])
    wmean_E, wrms_E = get_wrms(VRneu[:,0], SVRneu[:,0])
    wmean_H, wrms_H = get_wrms(np.vstack((VRneu[:,0], VRneu[:,1])).T, np.vstack((SVRneu[:,0], SVRneu[:,1])).T);
    
    # Compute rms %%% rjw
    rms_N = np.sqrt(np.mean(VRneu[:,1]**2))
    rms_E = np.sqrt(np.mean(VRneu[:,0]**2))
    rms_H = np.sqrt(np.mean(np.vstack((VRneu[:,0], VRneu[:,1]**2))))
    
    rmsNEH = np.array([[rms_N, rms_E, rms_H], [wrms_N, wrms_E, wrms_H]])
    
    # Display outputs
    print('-------------------------')
    print('ANGULAR VELOCITY:      Wx          Wy          Wz     (deg/My)')
    print('                    %7.4f     %7.4f     %7.4f' % tuple(OM))
    print('COVARIANCE ELMTS:    XX      XY      XZ      YY      YZ      ZZ   (deg/My)**2')
    print('                  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f' % tuple(COM))
    print('ROTATION PARAMETERS:  Lon (deg)   Lat (deg)   Ang (deg/My)')
    print('                      %7.4f     %7.4f     %7.4f' % tuple(Eu))
    print('STANDARD ERROR ELLIPSE: Semi-max  Semi-min   Azim (deg)  Ang (deg/My)')
    print('                        %7.4f   %7.4f    %7.4f   %7.4f' % tuple(EL))
    print('STATISTICS:')
    print('  CHI**2            : %.2f' % chi2)
    print('  DEGREES OF FREEDOM: %d' % dof)
    print('  REDUCED CHI**2    : %.2f' % chi2r)
    print('\nRESIDUALS (mm/yr):')
    print('   Lon       Lat    Ve    Vn    Sve   Svn    Corr   Site')                    
    for i in range(nsites):
        print('%8.3f %8.3f %5.2f %5.2f %5.2f %5.2f %8.5f  ' % tuple(R[i,:]))
    print('\nWEIGHTED MEAN OF RESIDUALS:')
    print('    North = %.1f +- %.1f mm/yr' % (wmean_N, wrms_N))
    print('    East  = %.1f +- %.1f mm/yr' % (wmean_E, wrms_E))
    print('    Horiz = %.1f +- %.1f mm/yr' % (wmean_H, wrms_H))

    # Rewrite covariance matrix
    COM = np.array([[COM[0], COM[1], COM[2]],
        [COM[1], COM[3], COM[4]],
        [COM[2], COM[4], COM[5]]])               
    
    
    return OM, COM, Eu, EL, chi2, chi2r, dof, rmsNEH, R

#-------------------------------------------------------------------------------

def wgs2xyz(lam, phi, h):
    '''
    Converts lam(longitude) phi(latitude) ellipsoidal coordinates
    from WGS-84 to ECEF cartesian coordinates
    lam, phi, h can be vectors

    Call: x, y, z = wgs2xyz(lam, phi, h)
        lon, lat in decimal degrees
        h in meters above ellipsoid
    '''
    
    # Semimajor and semiminor axis for WGS-84
    a = 6378137.0
    b = 6356752.314
    f = 1.0/298.257222101
    ee = 2*f - f**2
    
    # Degrees to radians
    lam = lam*np.pi/180
    phi = phi*np.pi/180

    # Radius of curvature in prime vertical
    N = a / np.sqrt(1-(np.sin(phi))**2*ee)

    x = np.cos(phi)*np.cos(lam)*(N+h)
    y = np.cos(phi)*np.sin(lam)*(N+h)
    z = np.sin(phi)*(N*(b**2/a**2) + h)
      
    return x, y, z

#-------------------------------------------------------------------------------

def xyz2wgs(S):
    '''
    Converts cartesian coordinates (x,y,z) into
    ellipsoidal coordinates (lat,lon,alt) on WGS-84
    according to a non iterative method (Bowring 76,
    see also GPS Theory and Practice, p. 258).

    Input:
      S = nx4 matrix with time, X, Y, Z
          A! first column of S is time but can be dummy.
    Output:
      R = nx4 matrix with time, lon (lam), lat (phi), elevation
           (lon,lat in degrees!)

     Usage:
       R = xyz2wgs(S)
    '''
    
    # Make sure S is a numpy matrix
    #S = np.matrix(S)
    
    # WGS-84 PARAMETERS
    # semimajor and semiminor axis
    a = 6378137.0
    b = 6356752.314
    # flattening
    f = 1.0/298.257222101
    # eccentricity
    eo = 2*f - f**2

    # Second numerical eccentricity
    e1 = (a**2 - b**2) / b**2

    # Read data
    t = S[:,0]
    x = S[:,1]
    y = S[:,2]
    z = S[:,3]
    
    # Auxiliary quantities
    p = np.sqrt(np.power(x,2) + np.power(y,2))
    theta = np.arctan2(z*a, p*b)
    
    # Longitude
    lam = np.arctan2(y, x)

    # Latitude
    phi = np.arctan2(z + np.power(np.sin(theta),3)*e1*b, p - np.power(np.cos(theta),3)*eo*a)

    # Radius of curvature in prime vertical
    N = a**2 / np.sqrt(np.power(np.cos(phi),2)*a**2 + np.power(np.sin(phi),2)*b**2); #this is matrix divide in the matlab script?

    # Geocentric (?) altitude
    alt_g = (p / np.cos(phi)) - N

    # Ellipsoidal altitude
    alt = np.multiply(p,np.cos(phi)) + np.multiply(z,np.sin(phi)) - a*np.sqrt(1.0 - eo*np.power(np.sin(phi),2))
    
    # Fill out result matrix
    #R = np.array((t, lam*(180.0/np.pi), phi*(180.0/np.pi), alt))
    R = np.vstack((t, lam*(180.0/np.pi), phi*(180.0/np.pi), alt)).T

    return R

#-------------------------------------------------------------------------------

def neu2xyz(O, V, SV, COR):
    '''
    Convert local topocentric into ECEF

    Input:
        O = origin vector in ellipsoidal coordinates (lat lon height)
        V = position or velocity vector in NEU frame (m or m/yr)
        SV = stdev in NEU frame (m or m/yr)
        COR = correlations, NE NU EU
        (NOTE: O, V, SV, COR can be n x 3 matrices, n = # of sites)

    Output:
        XYZ = output in ECEF Cartesian frame (m)
        CXYZ = associated covariance (m^2), format is:
            Cxx Cxy Cxz Cyy Cyz Czz
        (NOTE: XYZ and CXYZ will be matrices with n rows)

    Call: XYZ, CXYZ = neu2xyz(O, V, SV, COR)
    '''
    
    # If O is a single point, make it the same size as V
    if (O.shape[0] == 1):
        lat = np.ones(V.size) * O[0]
        lon = np.ones(V.size) * O[1]
        h = np.ones(V.size) * O[2]
    else:
        lat = O[:,0]; lon = O[:,1]; h = O[:,2]

    # Read rest of input
    vn = V[:,0]; ve = V[:,1]; vu = V[:,2]
    svn = SV[:,0]; sve = SV[:,1]; svu = SV[:,2]
    cne = COR[:,0]; cnu = COR[:,1]; ceu = COR[:,2]

    # Convert position(s) of origin to ECEF
    XR, YR, ZR = wgs2xyz(lon,lat,h)

    # Compute sines and cosines
    cp = np.cos(lon * np.pi/180); sp = np.sin(lon * np.pi/180); # longitude = phi
    cl = np.cos(lat * np.pi/180); sl = np.sin(lat * np.pi/180); # latitude = lam

    # For each site
    XYZ = np.array([])
    CXYZ = np.array([])
    for i in range(V.shape[0]):
      
        # Build the rotation matrix
        R = np.vstack((np.array([-sl[i]*cp[i], -sl[i]*sp[i], cl[i]]),
               np.array([-sp[i], cp[i], 0]),
               np.array([cl[i]*cp[i], cl[i]*sp[i], sl[i]])))
        
        # Apply the rotation
        XYZi = np.dot(R.T,np.array([[vn[i]], [ve[i]], [vu[i]]])).reshape(3) # Force back to array

        # svu cannot be zero or R'*CVi*R may return negative variances
        # if (svu(i)==0) svu(i)= mean([svn(i) sve(i)]); end;

        # Build covariance for that site
        CVi = np.vstack((
            np.array([svn[i]**2, cne[i]*svn[i]*sve[i],   cnu[i]*svn[i]*svu[i]]),
            np.array([cne[i]*svn[i]*sve[i],   sve[i]**2, ceu[i]*sve[i]*svu[i]]),
            np.array([cnu[i]*svn[i]*svu[i],   ceu[i]*sve[i]*svu[i],   svu[i]**2])))
                
        # Propagate covariance
        CXYZi = np.dot(np.dot(R.T, CVi), R)

        # Increment result matrices
        if i!=0:
            XYZ = np.vstack((XYZ,
                   np.array([XYZi[0], XYZi[1], XYZi[2]])))
            CXYZ = np.vstack((CXYZ,
                   np.array([CXYZi[0,0], CXYZi[0,1], CXYZi[0,2], CXYZi[1,1], CXYZi[1,2], CXYZi[2,2]])))
        else:
            XYZ = np.array([XYZi[0], XYZi[1], XYZi[2]]).T
            CXYZ = np.array([CXYZi[0,0], CXYZi[0,1], CXYZi[0,2], CXYZi[1,1], CXYZi[1,2], CXYZi[2,2]])
            
        #print(XYZ); print(' ')
    
    return XYZ, CXYZ


#-------------------------------------------------------------------------------

def vel2rot(P, V, CV):
    '''
    Estimate rotation parameters from site velocities.

    Input:
      - P = site positions (m) in Cartesian frame
      - V = site velocities (m/yr) in Cartesian frame
      - CV = velocity covariance in Cartesian frame, format is:
                        Cxx Cxy Cxz Cyy Cyz Czz  (m/yr)^2
      [use xyz2neu to convert velocities from NE to cartesian]
    Output:
      - OM = rotation vector (deg/Myr)
      - COM = associated covariance matrix (deg^2/Myr^2)
      - E = Euler parameters [lon,lat,ang] (deg/Myr)
      - EL = standard errors on Euler parameters:
        smaj = semi-major axis
        smin = semi-minor axis
        az = direction of semi-major axis
        sang = angular velocity
      - STAT = statistics:
        chi2 = chi-square
        chi2r = reduced chi-square
        dof = degres of freedom
      - RES = residuals:
            Vx Vy Vz Cxx Cxy Cxz Cyy Cyz Czz

    Call: OM, COM, Eu, EL, STAT, RES = vel2rot(P, V, CV)
    '''
    
    # Read site positions and convert to unit vector
    X = P[:,0]; Y = P[:,1]; Z = P[:,2]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    Xu = X/R; Yu = Y/R; Zu = Z/R

    # Read velocities and covariance
    VX = V[:,0]; VY = V[:,1]; VZ = V[:,2]
    CVxx = CV[:,0]; CVxy = CV[:,1]; CVxz = CV[:,2]
    CVyy = CV[:,3]; CVyz = CV[:,4]; CVzz = CV[:,5]
    
    # Number of data points
    n = Xu.size

    # Built observation vector, design matrix, and covariance
    L = np.array([])
    C = np.zeros((3*n,3*n))    

    for i in range(n):
        # Observation vector
        L = np.hstack((L, VX[i], VY[i], VZ[i])).T

        # Design matrix
        if i==0:
            A = np.vstack((np.array([0.0,      Zu[i],  -Yu[i]]),
                np.array([-Zu[i],     0.0,   Xu[i]]),
                 np.array([Yu[i],  -Xu[i],     0.0])))
        else:
            A = np.vstack((A,
                 np.array([0.0,      Zu[i],  -Yu[i]]),
                np.array([-Zu[i],     0.0,   Xu[i]]),
                 np.array([Yu[i],  -Xu[i],     0.0])))

        # Covariance matrix
        C[3*i:3*i+3,3*i:3*i+3] = np.array([[CVxx[i], CVxy[i], CVxz[i]],
                                  [CVxy[i], CVyy[i], CVyz[i]],
                                  [CVxz[i], CVyz[i], CVzz[i]]])
    
        
    # Test C for singularity
    if (1/np.linalg.cond(C) < 1e-3):
        print('WARNING: covariance matrix ill-conditioned in vel2rot')
        print('         removing off-diagonal terms')
        C = np.diag(np.diag(C))

    # Solve for OM -- result will differ according to solver...
    C = np.matrix(C); L = np.matrix(L).T

    OM = np.linalg.inv(A.T*(C.I)*A) * A.T * C.I * L
    
    # Compute covariance of unknowns
    # A! PROBLEM HERE -- C HAS ALWAYS POSITIVE DIAGONAL TERMS, BUT
    # COM SOMETIMES HAS NEGATIVE DIAGONAL TERMS... DON'T KNOW WHY.
    COM = np.linalg.inv(A.T*(C.I)*A)

    # Convert OM from m/y to deg/Myr
    Re = 6378137.0
    con = (Re * np.pi/180) * 1e-6
    OM = OM / con
    COM = COM / (con**2)
    
    # Make sure OM and COM and ndarrays not matrices
    OM = np.array(OM).flatten(); COM = np.array(COM)
    
    # Convert rotation vector into Euler pole and ang. velocity
    Eu, CE, EL = rot2eul(OM, COM)

    # Compute residuals and associated variance
    Vxyz, Venu = rotate(OM, COM, P)
    VX_res = VX - Vxyz[:,0]
    VY_res = VY - Vxyz[:,1]
    VZ_res = VZ - Vxyz[:,2]
    CV_res = CV + Vxyz[:,3:9]
    
    RES = np.vstack((VX_res, VY_res, VZ_res, CV_res.T)).T

    # Compute stats
    Vobs = np.hstack((VX, VY, VZ))
    Vmod = np.hstack((Vxyz[:,0], Vxyz[:,1], Vxyz[:,2]))
    Sigm = np.hstack((np.sqrt(CVxx), np.sqrt(CVyy), np.sqrt(CVzz)))
    chi2 = np.sum(((Vobs-Vmod)**2)/Sigm**2)
    dof = 3*n-3
    chi2r = chi2/dof
    STAT = np.array([chi2, chi2r, dof])
    
    return OM, COM, Eu, EL, STAT, RES


#-------------------------------------------------------------------------------

def xyz2neu(O, V, SV, COR):
    '''
    Convert ECEF into local topocentric

    Input:
       O = origin vector in ECEF frame (m)
       V = position or velocity vector in ECEF frame (m or m/yr)
       SV = stdev in ECEF frame (m or m/yr)
       COR = correlations, XY XZ YZ
       (NOTE: O, V, SV, COR can be n x 3 matrices, n = # of sites)

    Output:
       NEU = output on NEU frame (m)
       NEU = associated covariance (m^2), format is:
                        Cnn Cne Cnu Cee Ceu Cuu
                 (NOTE: NEU and CNEU will be matrices with n rows)

    Call: NEU, CNEU, E = xyz2neu(O, V, SV, COR)
    '''
    
    # If O is a single point, make it the same size as V
    if (O.shape[0] == 1):
        XR = np.ones(V.shape[0],1) * O[0]
        YR = np.ones(V.shape[0],1) * O[1]
        ZR = np.ones(V.shape[0],1) * O[2]
    else:
        XR = O[:,0]; YR = O[:,1]; ZR = O[:,2]
    

    # Read rest of input
    vx = V[:,0]; vy = V[:,1]; vz = V[:,2]
    svx = SV[:,0]; svy = SV[:,1]; svz = SV[:,2]
    cxy = COR[:,0]; cxz = COR[:,1]; cyz = COR[:,2]

    # Convert origin vector to ellipsoidal coordinates
    T = np.zeros(XR.shape[0])
    E = xyz2wgs(np.vstack((T, XR, YR, ZR)).T)

    # Compute sines and cosines
    cp = np.cos(E[:,1]*np.pi/180); sp = np.sin(E[:,1]*np.pi/180) # longitude
    cl = np.cos(E[:,2]*np.pi/180); sl = np.sin(E[:,2]*np.pi/180) # latitude
    
    # For each site
    for i in range(V.shape[0]):
        # Build the rotation matrix
        R = np.vstack((
            np.hstack(( -sl[i]*cp[i],   -sl[i]*sp[i],    cl[i])),
            np.hstack((-sp[i],            cp[i],           0.0)),
            np.hstack((cl[i]*cp[i],    cl[i]*sp[i],    sl[i]))))
        
        # Apply the rotation
        NEUi = np.dot(R, np.array([vx[i], vy[i], vz[i]]))
        
        # Build covariance for that site
        CVi = np.array([[svx[i]**2, cxy[i]*svx[i]*svy[i],   cxz[i]*svx[i]*svz[i]],
               [cxy[i]*svx[i]*svy[i],   svy[i]**2, cyz[i]*svy[i]*svz[i]],
               [cxz[i]*svx[i]*svz[i],   cyz[i]*svy[i]*svz[i],   svz[i]**2]])
        
        # Propagate covariance
        CNEUi = np.dot(np.dot(R, CVi), R.T)

        # Increment result matrices
        if i != 0:
            NEU = np.vstack((NEU,
                   np.array([NEUi[0], NEUi[1], NEUi[2]])))
            CNEU = np.vstack((CNEU,
                   np.array([CNEUi[0,0], CNEUi[0,1], CNEUi[0,2], CNEUi[1,1], CNEUi[1,2], CNEUi[2,2]])))
        else:
            NEU = np.array([NEUi[0], NEUi[1], NEUi[2]])
            CNEU = np.array([CNEUi[0,0], CNEUi[0,1], CNEUi[0,2], CNEUi[1,1], CNEUi[1,2], CNEUi[2,2]])

    
    return NEU, CNEU, E


#-------------------------------------------------------------------------------

def get_wrms(v, sv):
    '''
    Calculate weighted mean and scatter about weighted mean.
    
    Call: wmean, wrms = get_wrms(v, sv)
    '''
    
    n = v.size

    # Calculate weights
    w = 1/(sv**2)

    # Calculate weighted mean
    num = np.sum(v*w)
    sumw = np.sum(w)
    wmean = num/sumw

    # Calculate scatter w.r.t. weighted mean
    chi2 = np.sum( ((v-wmean)/sv)**2 )
    wrms = np.sqrt( (n/(n-1)) * chi2/sumw )
    
    return wmean, wrms


#-------------------------------------------------------------------------------

def rot2eul(OM, COM):
    '''
    Convert rotation vector into Euler parameters.

    Input:
    - OM = rotation vector in deg/My
    - COM = associated covariance matrix in (deg/My)^2
            [Cxx Cxy Cxz
             Cxy Cyy Cyz
             Cxz Cyz Czz]

    Output:
    - E = lon, lat, angular velocity (deg/Myr)
    - CE = associated covariance matrix
    - EL = [smaj,smin,az,sigmaz]

    COM can be [] (empty matrix) if no covariance info is
    available. In that case, CE and EL will also be [].

    Call: Eu, CE, EL = rot2eul(OM, COM)
    '''
    
    # Convert everything to radians
    OM = OM * (np.pi/180)
    COM = COM * (np.pi/180)**2

    x = OM[0]; y = OM[1]; z = OM[2]

    # Convert rotation vector from xyz to neu
    ang = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))
    lon = np.arctan2(y, x)
    Eu = np.array([lon*180/np.pi, lat*180/np.pi, ang*180/np.pi])

    # Convert covariance matrix from xyz to nea
    if (COM.size == 0):
        CE = np.array([])
        EL = np.array([])
    else:
        # Compute Jacobian
        # --> d(lat)/d(x,y,z)
        jcom = (x**2 + y**2) / (x**2 + y**2 + z**2)
        j11 = -x*z / (x**2 + y**2)**1.5;   j11 = j11*jcom
        j12 = -y*z / (x**2 + y**2)**1.5;   j12 = j12*jcom
        j13 = 1 / (x**2 + y**2)**0.5;      j13 = j13*jcom
        # --> d(lon)/d(x,y,z)
        jcom = x**2 / (x**2 + y**2)
        j21 = -y / x**2;  j21 = j21*jcom
        j22 =  1 / x;    j22 = j22*jcom
        j23 = 0;
        # --> d(ang)/d(x,y,z)
        jcom = 1 / (2*np.sqrt(x**2 + y**2 + z**2))
        j31 = 2*x;  j31 = j31*jcom
        j32 = 2*y;  j31 = j31*jcom
        j33 = 2*z;  j31 = j31*jcom

        J = np.vstack((np.hstack((j11, j12, j13)),
             np.hstack((j21, j22, j23)),
             np.hstack((j31, j32, j33))))

        # Propagate covariance
        CE = np.dot(np.dot(J,COM),J.T)

        # Extract lat-lon part
        CLL = CE[0:2,0:2] * (180/np.pi)**2
        
        # Compute error ellipse, 1 sigma => p=0.39346 (in 2D)
        smaj, smin, az = errell(CLL, 0.39346)

        # Compute uncertainty on angular velocity
        sigmaz = np.sqrt(CE[2,2])

        # Output
        EL = np.array([smaj, smin, az, sigmaz])
    
    
    return Eu, CE, EL


#-------------------------------------------------------------------------------

def rotate(OM, COM, P):
    '''
    Computes velocity at selected point given a rotation vector

    Input:
        - OM(Rx,Ry,Rz) = rotation vector, deg/my
        - COM = covariance matrix of rotation vector (3x3 matrix) 
        - P = site coordinates, XYZ cartesian frame, meters
          P can be a vector or a n x 3 matrix with n site positions

    Output:
      - Vxyz = result matrix in Cartesian frame. One line per site,
        with following format:
        vx vy vz  cxx cxy cxz cyy cyz czz
                    (m/yr)    (covariance m^2/yr^2)
      - Venu = result matrix in local frame. One line per site,
        with following format:
        ve vn sve svn cor
                       (mm/yr) (unitless)

    Call: Vxyz, Venu = rotate(OM, COM, P)
    '''

    # Convert OM and COM from deg/My to rad/yr
    Re = 6378137.0
    con = (Re * np.pi/180) * 1e-6
    OM = OM * con
    COM = COM * (con**2)

    # Compute unit position vector
    X = P[:,0]; Y = P[:,1]; Z = P[:,2]
    R = np.sqrt(np.power(X,2) + np.power(Y,2) + np.power(Z,2))
    Xu = X/R; Yu = Y/R; Zu = Z/R
    
    # Compute position in ellipsoidal coordinates
    T = np.zeros(X.shape[0])
    E = xyz2wgs(np.vstack((T, X, Y, Z)).T)
    lon = E[:,1]*np.pi/180
    lat = E[:,2]*np.pi/180
    cl = np.cos(lat); sl = np.sin(lat)
    cp = np.cos(lon); sp = np.sin(lon)
    
    # Build cross product transformation matrix A
    A = np.matrix([]);
        
    for i in range(Xu.size):
        if i == 0:
            A = np.vstack((np.array([   float(0),   Zu[i], -Yu[i]]),
               np.array([-Zu[i],    float(0),   Xu[i]]),
               np.array([ Yu[i], -Xu[i],    float(0)])))
        else:
            A = np.vstack((A,
               np.array([  float(0),   Zu[i], -Yu[i]]),
               np.array([-Zu[i],    float(0),   Xu[i]]),
               np.array([ Yu[i], -Xu[i],    float(0)])))
    
    # Compute cross product
    #v = np.matrix(A) * np.matrix(OM).T
    v = np.dot(A,OM)
    # Propagate rotation vector covariance
    #CV = A * COM * A.T
    CV = np.dot(np.dot(A,COM),A.T)

    # For each site
    Vxyz = np.matrix([])
    Venu = np.matrix([])
    n = 0
    
    for i in np.arange(0,v.size,3):

        # Extract velocities in cartesian frame
        vx = v[i]; vy = v[i+1]; vz = v[i+2]

        # Extract corresponding covariance
        Cxyz = CV[i:i+3,i:i+3]

        # Build rotation matrix from cartesian to local coordinates
        R = np.vstack((np.array([ -sl[n]*cp[n],   -sl[n]*sp[n],    cl[n]]).reshape((3)), #reshape force dimensions to be the same
                       np.array([ -sp[n],            cp[n],     np.array([0])]).reshape((3)),
                       np.array([-cl[n]*cp[n],   -cl[n]*sp[n],   -sl[n]]).reshape((3))))

        # Rotate velocities from cartesian to local frame, convert to mm/yr
        tmp = np.dot(R,np.array([vx, vy, vz]))
        tmp = tmp * 1e3
        vn = tmp[0]; ve = tmp[1]; vu = tmp[2]
        
        # Rotate covariance from cartesian to local frame
        Cneu = np.dot(np.dot(R,Cxyz),R.T)

        # Extract stdev, correlation, convert to mm/yr
        svn = np.sqrt(Cneu[0,0]) * 1e3
        sve = np.sqrt(Cneu[1,1]) * 1e3
        svu = np.sqrt(Cneu[2,2]) * 1e3
        cor_ne = Cneu[0,1]/(Cneu[0,0]**0.5 * Cneu[1,1]**0.5)
        cor_nu = Cneu[0,2]/(Cneu[0,0]**0.5 * Cneu[2,2]**0.5)
        cor_eu = Cneu[1,2]/(Cneu[1,1]**0.5 * Cneu[2,2]**0.5) # correlation, unitless
        
        # Increment result matrices
        if i == 0:
            Vxyz = np.hstack((vx, vy, vz, Cxyz[0,0], Cxyz[0,1], Cxyz[0,2], Cxyz[1,1], Cxyz[1,2], Cxyz[2,2]))
            Venu = np.hstack((ve, vn, sve, svn, cor_ne))
        else:
            Vxyz = np.vstack((Vxyz, np.hstack((vx, vy, vz, Cxyz[0,0], Cxyz[0,1], Cxyz[0,2], Cxyz[1,1], Cxyz[1,2], Cxyz[2,2]))))
            Venu = np.vstack((Venu, np.hstack((ve, vn, sve, svn, cor_ne))))

        # Increment site number
        n+=1
    
    return Vxyz, Venu


#-------------------------------------------------------------------------------

def errell(C, p):
    '''
    ERRELL	Compute semi-major axis, semi-minor axis, and orientatiion
    of an error ellipse from a 2 x 2 covariance matrix

    Input:
    - C = 2 x 2 covariance matrix (deg^2)
    - p = confidence level (ex: 0.95)
                     for one-sigma (in 2D) use p=0.39346;

    Output:
    - smaj = semi-major axis (deg)
    - smin = semi-minor axis (deg)
    - azim = azimuth of semi-major axis
    (degrees clockwise from north)
    '''
    
    # in case matrix larger than 2x2
    C = C[0:2,0:2]

    # compute eigenvectors and eigenvalues
    eigval, eigvec = np.linalg.eig(C)

    # find maj and min
    smaj = np.amax(eigval)
    smaj_i = np.argmax(eigval)
    smajvec = eigvec[:,smaj_i]
    smin = np.amin(eigval)

    # compute azimuth of semi-major axis
    azim = np.arctan2(smajvec[1],smajvec[0]) * (180/np.pi)
    azim = 90 - azim

    # compute and apply chi2 (dof=2)
    chi2 = st.chi2.ppf(p,2)
    smaj = np.sqrt(chi2*smaj)
    smin = np.sqrt(chi2*smin)
    
    return smaj, smin, azim

 