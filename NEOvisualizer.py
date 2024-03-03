#!/usr/bin/env python

# Object passed to construct an NEOvisualizer object must have the following attributes: a (AU), e (1), i (deg), om (deg), w (deg), ma (deg), n (deg/JDTDB), epoch (JDTDB)

# Use matplotlib widget magic to make plots interactive on jupyter notebook

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from mpl_toolkits.mplot3d import axes3d
import pandas as pd

class NEOvisualizer:
    """
        docstring goes here 
        
    Attributes:
    
    """
    def __init__(self, obj):
        # self.elem = elem
        self.a = obj.a
        self.e = obj.e
        self.i = np.radians(obj.i)
        self.om = np.radians(obj.om)
        self.w = np.radians(obj.w)
        self.ma = np.radians(obj.ma)
        self.n = np.radians(obj.n) # mean motion in rad/d
        self.epoch = Time(obj.epoch, format='jd', scale='tdb').jd # dtype:float. epoch must be in JD (TDB)
        self.init_time = Time.now() # dtype: Time
        self.time = self.init_time.tdb.jd # dtype: float
        
        
    def getOrbit(self):
        # return x,y,z for all pos corresponding to each step in true anomaly from -pi to pi or 0 to 2 pi or whatever
        ecc = np.linspace(-np.pi, np.pi, 500)
        trueAnom = self.ecc2trueAnom(ecc)
        x,y,z = self.trueAnom2pos(trueAnom)
        return x,y,z
    
    def getPositionAt(self, time):
        # return current (at init_time) x,y,z position
        days_since_epoch = np.ceil(self.time - self.epoch)
        timeline = np.linspace(self.epoch, self.time, int(days_since_epoch))
        trueAnom = self.ecc2trueAnom(self.getEccAnom(timeline))
        x,y,z = self.trueAnom2pos(trueAnom)
        return x[-1], y[-1], z[-1] #current position
    
    def plotOrbit(self, ax, showNEO=True):
        # plot orbit
        x,y,z = self.getOrbit()
        ax.plot(x,y,z, color='#ced3db', linestyle='-')
        if showNEO:
            # plot NEO in today's position
            x,y,z = self.getPositionAt(self.init_time)
            ax.scatter(x,y,z, color='cyan', s=20)
            
        
    def plotNEO(self, ax):
        # plot NEO in today's position
        x,y,z = self.getPositionAt(self.init_time)
        ax.scatter(x,y,z, color='cyan', s=20)
        
            
    def getEccAnom(self, time, tol=0.000000001):
        # converts mean anomaly to eccentric anomaly using Newton-Rhapson root finding method
        M = self.ma + self.n*(time - self.epoch)
        E = M # Eccentric anomaly, initial value set equal to mean anomaly
        err = 1;
        while np.all(err >= tol):
            temp = E
            E = E - ((E-self.e*np.sin(E)-M)/(1-self.e*np.cos(E)))
            err = abs(temp - E)
        return E
    
    def ecc2trueAnom(self, eccAnom):
        # converts eccentric anomaly to true anomaly using arctan2 (2 arguments)
        trueAnom = 2*np.arctan2(np.sqrt(1+self.e)*np.sin(eccAnom/2), np.sqrt(1-self.e)*np.cos(eccAnom/2))
        return trueAnom
    
    def trueAnom2pos(self, trueAnom):
        # true anomaly to radius
        r = self.a*((1-self.e**2)/(1+self.e*np.cos(trueAnom)))
        # pos in perifocal frame (heliocentric) z comp is 0.
        x = r*np.cos(trueAnom)
        y = r*np.sin(trueAnom) 
        # pos in J2000 ecliptic reference frame (heliocentric)
        X= x*(np.cos(self.om)*np.cos(self.w)-np.sin(self.om)*np.cos(self.i)*np.sin(self.w)) - y*(np.cos(self.om)*np.sin(self.w)+np.sin(self.om)*np.cos(self.i)*np.cos(self.w))
        Y= x*(np.sin(self.om)*np.cos(self.w)+np.cos(self.om)*np.cos(self.i)*np.sin(self.w)) - y*(np.sin(self.om)*np.sin(self.w)-np.cos(self.om)*np.cos(self.i)*np.cos(self.w))
        Z= x*(np.sin(self.i)*np.sin(self.w)) + y*(np.sin(self.i)*np.cos(self.w))
        return X,Y,Z
        

        
def plotEarth(ax, showOrbit=True):
    # Orbital Elements
    a = 1.00000011 # AU
    e = 0.01671022
    i = 0
    om = np.radians(-11.26064)
    w = np.radians(102.94719)
    la = np.radians(100.46435)
    ma = la - om - w
    epoch = Time(2451545.0, format='jd', scale='tt').tdb
    n = 0.98560912 # deg/d
    
    elems = {'a':a, 'e':e, 'i':i, 'om':om, 'w':w, 'ma':ma, 'epoch':epoch, 'n':n}
    ser = pd.Series(elems)
    earth = NEOvisualizer(ser)
    
    px,py,pz = earth.getPositionAt(earth.time)
    ax.scatter(px,py,pz, c='green', s=25)
        
    if showOrbit:
        # get orbit
        ecc = np.linspace(-np.pi, np.pi, 500)
        trueAnom = earth.ecc2trueAnom(ecc)
        x,y,z = earth.trueAnom2pos(trueAnom)
        # plot orbit
        ax.plot(x,y,z, color='#89ad93', linestyle='-')