#!/usr/bin/env python

# Object passed to construct an NEOvisualizer object must have the following attributes: a (AU), e (1), i (deg), om (deg), w (deg), ma (deg), n (deg/JDTDB), epoch (JDTDB)

# Use matplotlib widget magic to make plots interactive on jupyter notebook

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import secrets

class NEOvisualizer:
    """
    NEOvisualizer class holds orbital elements. NEOvisulizer object takes in matplotlib axes and plots its approximate orbit and/or position at a given time.
        
    Attributes:
    a : float
        Semimajor Axis (AU)
    e : float
        Eccentricity (dimensionless)
    i : float
        Inclination with respect to J2000 ecliptic (deg)
    om : float
        Longitude of the Ascending Node with respect to J2000 ecliptic (deg)
    w : float
        Argument of perihelion (periapsis) with respect to J2000 ecliptic (deg)
    ma : float
        Mean Anomaly at Epoch (deg)
    n : float
        Mean Orbital Motion (deg/d)
    epoch : float
        Reference Epoch (in TDB scale Julian Days)
    init_time : astropy.time Time object
        Time at which the object is instantiated.
    time : float
        init_time value in TDB scale Julian Days.
    
    """
    
    def __init__(self, obj):
        """
        Parameters
        ----------
        obj : Pandas.Series
            Series object containing columns a, e, i, om, w, ma, n, epoch (in AU, deg, JD TDB)
        """
        
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
        """
        Returns x, y, z components (AU) of one cycle of the orbit that can be plotted.
        """
        
        # return x,y,z for all pos corresponding to each step in true anomaly from -pi to pi or 0 to 2 pi
        ecc = np.linspace(-np.pi, np.pi, 500)
        trueAnom = self.ecc2trueAnom(ecc)
        x,y,z = self.trueAnom2pos(trueAnom)
        return x,y,z
    
    def getPositionAt(self, time=None):
        """
        Returns x, y, z components (AU) of the NEO's position at either time of instance or user specified time.
        
        Parameters
        ----------
        time (Optional) : float
            User specified time for the object's position
            defaults to current time if not specified.
            
        Returns
        -------
        x : float
            x component of position at given time (AU)
        y : float
            y component of position at given time (AU)
        z : float
            z component of position at given time (AU)
        """
        
        # return x,y,z position at either init time or user specified time
        # time must be in JD TDB
        if time is None:
            days_since_epoch = np.ceil(self.time - self.epoch)
            timeline = np.linspace(self.epoch, self.time, int(days_since_epoch))
        else: # if time is specified
            days_since_epoch = np.ceil(time - self.epoch)
            timeline = np.linspace(self.epoch, time, int(days_since_epoch))
        trueAnom = self.ecc2trueAnom(self.getEccAnom(timeline))
        x,y,z = self.trueAnom2pos(trueAnom)
        return x[-1], y[-1], z[-1] #current position
    
    def plotOrbit(self, ax, showNEO=True):
        """
        Takes in matplotlib 3d subplot axes onto which the orbit of the NEO is plotted. Also plots the current position of the NEO unless otherwise specified.
        
        Parameters
        ----------
        ax : matplotlib.axes._subplots.Axes3DSubplot
            Axes from a 3d subplot for orbit and position to be plotted on.
        showNEO (Optional) : bool
            Whether the current position of the object to be displayed. Defaults to True unless otherwise specified.
        """
        
        # plot orbit
        x,y,z = self.getOrbit()
        ax.plot(x,y,z, color='#ced3db', linestyle='-', alpha=0.5)
        if showNEO:
            # plot NEO in today's position
            x,y,z = self.getPositionAt(self.time)
            colors = ['xkcd:aquamarine', 'xkcd:pink', 'xkcd:goldenrod', 'xkcd:lavender', 'xkcd:khaki']
            ax.scatter(x,y,z, color=secrets.choice(colors), s=20)
            
        
    def plotNEO(self, ax, colorStr=None, time=None):
        """
        Takes in matplotlib 3d subplot axes onto which the position of the NEO is plotted. Supports user specified (1) color for NEO, and (2) particular time of interest for position.
        
        Parameters
        ----------
        ax : matplotlib.axes._subplots.Axes3DSubplot
            Axes from a 3d subplot for the position to be plotted on.
        colorStr (Optional) : string
            matplotlib readable color string. Defaults to colors['xkcd:aquamarine', 'xkcd:pink', 'xkcd:goldenrod', 'xkcd:lavender', 'xkcd:khaki'] unless otherwise specified
        time (Optional) : float
            User specified time of interest. Defaults to current time unless otherwise specified.
        """
        
        # Time condition
        if time is None:
            x,y,z = self.getPositionAt(self.time)
        else:
            x,y,z = self.getPositionAt(time)
        
        # Color condition
        if colorStr is None:
            colors = ['xkcd:aquamarine', 'xkcd:pink', 'xkcd:goldenrod', 'xkcd:lavender', 'xkcd:khaki']
            ax.scatter(x,y,z, color=secrets.choice(colors), s=20)
        else:
            ax.scatter(x,y,z, color=colorStr, s=20)
        
            
    def getEccAnom(self, timeline, tol=0.00000001):
        """
        Returns eccentric anomaly:
        Computes mean anomaly from mean anomaly at epoch and epoch, then computes eccentric anomaly from mean anomaly using numerical root finding method
        
        Parameters
        ----------
        timeline : array-like
            Considered time (in JD TDB) to be converted to the corresponding mean anomaly.
            
        Returns
        -------
        E : array-like
            Eccentric anomaly. Same data type and dimensions as the input.
        """
        
        # converts mean anomaly to eccentric anomaly using Newton-Rhapson root finding method
        M = self.ma + self.n*(timeline - self.epoch)
        E = M # Eccentric anomaly, initial value set equal to mean anomaly
        err = 1;
        while np.all(err >= tol):
            temp = E
            E = E - ((E-self.e*np.sin(E)-M)/(1-self.e*np.cos(E)))
            err = abs(temp - E)
        return E
    
    def ecc2trueAnom(self, eccAnom):
        """
        Converts eccentric anomaly to true anomaly.
        
        Parameters
        ----------
        eccAnom : array-like
            Eccentric anomaly. 
            
        Returns
        -------
        trueAnom : array-like
            True anomaly. Same data type and dimensions as the input.
        """
        # converts eccentric anomaly to true anomaly using arctan2 (2 arguments)
        trueAnom = 2*np.arctan2(np.sqrt(1+self.e)*np.sin(eccAnom/2), np.sqrt(1-self.e)*np.cos(eccAnom/2))
        return trueAnom
    
    def trueAnom2pos(self, trueAnom):
        """
        Returns x, y, z components of positions (in J2000 ecliptic reference frame) corresponding to the input true anomaly.
        
        Parameters
        ----------
        trueAnom : array-like
            True anomaly of interest.
        
        Returns
        -------
        X,Y,Z : array-like
            Cartesian position vector components in J2000 ecliptic reference frame. Each has the same datatype and dimensions as the input.
        """
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
        
############## Non-class Function ##############
        
def plotEarth(ax, showOrbit=True):
    """
    Takes in matplotlib 3d subplot axes onto which the Earth is plotted. Also plots the orbit of the Earth unless otherwised specified.
    
    Parameters
    ----------
    ax : matplotlib.axes._subplots.Axes3DSubplot
        Axes from a 3d subplot for the position to be plotted on.
    showOrbit (Optional) : bool
        whether the orbit of the Earth to be displayed. Defaults to True unless otherwise specified.
    """
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
    ser = pd.Series(elems) # construct dummy object to pass in NEOvis constructor
    earth = NEOvisualizer(ser)
    
    px,py,pz = earth.getPositionAt(earth.time)
    ax.scatter(px,py,pz, c='green', s=27)
        
    if showOrbit:
        # get orbit
        ecc = np.linspace(-np.pi, np.pi, 500)
        trueAnom = earth.ecc2trueAnom(ecc)
        x,y,z = earth.trueAnom2pos(trueAnom)
        # plot orbit
        ax.plot(x,y,z, color='#89ad93', linestyle='-')