#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.table import QTable
from astropy.time import Time

def detector(filename, objID):
    '''
    Takes sorcha data and desired object and finds out if it would be confirmed as detected.
    Inputs:
        filename - the name of the file with the observation data
        objID - the ID of the desired object
    Outputs:
        output - boolean value on if the object was detected
        output_text - message including the boolean value and the dates of the second observation in a detection pair
        
    Doesn't sort by date yet, need to create test suite.
        '''
    
    data = QTable.read(filename)
    
    data = data[[data['ObjID'] == objID]]
    
    data.sort(['FieldMJD_TAI'])
    
    TIME = Time(data['FieldMJD_TAI'], format = 'mjd')
    # obj = data['ObjID'][data['ObjID'] == objID]
    
    detections = []
    for i in range(len(data) - 1):
        deltaT = (TIME[i+1] - TIME[i]).to_value("sec") / 60 #converting to minutes between observations
        if deltaT < 90:
            detections.append(TIME[i])
    
    detected = []
    for j in range(len(detections) - 1):
        deltaT = (detections[j+1] - detections[j]).to_value("sec") / 604800 #converting to weeks between detections
        if deltaT < 2:
            detected.append(detections[j+1].to_value("iso"))
    
    if len(detected) >= 3:
        output_text = f'Object {objID} was successfully detected. Paired observations occured: {detected}.'
        output = True
    else:
        output_text = f'Object {objID} was not detected.'
        output = False
    
    return output, output_text