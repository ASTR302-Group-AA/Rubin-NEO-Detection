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
        '''
    
    data = QTable.read(filename)
    
    #names = data['ObjID']
    
    right_name = []
    for i in range(len(data)):
        if str(data['ObjID'][i]) == objID:
            right_name.append(True)
            
        else:
            right_name.append(False)
    
    return right_name
    
    #data = data[data['ObjID'] == objID]
    
    data.sort(['FieldMJD_TAI'])
    
    TIME = Time(data['FieldMJD_TAI'], format = 'mjd') - 0.67 #changing all times to noon local time
    output = False
    
    pairs = []
    for i in range(len(data) - 1):
        deltaT = (TIME[i+1] - TIME[i])
        if deltaT < .0625: #90 minutes = .0625 days
            pairs.append(TIME[i])
    if len(pairs) < 2:
        output_text = f'Object {objID} was not detected. Error: No observation pairs.'
        return output, output_text
    
    days = []

    for k in range(len(pairs)):
        days.append(pairs[k])
        if len(days)>1:
            for l in range(len(days) - 1):
                if int(pairs[k].value) - int(days[l].value) == 0:
                    days.pop(-1)
    if len(days) < 3:
        output_text = f'Object {objID} was not detected. Error: Not 3 days.'
        return output, output_text
    
    detected = False
    j = 0
    while np.logical_and(j < len(days) - 2, detected == False):
        if (int(days[j+2].value) - int(days[j].value)) <= 14:
            detected = True
        j+=1
        
    if (detected):
        output_text = f'Object {objID} was successfully detected. Third pair of observations occured: {days[j+1].to_value("iso")}.'
        output = True
    else:
        output_text = f'Object {objID} was not detected. Error: More than 2 weeks.'
        output = False
    
    return output, output_text