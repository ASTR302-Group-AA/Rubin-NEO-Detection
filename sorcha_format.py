#!/usr/bin/env python

import pandas as pd

def sorcha_format(filename):
    
    sorcha_data = pd.read_csv(filename)
    
    renamed_columns = sorcha_data.rename(columns={"pdes":"ObjID",
                                                  "epoch_mjd":"epochMJD_TDB", "i":"inc",
                                                  "om":"node", "w":"argPeri", "H":"H_r"})
    
    orbit_file = renamed_columns.iloc[:, 2:8]
    orbit_file.insert(0, 'ObjID', renamed_columns['ObjID'], allow_duplicates = False)
    orbit_file.insert(1, 'FORMAT', 'KEP', allow_duplicates = False)
    
    #take the following generalized values for an asteroid's color profile: 
    #u-r g-r i-r z-r y-r GS
    #1.72 0.48 -0.11 -0.12 -0.12 0.15
    
    pparams_file = renamed_columns.iloc[:,0:2]
    pparams_file.insert(2, 'u-r', 1.72)
    pparams_file.insert(3, 'g-r', 0.48)
    pparams_file.insert(4, 'i-r', -0.11)
    pparams_file.insert(5, 'z-r', -0.12)
    pparams_file.insert(6, 'y-r', -0.12)
    pparams_file.insert(7, 'GS', 0.15)   
       
    return orbit_file, pparams_file
       
orbit_file, pparams_file = sorcha_format('test_data_final.csv')

orbit_file.to_csv('orbit_file.des',index=False, sep=',')
pparams_file.to_csv('pparams_file.txt',index=False, sep=',')

#------------------------------------------------------------------------------------

#import pandas as pd

# def sorcha_format(filename):
    
#     sorcha_data = pd.read_csv(filename)
    
    
    
#     renamed_columns = sorcha_data.rename(columns={"pdes":"ObjID", "epoch_mjd":"epochMJD_TDB", "i":"inc", "om":"node", "w":"argPeri", "H":"H_r"})
    
#     # u-r g-r i-r z-r y-r GS
#     #1.72 0.48 -0.11 -0.12 -0.12 0.15
    
#     orbit_file = renamed_columns.iloc[:,0:[2,7]]
        
#     sorcha_data.insert(1, 'FORMAT', 'KEP', allow_duplicates = False)
    
#     sorcha_renamed.insert(2, 'u-r', 1.72)
#     sorcha_renamed.insert(3, 'g-r', 0.48)
#     sorcha_renamed.insert(4, 'i-r', -0.11)
#     sorcha_renamed.insert(5, 'z-r', -0.12)
#     sorcha_renamed.insert(6, 'y-r', -0.12)
#     sorcha_renamed.insert(7, 'GS', 0.15)
    
    
#     return sorcha_renamed
       
# sorcha_format('test_data_final.csv')
