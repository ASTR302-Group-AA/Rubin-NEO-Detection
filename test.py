#!/usr/bin/env python
import warnings
warnings.filterwarnings('ignore')
import orbit_detection as od
filename = 'testing.csv'
obj_list = ['91mins', 'not3days', '15days', 'all']
expected_output = [False, False, False, True]
for i in range(len(obj_list)):
    print(od.detector(filename, obj_list[i]))
    print(expected_output[i])